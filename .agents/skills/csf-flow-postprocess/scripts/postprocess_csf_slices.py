#!/usr/bin/env python3
"""Extract CSF slice flow, area, pressure, and tissue-displacement metrics."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

from slice_common import (
    FenicsSeriesReader,
    discover_xdmf_series,
    ensure_output_dir,
    flow_slice_grid,
    integrate_surface,
    max_vector_magnitude_on_surface,
    mesh_grid_from_h5,
    numeric_h5_groups,
    parse_chunk,
    parse_indices,
    parse_xdmf,
    read_xdmf_grid,
    result_index,
    slice_at_z,
    source_signature,
    threshold_domain,
    warp_if_available,
)

METRICS = ("flow-rate", "area-deformation", "pressure", "max-displacement")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--results-dir", type=Path, help="Folder containing result*.xdmf")
    src.add_argument("--snapshots-h5", type=Path, help="Compact FEniCS function snapshots")
    p.add_argument("--mesh-h5", type=Path, help="FEniCS mesh HDF5 (required with snapshots)")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--metrics", nargs="+", choices=METRICS, default=list(METRICS))
    p.add_argument("--first", type=int)
    p.add_argument("--last", type=int)
    p.add_argument("--stride", type=int, default=1)
    p.add_argument("--indices", help="HDF5 indices, e.g. 0,3,6 or 0:100:3")
    p.add_argument("--snapshot-chunk", help="Collision-safe serial chunk K/N")
    p.add_argument("--dt-s", type=float)
    p.add_argument("--period-s", type=float)
    p.add_argument("--slice-spacing-mm", type=float, default=1.0)
    p.add_argument("--fm-z-mm", type=float, help="Override first-snapshot fluid z_top")
    p.add_argument("--max-slices", type=int)
    p.add_argument("--velocity-name", default="Velocity")
    p.add_argument("--pressure-name", default="Pressure")
    p.add_argument("--displacement-name", default="Displacement")
    p.add_argument("--domain-name", default="Sub-domain")
    p.add_argument("--velocity-group", default="velocity")
    p.add_argument("--pressure-group", default="pressure")
    p.add_argument("--displacement-group", default="displacement")
    p.add_argument("--fluid-tag", type=float)
    p.add_argument("--solid-tag", type=float)
    p.add_argument("--pressure-scale", type=float, default=1.0)
    p.add_argument("--pressure-unit", default="solver_units")
    p.add_argument("--waveform-distances-mm", help="Comma-separated distances")
    p.add_argument("--profile-times-s", help="Comma-separated physical times")
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def mpi_context() -> tuple[object | None, int, int]:
    try:
        from mpi4py import MPI

        return MPI.COMM_WORLD, MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size
    except ImportError:
        return None, 0, 1


def require_fields(grid, names: list[str]) -> None:
    missing = [name for name in names if name not in grid.point_data]
    if missing:
        raise KeyError(f"Missing required point field(s): {', '.join(missing)}")


def source_setup(args: argparse.Namespace):
    """Return descriptors plus closures that load a full snapshot grid."""
    if args.results_dir:
        paths = discover_xdmf_series(args.results_dir, args.first, args.last, args.stride)
        if not paths:
            raise FileNotFoundError(f"No selected result*.xdmf files in {args.results_dir}")
        descriptors = []
        for path in paths:
            meta = parse_xdmf(path)
            if meta["time_s"] is None:
                raise ValueError(f"Embedded physical time is missing in {path}")
            descriptors.append((result_index(path), float(meta["time_s"]), path))

        def load(desc):
            return read_xdmf_grid(desc[2])

        return "xdmf-series", descriptors, load, paths

    if args.mesh_h5 is None:
        raise ValueError("--mesh-h5 is required with --snapshots-h5")
    if args.dt_s is None and args.period_s is None:
        raise ValueError("HDF5 processing requires --dt-s or --period-s")
    available = numeric_h5_groups(args.snapshots_h5, args.velocity_group)
    if not available:
        raise KeyError(f"No numeric snapshots in /{args.velocity_group}")
    indices = parse_indices(args.indices, available)
    indices = [i for i in indices if (args.first is None or i >= args.first) and (args.last is None or i <= args.last)]
    if args.stride > 1 and indices:
        start = indices[0]
        indices = [i for i in indices if (i - start) % args.stride == 0]
    if not indices:
        raise ValueError("No HDF5 snapshots selected")
    if args.dt_s is not None:
        times = [(i - indices[0]) * args.dt_s for i in indices]
    else:
        times = np.linspace(0.0, args.period_s, len(indices), endpoint=False).tolist()
    base, _, _, domains = mesh_grid_from_h5(args.mesh_h5)
    reader = FenicsSeriesReader(args.mesh_h5, args.snapshots_h5)
    velocity_set = set(numeric_h5_groups(args.snapshots_h5, args.velocity_group))
    pressure_set = set(numeric_h5_groups(args.snapshots_h5, args.pressure_group))
    displacement_set = set(numeric_h5_groups(args.snapshots_h5, args.displacement_group))
    descriptors = [(i, float(t), None) for i, t in zip(indices, times)]

    def load(desc):
        index = desc[0]
        grid = base.copy(deep=True)
        if index not in velocity_set:
            raise KeyError(f"Missing /{args.velocity_group}/{index}")
        grid.point_data[args.velocity_name] = reader.read_vector(f"/{args.velocity_group}/{index}")
        if index in pressure_set:
            grid.point_data[args.pressure_name] = reader.read_scalar(f"/{args.pressure_group}/{index}")
        if index in displacement_set:
            grid.point_data[args.displacement_name] = reader.read_vector(f"/{args.displacement_group}/{index}")
        if domains is None:
            grid.cell_data.pop(args.domain_name, None)
        elif args.domain_name != "Sub-domain":
            grid.cell_data[args.domain_name] = grid.cell_data.pop("Sub-domain")
        return grid

    return "fenics-h5", descriptors, load, [args.mesh_h5, args.snapshots_h5]


def domain_parts(grid, args, reader_kind: str):
    has_domain = args.domain_name in grid.cell_data or args.domain_name in grid.point_data
    fluid_default = 0.0 if reader_kind == "xdmf-series" else 258514.0
    solid_default = 1.0 if reader_kind == "xdmf-series" else 258515.0
    if has_domain:
        fluid = threshold_domain(grid, args.domain_name, args.fluid_tag if args.fluid_tag is not None else fluid_default, required=True)
        solid = threshold_domain(grid, args.domain_name, args.solid_tag if args.solid_tag is not None else solid_default, required=False)
    else:
        fluid, solid = grid, None
    return fluid, solid, has_domain


def fixed_planes(args, reader_kind, descriptors, load):
    first_grid = load(descriptors[0])
    fluid, _, _ = domain_parts(first_grid, args, reader_kind)
    physical, warped = warp_if_available(fluid, args.displacement_name)
    z_top = float(args.fm_z_mm if args.fm_z_mm is not None else physical.bounds[5])
    z_bottom = float(physical.bounds[4])
    distance, edges, planes = flow_slice_grid(z_top, z_bottom, args.slice_spacing_mm, args.max_slices)
    return distance, edges, planes, z_top, z_bottom, warped


def process_one(desc, load, args, reader_kind, planes, distances):
    index, time_s, path = desc
    grid = load(desc)
    fluid, solid, has_domain = domain_parts(grid, args, reader_kind)
    need_velocity = "flow-rate" in args.metrics
    need_pressure = "pressure" in args.metrics
    require_fields(fluid, ([args.velocity_name] if need_velocity else []) + ([args.pressure_name] if need_pressure else []))
    fluid_physical, fluid_warped = warp_if_available(fluid, args.displacement_name)
    solid_physical = None
    if solid is not None:
        solid_physical, _ = warp_if_available(solid, args.displacement_name)
    rows = []
    for distance, z_mm in zip(distances, planes):
        fluid_slice = slice_at_z(fluid_physical, z_mm)
        integrals = integrate_surface(
            fluid_slice,
            [args.pressure_name] if need_pressure else [],
            args.velocity_name if need_velocity else None,
        )
        area = integrals["area"]
        flow = integrals.get(f"integral:{args.velocity_name}[2]", math.nan)
        p_int = integrals.get(f"integral:{args.pressure_name}", math.nan)
        p_mean = p_int / area * args.pressure_scale if area > 0 and np.isfinite(p_int) else math.nan
        max_disp = math.nan
        if "max-displacement" in args.metrics and solid_physical is not None and args.displacement_name in solid_physical.point_data:
            max_disp = max_vector_magnitude_on_surface(slice_at_z(solid_physical, z_mm), args.displacement_name)
        rows.append(
            {
                "snapshot_index": index,
                "time_s": time_s,
                "distance_from_FM_mm": float(distance),
                "z_plane_mm": float(z_mm),
                "area_mm2": area,
                "flow_rate_mm3_per_s": flow,
                "relative_area_deformation": math.nan,
                "mean_pressure": p_mean,
                "max_tissue_displacement_mm": max_disp,
                "geometry_warped": fluid_warped,
                "source_file": "" if path is None else str(path),
                "domain_field_present": has_domain,
            }
        )
    return rows


def rows_to_arrays(rows, selected, baseline_area=None):
    distances = np.array(sorted({row["distance_from_FM_mm"] for row in rows}), dtype=float)
    order = {desc[0]: i for i, desc in enumerate(selected)}
    rows.sort(key=lambda row: (order[row["snapshot_index"]], row["distance_from_FM_mm"]))
    times = np.array([desc[1] for desc in selected], dtype=float)
    shape = (len(selected), len(distances))
    arrays = {}
    columns = ["area_mm2", "flow_rate_mm3_per_s", "relative_area_deformation", "mean_pressure", "max_tissue_displacement_mm"]
    for column in columns:
        arrays[column] = np.asarray([row[column] for row in rows], dtype=float).reshape(shape)
    baseline = arrays["area_mm2"][0] if baseline_area is None else np.asarray(baseline_area, dtype=float)
    relative = np.full(shape, np.nan)
    valid = np.isfinite(baseline) & (baseline > 0)
    relative[:, valid] = (arrays["area_mm2"][:, valid] - baseline[valid]) / baseline[valid]
    arrays["relative_area_deformation"] = relative
    for i, row in enumerate(rows):
        row["relative_area_deformation"] = relative.reshape(-1)[i]
    return times, distances, arrays


def write_outputs(args, reader_kind, selected, rows, times, distances, planes, arrays, metadata):
    ensure_output_dir(args.output_dir, args.overwrite)
    fields = list(rows[0])
    with (args.output_dir / "csf_slice_metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    with (args.output_dir / "selected_snapshots.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream); writer.writerow(["snapshot_index", "time_s", "source_file"])
        writer.writerows([[i, t, "" if p is None else str(p)] for i, t, p in selected])
    np.savez_compressed(args.output_dir / "csf_slice_metrics.npz", times_s=times, distance_from_FM_mm=distances, z_plane_mm=planes, **arrays)
    (args.output_dir / "csf_slice_metrics_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    reader_kind, descriptors, load, source_paths = source_setup(args)
    comm, rank, size = mpi_context() if reader_kind == "xdmf-series" else (None, 0, 1)
    distances, edges, planes, z_top, z_bottom, reference_warped = fixed_planes(args, reader_kind, descriptors, load)
    local_positions = np.arange(len(descriptors))[rank::size]
    if args.snapshot_chunk:
        if size > 1:
            raise ValueError("Do not combine MPI and --snapshot-chunk")
        local_positions = parse_chunk(args.snapshot_chunk, len(descriptors))
    local_rows = []
    for pos in local_positions:
        local_rows.extend(process_one(descriptors[int(pos)], load, args, reader_kind, planes, distances))
    if comm is not None and size > 1:
        gathered = comm.gather(local_rows, root=0)
        if rank != 0:
            return
        local_rows = [row for part in gathered for row in part]
    selected = [descriptors[int(pos)] for pos in local_positions] if args.snapshot_chunk else descriptors
    if not local_rows:
        raise RuntimeError("Selected chunk contains no snapshots")
    baseline_area = None
    if args.snapshot_chunk and selected[0][0] != descriptors[0][0]:
        reference_rows = process_one(descriptors[0], load, args, reader_kind, planes, distances)
        baseline_area = [row["area_mm2"] for row in reference_rows]
    times, distances, arrays = rows_to_arrays(local_rows, selected, baseline_area)
    unavailable = []
    if "max-displacement" in args.metrics and np.isnan(arrays["max_tissue_displacement_mm"]).all():
        unavailable.append("max tissue displacement: no solid slice with displacement was available")
    metadata = {
        "schema_version": 1,
        "reader": reader_kind,
        "source_signature": source_signature(source_paths, {
            "reader": reader_kind, "metrics": args.metrics, "velocity_name": args.velocity_name,
            "pressure_name": args.pressure_name, "displacement_name": args.displacement_name,
            "domain_name": args.domain_name, "velocity_group": args.velocity_group,
            "pressure_group": args.pressure_group, "displacement_group": args.displacement_group,
            "fluid_tag": args.fluid_tag, "solid_tag": args.solid_tag,
            "pressure_scale": args.pressure_scale, "slice_spacing_mm": args.slice_spacing_mm,
            "fm_z_mm": args.fm_z_mm, "indices": args.indices, "first": args.first,
            "last": args.last, "stride": args.stride, "dt_s": args.dt_s, "period_s": args.period_s,
        }),
        "metrics_requested": args.metrics,
        "definitions": {
            "flow_rate_mm3_per_s": "integral_A v_z dA; positive v_z points toward increasing global z",
            "relative_area_deformation": "(A(d,t)-A(d,t0))/A(d,t0), with t0 the first selected snapshot",
            "mean_pressure": "integral_A p dA / A, multiplied by pressure_scale",
            "max_tissue_displacement_mm": "maximum displacement magnitude on the solid-domain plane slice",
        },
        "units": {"coordinates": "mm", "velocity": "mm/s", "area": "mm^2", "flow_rate": "mm^3/s", "pressure": args.pressure_unit, "displacement": "mm"},
        "pressure_scale": args.pressure_scale,
        "slice_convention": {"spacing_mm": args.slice_spacing_mm, "first_center_mm": 0.5 * args.slice_spacing_mm, "z_top_mm": z_top, "z_bottom_mm": z_bottom, "distance_formula": "z_top(first selected physical fluid geometry) - z_plane"},
        "geometry": {"reference_snapshot_warped": reference_warped, "rigid_wall_rule": "stored geometry used when displacement is absent"},
        "unavailable_metrics": unavailable,
        "snapshot_chunk": args.snapshot_chunk,
        "field_paths": {"velocity": args.velocity_name if reader_kind == "xdmf-series" else f"/{args.velocity_group}/INDEX", "pressure": args.pressure_name if reader_kind == "xdmf-series" else f"/{args.pressure_group}/INDEX", "displacement": args.displacement_name if reader_kind == "xdmf-series" else f"/{args.displacement_group}/INDEX", "domain": args.domain_name},
    }
    if args.snapshot_chunk:
        tag = args.snapshot_chunk.replace("/", "_of_")
        args.output_dir = args.output_dir / f"partial_{tag}"
    write_outputs(args, reader_kind, selected, local_rows, times, distances, planes, arrays, metadata)
    if not args.no_plots:
        from plot_csf_metrics import plot_all

        plot_all(args.output_dir, args.waveform_distances_mm, args.profile_times_s)
    print(f"Wrote {len(local_rows)} slice records to {args.output_dir}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
