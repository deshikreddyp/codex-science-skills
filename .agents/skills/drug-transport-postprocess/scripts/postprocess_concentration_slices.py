#!/usr/bin/env python3
"""Extract cross-sectional c_net and c_avg from drug-transport simulations."""

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
    integrate_surface,
    mesh_grid_from_h5,
    numeric_h5_groups,
    parse_chunk,
    parse_indices,
    parse_xdmf,
    read_xdmf_grid,
    result_index,
    source_signature,
    threshold_domain,
    transport_slice_grid,
    warp_if_available,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__); src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--results-dir", type=Path); src.add_argument("--snapshots-h5", type=Path)
    p.add_argument("--mesh-h5", type=Path); p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--first", type=int); p.add_argument("--last", type=int); p.add_argument("--stride", type=int, default=1); p.add_argument("--indices")
    p.add_argument("--snapshot-chunk", help="Collision-safe serial chunk K/N")
    p.add_argument("--dt-s", type=float); p.add_argument("--period-s", type=float)
    p.add_argument("--slice-spacing-mm", type=float, default=1.0); p.add_argument("--fm-z-mm", type=float); p.add_argument("--z-bottom-mm", type=float); p.add_argument("--max-slices", type=int)
    p.add_argument("--concentration-name", default="Concentration"); p.add_argument("--displacement-name", default="Displacement"); p.add_argument("--domain-name", default="Sub-domain")
    p.add_argument("--concentration-group", default="concentration"); p.add_argument("--displacement-group", default="displacement")
    p.add_argument("--fluid-tag", type=float); p.add_argument("--concentration-unit", default="solver_units")
    p.add_argument("--profile-times-s", help="Comma-separated physical times")
    p.add_argument("--cavg-limits", nargs=2, type=float); p.add_argument("--cnet-limits", nargs=2, type=float)
    p.add_argument("--no-plots", action="store_true"); p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def mpi_context():
    try:
        from mpi4py import MPI
        return MPI.COMM_WORLD, MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size
    except ImportError:
        return None, 0, 1


def source_setup(a):
    if a.results_dir:
        paths = discover_xdmf_series(a.results_dir, a.first, a.last, a.stride)
        if not paths: raise FileNotFoundError(f"No selected result*.xdmf files in {a.results_dir}")
        desc = []
        for path in paths:
            meta = parse_xdmf(path)
            if meta["time_s"] is None: raise ValueError(f"Embedded physical time is missing in {path}")
            desc.append((result_index(path), float(meta["time_s"]), path))
        return "xdmf-series", desc, lambda d: read_xdmf_grid(d[2]), paths
    if a.mesh_h5 is None: raise ValueError("--mesh-h5 is required with --snapshots-h5")
    if a.dt_s is None and a.period_s is None: raise ValueError("HDF5 processing requires --dt-s or --period-s")
    available = numeric_h5_groups(a.snapshots_h5, a.concentration_group)
    if not available: raise KeyError(f"No numeric snapshots in /{a.concentration_group}")
    indices = parse_indices(a.indices, available)
    indices = [i for i in indices if (a.first is None or i >= a.first) and (a.last is None or i <= a.last)]
    if a.stride > 1 and indices:
        origin = indices[0]; indices = [i for i in indices if (i - origin) % a.stride == 0]
    if not indices: raise ValueError("No HDF5 snapshots selected")
    times = [(i - indices[0]) * a.dt_s for i in indices] if a.dt_s is not None else np.linspace(0, a.period_s, len(indices), endpoint=False).tolist()
    base, domains = mesh_grid_from_h5(a.mesh_h5); reader = FenicsSeriesReader(a.mesh_h5, a.snapshots_h5)
    displacements = set(numeric_h5_groups(a.snapshots_h5, a.displacement_group)); desc = [(i, float(t), None) for i, t in zip(indices, times)]
    def load(d):
        index = d[0]; grid = base.copy(deep=True); grid.point_data[a.concentration_name] = reader.read_scalar(f"/{a.concentration_group}/{index}")
        if index in displacements: grid.point_data[a.displacement_name] = reader.read_vector(f"/{a.displacement_group}/{index}")
        if domains is None: grid.cell_data.pop(a.domain_name, None)
        elif a.domain_name != "Sub-domain": grid.cell_data[a.domain_name] = grid.cell_data.pop("Sub-domain")
        return grid
    return "fenics-h5", desc, load, [a.mesh_h5, a.snapshots_h5]


def prepare(grid, a, reader_kind):
    has_domain = a.domain_name in grid.cell_data or a.domain_name in grid.point_data
    if has_domain:
        default_tag = 0.0 if reader_kind == "xdmf-series" else 258514.0
        grid = threshold_domain(grid, a.domain_name, a.fluid_tag if a.fluid_tag is not None else default_tag, required=True)
    if a.concentration_name not in grid.point_data: raise KeyError(f"Missing point scalar {a.concentration_name!r}")
    physical, warped = warp_if_available(grid, a.displacement_name)
    return physical, warped, has_domain


def process_one(desc, load, a, reader_kind, distances, planes):
    grid, warped, has_domain = prepare(load(desc), a, reader_kind); rows = []
    for distance, z in zip(distances, planes):
        slc = grid.slice(normal=(0, 0, 1), origin=(0, 0, float(z)))
        if slc.n_points == 0:
            area = cnet = cavg = math.nan
        else:
            area, cnet = integrate_surface(slc, a.concentration_name); cavg = cnet / area if area > 0 else math.nan
        rows.append({"snapshot_index": desc[0], "time_s": desc[1], "distance_from_FM_mm": float(distance), "z_plane_mm": float(z), "area_mm2": area, "c_net": cnet, "c_avg": cavg, "geometry_warped": warped, "source_file": "" if desc[2] is None else str(desc[2]), "domain_field_present": has_domain})
    return rows


def arrange(rows, selected):
    d = np.array(sorted({r["distance_from_FM_mm"] for r in rows}), dtype=float); order = {x[0]: i for i, x in enumerate(selected)}
    rows.sort(key=lambda r: (order[r["snapshot_index"]], r["distance_from_FM_mm"])); shape = (len(selected), len(d))
    return np.array([x[1] for x in selected]), d, {key: np.asarray([r[key] for r in rows], dtype=float).reshape(shape) for key in ("area_mm2", "c_net", "c_avg")}


def wide_csv(path, times, distances, values):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream); writer.writerow(["time_seconds", *[f"{t:.17g}" for t in times]])
        for j, distance in enumerate(distances): writer.writerow([f"{distance:.17g}", *[f"{v:.17g}" for v in values[:, j]]])


def write_outputs(a, selected, rows, times, distances, planes, arrays, metadata):
    ensure_output_dir(a.output_dir, a.overwrite)
    with (a.output_dir / "concentration_slice_metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    wide_csv(a.output_dir / "net_concentration.csv", times, distances, arrays["c_net"]); wide_csv(a.output_dir / "avg_concentration.csv", times, distances, arrays["c_avg"])
    with (a.output_dir / "selected_timesteps.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream); writer.writerow(["snapshot_index", "time_s", "source_file"]); writer.writerows([[i, t, "" if p is None else str(p)] for i, t, p in selected])
    np.savez_compressed(a.output_dir / "concentration_slice_metrics.npz", times_s=times, distance_from_FM_mm=distances, z_plane_mm=planes, **arrays)
    (a.output_dir / "concentration_slice_metrics_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def main():
    a = parse_args(); reader_kind, desc, load, paths = source_setup(a); comm, rank, size = mpi_context() if reader_kind == "xdmf-series" else (None, 0, 1)
    reference, reference_warped, _ = prepare(load(desc[0]), a, reader_kind)
    fm_z = float(a.fm_z_mm if a.fm_z_mm is not None else math.floor(reference.bounds[5])); z_bottom = float(a.z_bottom_mm if a.z_bottom_mm is not None else math.ceil(reference.bounds[4]))
    distances, planes = transport_slice_grid(fm_z, z_bottom, a.slice_spacing_mm, a.max_slices)
    positions = np.arange(len(desc))[rank::size]
    if a.snapshot_chunk:
        if size > 1: raise ValueError("Do not combine MPI and --snapshot-chunk")
        positions = parse_chunk(a.snapshot_chunk, len(desc))
    rows = []
    for pos in positions: rows.extend(process_one(desc[int(pos)], load, a, reader_kind, distances, planes))
    if comm is not None and size > 1:
        gathered = comm.gather(rows, root=0)
        if rank != 0: return
        rows = [r for part in gathered for r in part]
    selected = [desc[int(pos)] for pos in positions] if a.snapshot_chunk else desc
    if not rows: raise RuntimeError("Selected chunk contains no snapshots")
    times, distances, arrays = arrange(rows, selected)
    signature_parameters = {"reader": reader_kind, "concentration_name": a.concentration_name, "displacement_name": a.displacement_name, "domain_name": a.domain_name, "concentration_group": a.concentration_group, "displacement_group": a.displacement_group, "fluid_tag": a.fluid_tag, "slice_spacing_mm": a.slice_spacing_mm, "fm_z_mm": a.fm_z_mm, "z_bottom_mm": a.z_bottom_mm, "indices": a.indices, "first": a.first, "last": a.last, "stride": a.stride, "dt_s": a.dt_s, "period_s": a.period_s}
    metadata = {"schema_version": 1, "reader": reader_kind, "source_signature": source_signature(paths, signature_parameters), "definitions": {"c_net": "cross-sectional concentration integral integral_A c dA; this is not global 3D mass", "c_avg": "c_net/A on the physical fluid-domain plane slice"}, "units": {"coordinates": "mm", "area": "mm^2", "concentration": a.concentration_unit, "c_net": f"{a.concentration_unit}*mm^2"}, "slice_convention": {"spacing_mm": a.slice_spacing_mm, "fm_z_mm": fm_z, "z_bottom_mm": z_bottom, "distances": "integer multiples of spacing increasing caudally from the FM", "zero_at_plot_top": True}, "geometry": {"reference_snapshot_warped": reference_warped, "rigid_wall_rule": "stored geometry used when displacement is absent"}, "snapshot_chunk": a.snapshot_chunk, "field_paths": {"concentration": a.concentration_name if reader_kind == "xdmf-series" else f"/{a.concentration_group}/INDEX", "displacement": a.displacement_name if reader_kind == "xdmf-series" else f"/{a.displacement_group}/INDEX", "domain": a.domain_name}}
    if a.snapshot_chunk:
        a.output_dir = a.output_dir / f"partial_{a.snapshot_chunk.replace('/', '_of_')}"
    write_outputs(a, selected, rows, times, distances, planes, arrays, metadata)
    if not a.no_plots:
        from plot_concentration_profiles import plot_all
        plot_all(a.output_dir, a.profile_times_s, a.cavg_limits, a.cnet_limits)
    print(f"Wrote {len(rows)} slice records to {a.output_dir}")


if __name__ == "__main__":
    try: main()
    except Exception as exc: print(f"ERROR: {exc}", file=sys.stderr); raise
