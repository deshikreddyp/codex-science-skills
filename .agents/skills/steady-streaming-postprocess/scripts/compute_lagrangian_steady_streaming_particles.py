#!/usr/bin/env python3
"""Compute phase-averaged Lagrangian steady streaming by particle tracking."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from snapshot_series import make_reader
from steady_streaming_common import (
    build_pyvista_grid,
    compact_tetra_data,
    parse_indices,
    read_mesh,
    sample_vector,
    select_cells,
    write_field_h5_xdmf,
    write_json,
)


NODE_FIELDS = [
    "local_vertex",
    "global_vertex",
    "x_mm",
    "y_mm",
    "z_mm",
    "vss_x_mm_s",
    "vss_y_mm_s",
    "vss_z_mm_s",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute the Parras/Sanchez-style Lagrangian steady-streaming field: seed fluid "
            "nodes at multiple starting phases, advect each set for one cycle, assign net drift "
            "to trajectory-mean positions, interpolate to a common mesh, and average phases."
        )
    )
    parser.add_argument("--mesh-h5", type=Path, required=True)
    parser.add_argument("--input", required=True, help="HDF5 file or filename template containing {index}.")
    parser.add_argument("--indices", required=True, help="Ordered snapshots over one cycle, e.g. 0:64.")
    parser.add_argument("--period-s", type=float, default=0.8)
    parser.add_argument("--reader", choices=("h5-array", "fenics-function"), default="fenics-function")
    parser.add_argument("--velocity-path", default="/velocity/{index}")
    parser.add_argument("--displacement-path", default=None)
    parser.add_argument("--fluid-tag", type=int, default=258514)
    parser.add_argument("--all-cells", action="store_true")
    parser.add_argument("--geometry-path", default=None)
    parser.add_argument("--topology-path", default=None)
    parser.add_argument("--cell-tags-path", default="/domains/values")
    parser.add_argument("--integrator", choices=("euler", "heun"), default="heun")
    parser.add_argument("--start-count", type=int, default=20)
    parser.add_argument(
        "--start-run",
        type=int,
        default=None,
        help="Compute only one phase number in [0,start-count); use for a SLURM array.",
    )
    parser.add_argument("--idw-neighbors", type=int, default=20)
    parser.add_argument(
        "--minimum-valid-fraction",
        type=float,
        default=0.0,
        help="Fail a phase when fewer than this fraction of seeded particles remains valid.",
    )
    parser.add_argument(
        "--maximum-nearest-distance-mm",
        type=float,
        default=None,
        help="Fail a phase when its maximum nearest-source IDW distance exceeds this value.",
    )
    parser.add_argument(
        "--phase-output-dir",
        type=Path,
        default=None,
        help="Required with --start-run; optional directory for per-phase fields in a full run.",
    )
    output_help = "Required for a full run; derived from --phase-output-dir for --start-run workers."
    parser.add_argument("--output-csv", type=Path, default=None, help=output_help)
    parser.add_argument("--output-vtp", type=Path, default=None, help=output_help)
    parser.add_argument("--field-h5", type=Path, default=None, help=output_help)
    parser.add_argument("--field-xdmf", type=Path, default=None, help=output_help)
    parser.add_argument("--field-name", default="lagrangian_cycle_mean_velocity")
    parser.add_argument("--metadata", type=Path, default=None)
    return parser.parse_args()


def resolve_outputs(args: argparse.Namespace) -> None:
    if args.start_run is not None:
        if args.phase_output_dir is None:
            raise ValueError("--start-run requires --phase-output-dir for collision-safe phase output")
        stem = args.phase_output_dir / f"lss_phase_{args.start_run:03d}"
        args.output_csv = args.output_csv or stem.with_name(stem.name + "_nodes.csv")
        args.output_vtp = args.output_vtp or stem.with_suffix(".vtp")
        args.field_h5 = args.field_h5 or stem.with_suffix(".h5")
        args.field_xdmf = args.field_xdmf or stem.with_suffix(".xdmf")
        args.metadata = args.metadata or stem.with_name(stem.name + "_metadata.json")
        return
    missing = [
        name
        for name in ("output_csv", "output_vtp", "field_h5", "field_xdmf")
        if getattr(args, name) is None
    ]
    if missing:
        flags = ", ".join("--" + name.replace("_", "-") for name in missing)
        raise ValueError(f"The full phase-average run requires {flags}")


def phase_h5_attributes(
    args: argparse.Namespace,
    phase: int,
    phase_meta: dict[str, float | int],
    indices: list[int],
) -> dict[str, object]:
    return {
        "steady_streaming_kind": "lagrangian_phase",
        "method_version": 1,
        "phase_index": phase,
        "start_snapshot_offset": int(phase_meta["start_snapshot_offset"]),
        "start_snapshot_index": int(phase_meta["start_snapshot_index"]),
        "start_count": args.start_count,
        "period_s": args.period_s,
        "mesh_h5": str(args.mesh_h5.resolve()),
        "input_template": args.input,
        "reader": args.reader,
        "snapshot_indices_json": json.dumps(indices, separators=(",", ":")),
        "velocity_path": args.velocity_path,
        "displacement_path": args.displacement_path or "",
        "integrator": args.integrator,
        "idw_neighbors": args.idw_neighbors,
        "valid_fraction": float(phase_meta["valid_fraction"]),
        "nearest_distance_max_mm": float(phase_meta["nearest_distance_max_mm"]),
        "nearest_distance_median_mm": float(phase_meta["nearest_distance_median_mm"]),
        "heun_euler_fallback_updates": int(phase_meta["heun_euler_fallback_updates"]),
    }


def read_vector(reader, index: int, path: str, n_points: int, name: str) -> np.ndarray:
    values = reader.read_vector(index, path)
    if values.shape != (n_points, 3):
        raise ValueError(f"{name} shape {values.shape} does not match mesh ({n_points}, 3)")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} snapshot {index} contains nonfinite values")
    return values


def snapshot_grid(
    reader,
    mesh_points: np.ndarray,
    topology_global: np.ndarray,
    expected_used: np.ndarray,
    index: int,
    velocity_path: str,
    displacement_path: str | None,
):
    n_points = len(mesh_points)
    velocity = read_vector(reader, index, velocity_path, n_points, "Velocity")
    displacement = (
        np.zeros_like(mesh_points)
        if displacement_path is None
        else read_vector(reader, index, displacement_path, n_points, "Displacement")
    )
    points, topology, fields, used = compact_tetra_data(
        mesh_points + displacement, topology_global, {"velocity": velocity}
    )
    if not np.array_equal(used, expected_used):
        raise RuntimeError("Fluid vertex selection changed unexpectedly")
    return build_pyvista_grid(points, topology, fields), points, topology


def interpolate_idw(
    source_points: np.ndarray,
    source_vectors: np.ndarray,
    target_points: np.ndarray,
    neighbors: int,
) -> tuple[np.ndarray, dict[str, float]]:
    try:
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise RuntimeError("Phase remapping requires scipy.spatial.cKDTree") from exc
    if len(source_points) == 0:
        raise ValueError("No valid particle trajectories are available for phase remapping")
    k = min(neighbors, len(source_points))
    tree = cKDTree(source_points)
    distances, indices = tree.query(target_points, k=k, workers=-1)
    if k == 1:
        distances = distances[:, None]
        indices = indices[:, None]
    vectors = source_vectors[indices]
    mapped = np.empty((len(target_points), 3), dtype=np.float64)
    exact = distances[:, 0] <= 1.0e-12
    mapped[exact] = vectors[exact, 0]
    nonexact = ~exact
    if np.any(nonexact):
        weights = 1.0 / np.maximum(distances[nonexact] ** 2, 1.0e-24)
        mapped[nonexact] = np.sum(weights[:, :, None] * vectors[nonexact], axis=1) / np.sum(
            weights, axis=1
        )[:, None]
    return mapped, {
        "nearest_distance_max_mm": float(np.max(distances[:, 0])),
        "nearest_distance_median_mm": float(np.median(distances[:, 0])),
        "idw_neighbors": int(k),
    }


def advect_phase(
    args: argparse.Namespace,
    reader,
    mesh_points: np.ndarray,
    topology_global: np.ndarray,
    used_vertices: np.ndarray,
    indices: list[int],
    start_offset: int,
    target_points: np.ndarray,
) -> tuple[np.ndarray, dict[str, float | int]]:
    schedule = [indices[(start_offset + step) % len(indices)] for step in range(len(indices))]
    current_grid, start_points, _ = snapshot_grid(
        reader,
        mesh_points,
        topology_global,
        used_vertices,
        schedule[0],
        args.velocity_path,
        args.displacement_path,
    )
    particles = start_points.copy()
    initial = start_points.copy()
    accumulated = start_points.copy()
    active = np.ones(len(particles), dtype=bool)
    heun_fallback_updates = 0
    dt = args.period_s / float(len(indices))

    for step, index in enumerate(schedule):
        next_index = schedule[(step + 1) % len(schedule)]
        if step + 1 == len(schedule):
            next_grid, _, _ = snapshot_grid(
                reader,
                mesh_points,
                topology_global,
                used_vertices,
                schedule[0],
                args.velocity_path,
                args.displacement_path,
            )
        else:
            next_grid, _, _ = snapshot_grid(
                reader,
                mesh_points,
                topology_global,
                used_vertices,
                next_index,
                args.velocity_path,
                args.displacement_path,
            )
        active_ids = np.flatnonzero(active)
        if len(active_ids) == 0:
            break
        velocity_now, valid_now = sample_vector(current_grid, particles[active_ids], "velocity")
        step_ids = active_ids[valid_now]
        active[active_ids[~valid_now]] = False
        current_velocity = velocity_now[valid_now]
        if args.integrator == "euler":
            particles[step_ids] += dt * current_velocity
        else:
            predictor = particles[step_ids] + dt * current_velocity
            velocity_next, valid_next = sample_vector(next_grid, predictor, "velocity")
            accepted_ids = step_ids[valid_next]
            particles[accepted_ids] += 0.5 * dt * (
                current_velocity[valid_next] + velocity_next[valid_next]
            )
            fallback_ids = step_ids[~valid_next]
            if len(fallback_ids):
                fallback_positions = particles[fallback_ids] + dt * current_velocity[~valid_next]
                _, valid_fallback = sample_vector(next_grid, fallback_positions, "velocity")
                particles[fallback_ids[valid_fallback]] = fallback_positions[valid_fallback]
                active[fallback_ids[~valid_fallback]] = False
                heun_fallback_updates += int(np.count_nonzero(valid_fallback))
        accumulated[active] += particles[active]
        current_grid = next_grid
        if step == 0 or (step + 1) % 8 == 0 or step + 1 == len(schedule):
            print(
                f"  step {step+1}/{len(schedule)} ({index} -> {next_index}), "
                f"active={np.mean(active):.2%}",
                flush=True,
            )

    mean_positions = accumulated / float(len(schedule) + 1)
    drift = (particles - initial) / args.period_s
    mapped, interpolation_meta = interpolate_idw(
        mean_positions[active], drift[active], target_points, args.idw_neighbors
    )
    meta: dict[str, float | int] = {
        "start_snapshot_offset": int(start_offset),
        "start_snapshot_index": int(schedule[0]),
        "valid_particles": int(np.count_nonzero(active)),
        "total_particles": int(len(active)),
        "valid_fraction": float(np.mean(active)),
        "heun_euler_fallback_updates": int(heun_fallback_updates),
    }
    meta.update(interpolation_meta)
    return mapped, meta


def write_nodes_csv(
    path: Path,
    used_vertices: np.ndarray,
    points: np.ndarray,
    velocity: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NODE_FIELDS)
        writer.writeheader()
        for local_index in range(len(points)):
            writer.writerow(
                {
                    "local_vertex": local_index,
                    "global_vertex": int(used_vertices[local_index]),
                    "x_mm": float(points[local_index, 0]),
                    "y_mm": float(points[local_index, 1]),
                    "z_mm": float(points[local_index, 2]),
                    "vss_x_mm_s": float(velocity[local_index, 0]),
                    "vss_y_mm_s": float(velocity[local_index, 1]),
                    "vss_z_mm_s": float(velocity[local_index, 2]),
                }
            )


def main() -> None:
    args = parse_args()
    if args.period_s <= 0.0:
        raise ValueError("--period-s must be positive")
    indices = parse_indices(args.indices)
    if args.start_count < 1 or args.start_count > len(indices):
        raise ValueError("--start-count must be between 1 and the number of snapshots")
    if args.idw_neighbors < 1:
        raise ValueError("--idw-neighbors must be positive")
    if not 0.0 <= args.minimum_valid_fraction <= 1.0:
        raise ValueError("--minimum-valid-fraction must be in [0, 1]")
    if args.maximum_nearest_distance_mm is not None and args.maximum_nearest_distance_mm <= 0.0:
        raise ValueError("--maximum-nearest-distance-mm must be positive")
    if args.start_run is not None and not 0 <= args.start_run < args.start_count:
        raise ValueError("--start-run must satisfy 0 <= start-run < start-count")
    resolve_outputs(args)

    mesh = read_mesh(args.mesh_h5, args.geometry_path, args.topology_path, args.cell_tags_path)
    selected_topology, selected_tags = select_cells(mesh, None if args.all_cells else args.fluid_tag)
    _, local_topology, _, used_vertices = compact_tetra_data(mesh.points, selected_topology)
    reader = make_reader(args.reader, args.mesh_h5, args.input)
    _, target_points, target_topology = snapshot_grid(
        reader,
        mesh.points,
        selected_topology,
        used_vertices,
        indices[0],
        args.velocity_path,
        args.displacement_path,
    )
    if not np.array_equal(local_topology, target_topology):
        raise RuntimeError("Local topology mismatch")

    start_offsets = [phase * len(indices) // args.start_count for phase in range(args.start_count)]
    exact_phase_spacing = len(indices) % args.start_count == 0
    if not exact_phase_spacing:
        print(
            f"warning: {len(indices)} snapshots are not divisible by {args.start_count} starting phases; "
            "using approximately uniform floor-spaced offsets",
            flush=True,
        )
    selected_phases = [args.start_run] if args.start_run is not None else list(range(args.start_count))
    accumulator = np.zeros((len(target_points), 3), dtype=np.float64)
    phase_metadata: list[dict[str, float | int]] = []
    for count, phase in enumerate(selected_phases, start=1):
        print(
            f"phase {count}/{len(selected_phases)}: phase={phase}, "
            f"start offset={start_offsets[phase]}",
            flush=True,
        )
        mapped, phase_meta = advect_phase(
            args,
            reader,
            mesh.points,
            selected_topology,
            used_vertices,
            indices,
            start_offsets[phase],
            target_points,
        )
        phase_meta["phase"] = phase
        if float(phase_meta["valid_fraction"]) < args.minimum_valid_fraction:
            raise RuntimeError(
                f"Phase {phase} valid fraction {phase_meta['valid_fraction']:.2%} is below "
                f"{args.minimum_valid_fraction:.2%}"
            )
        if (
            args.maximum_nearest_distance_mm is not None
            and float(phase_meta["nearest_distance_max_mm"]) > args.maximum_nearest_distance_mm
        ):
            raise RuntimeError(
                f"Phase {phase} nearest remap distance {phase_meta['nearest_distance_max_mm']:.6g} mm "
                f"exceeds {args.maximum_nearest_distance_mm:.6g} mm"
            )
        phase_metadata.append(phase_meta)
        accumulator += mapped
        if args.phase_output_dir is not None:
            args.phase_output_dir.mkdir(parents=True, exist_ok=True)
            phase_h5 = args.phase_output_dir / f"lss_phase_{phase:03d}.h5"
            phase_xdmf = args.phase_output_dir / f"lss_phase_{phase:03d}.xdmf"
            if phase_h5.resolve() != args.field_h5.resolve():
                write_field_h5_xdmf(
                    phase_h5,
                    phase_xdmf,
                    target_points,
                    local_topology,
                    mapped,
                    args.field_name,
                    selected_tags,
                    phase_h5_attributes(args, phase, phase_meta, indices),
                )

    mean_lss = accumulator / float(len(selected_phases))
    write_field_h5_xdmf(
        args.field_h5,
        args.field_xdmf,
        target_points,
        local_topology,
        mean_lss,
        args.field_name,
        selected_tags,
        (
            phase_h5_attributes(args, int(args.start_run), phase_metadata[0], indices)
            if args.start_run is not None
            else {
                "steady_streaming_kind": "lagrangian_phase_average",
                "method_version": 1,
                "phase_count": len(selected_phases),
                "start_count": args.start_count,
                "period_s": args.period_s,
                "mesh_h5": str(args.mesh_h5.resolve()),
                "input_template": args.input,
                "reader": args.reader,
                "snapshot_indices_json": json.dumps(indices, separators=(",", ":")),
                "velocity_path": args.velocity_path,
                "displacement_path": args.displacement_path or "",
                "integrator": args.integrator,
                "idw_neighbors": args.idw_neighbors,
            }
        ),
    )
    write_nodes_csv(args.output_csv, used_vertices, target_points, mean_lss)

    import pyvista as pv

    cloud = pv.PolyData(target_points)
    cloud.point_data[args.field_name] = mean_lss
    cloud.point_data["global_vertex"] = used_vertices
    args.output_vtp.parent.mkdir(parents=True, exist_ok=True)
    cloud.save(args.output_vtp)

    metadata_path = args.metadata or args.output_csv.with_name(args.output_csv.stem + "_metadata.json")
    write_json(
        metadata_path,
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "definition": (
                "phase-average of one-cycle particle drift, remapped from each trajectory-mean "
                "position to the first-snapshot fluid mesh by local inverse-distance weighting"
            ),
            "method_reference": "/home/dputluru/fenics_run/steady_streaming_better_toy",
            "mesh_h5": str(args.mesh_h5),
            "input": args.input,
            "reader": args.reader,
            "indices": indices,
            "period_s": args.period_s,
            "dt_s": args.period_s / len(indices),
            "integrator": args.integrator,
            "velocity_path": args.velocity_path,
            "displacement_path": args.displacement_path,
            "fluid_tag": None if args.all_cells else args.fluid_tag,
            "start_count_requested": args.start_count,
            "start_run": args.start_run,
            "start_offsets": start_offsets,
            "exact_uniform_phase_spacing": exact_phase_spacing,
            "phases_computed": selected_phases,
            "phase_metadata": phase_metadata,
            "minimum_valid_fraction": args.minimum_valid_fraction,
            "maximum_nearest_distance_mm": args.maximum_nearest_distance_mm,
            "field_name": args.field_name,
            "output_csv": str(args.output_csv),
            "output_vtp": str(args.output_vtp),
            "field_h5": str(args.field_h5),
            "field_xdmf": str(args.field_xdmf),
        },
    )
    print(f"wrote {args.field_xdmf}", flush=True)
    print(f"wrote {args.output_csv}", flush=True)
    print(f"wrote {args.output_vtp}", flush=True)


if __name__ == "__main__":
    main()
