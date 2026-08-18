#!/usr/bin/env python3
"""Compute a fixed-spatial Eulerian or fixed-node ALE cycle-mean velocity."""

from __future__ import annotations

import argparse
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute steady streaming as either a true Eulerian cycle mean at fixed physical "
            "points or an ALE nodal mean at fixed mesh identities. Snapshot intervals are "
            "assumed uniform and the periodic endpoint must not be duplicated."
        )
    )
    parser.add_argument("--mesh-h5", type=Path, required=True)
    parser.add_argument(
        "--input",
        required=True,
        help="Snapshot HDF5 file or filename template containing {index}; quote shell braces.",
    )
    parser.add_argument("--indices", required=True, help="For example 0:64 or 0:640:3")
    parser.add_argument("--reader", choices=("h5-array", "fenics-function"), default="fenics-function")
    parser.add_argument("--velocity-path", default="/velocity/{index}")
    parser.add_argument(
        "--displacement-path",
        default=None,
        help="Optional nodal displacement path/template. Omit for rigid-wall meshes.",
    )
    parser.add_argument("--frame", choices=("eulerian", "ale"), default="ale")
    parser.add_argument(
        "--reference-index",
        type=int,
        default=None,
        help="Snapshot geometry used as the fixed Eulerian grid and output geometry; default is first index.",
    )
    parser.add_argument("--fluid-tag", type=int, default=258514)
    parser.add_argument("--all-cells", action="store_true", help="Do not filter cells by subdomain tag.")
    parser.add_argument("--geometry-path", default=None)
    parser.add_argument("--topology-path", default=None)
    parser.add_argument("--cell-tags-path", default="/domains/values")
    parser.add_argument("--field-name", default=None)
    parser.add_argument("--minimum-valid-fraction", type=float, default=1.0)
    parser.add_argument("--output-h5", type=Path, default=None)
    parser.add_argument("--output-xdmf", type=Path, default=None)
    parser.add_argument("--metadata", type=Path, default=None)
    parser.add_argument(
        "--snapshot-chunk",
        default=None,
        help="Parallel ALE partial as K/N, selecting indices[K::N]. Requires --partial-output.",
    )
    parser.add_argument("--partial-output", type=Path, default=None)
    parser.add_argument(
        "--combine-partials",
        nargs="+",
        type=Path,
        default=None,
        help="Combine ALE sum/count NPZ files made with --snapshot-chunk.",
    )
    return parser.parse_args()


def read_displacement(reader, index: int, path: str | None, n_points: int) -> np.ndarray:
    if path is None:
        return np.zeros((n_points, 3), dtype=np.float64)
    displacement = reader.read_vector(index, path)
    if displacement.shape != (n_points, 3):
        raise ValueError(f"Displacement shape {displacement.shape} does not match mesh ({n_points}, 3)")
    if not np.all(np.isfinite(displacement)):
        raise ValueError(f"Displacement snapshot {index} contains nonfinite values")
    return displacement


def read_velocity(reader, index: int, path: str, n_points: int) -> np.ndarray:
    velocity = reader.read_vector(index, path)
    if velocity.shape != (n_points, 3):
        raise ValueError(f"Velocity shape {velocity.shape} does not match mesh ({n_points}, 3)")
    if not np.all(np.isfinite(velocity)):
        raise ValueError(f"Velocity snapshot {index} contains nonfinite values")
    return velocity


def parse_chunk(spec: str) -> tuple[int, int]:
    try:
        rank_text, count_text = spec.split("/", maxsplit=1)
        rank, count = int(rank_text), int(count_text)
    except (ValueError, AttributeError) as exc:
        raise ValueError("--snapshot-chunk must have form K/N") from exc
    if count < 1 or rank < 0 or rank >= count:
        raise ValueError("--snapshot-chunk requires 0 <= K < N")
    return rank, count


def ale_sum(reader, indices: list[int], velocity_path: str, n_points: int) -> np.ndarray:
    total = np.zeros((n_points, 3), dtype=np.float64)
    for count, index in enumerate(indices, start=1):
        total += read_velocity(reader, index, velocity_path, n_points)
        if count == 1 or count % 16 == 0 or count == len(indices):
            print(f"velocity snapshot {count}/{len(indices)} (index {index})", flush=True)
    return total


def main() -> None:
    args = parse_args()
    if not 0.0 < args.minimum_valid_fraction <= 1.0:
        raise ValueError("--minimum-valid-fraction must be in (0, 1]")
    if args.snapshot_chunk is not None and args.combine_partials:
        raise ValueError("--snapshot-chunk and --combine-partials are mutually exclusive")
    if args.snapshot_chunk is None and (args.output_h5 is None or args.output_xdmf is None):
        raise ValueError("Provide --output-h5 and --output-xdmf unless writing a snapshot partial")
    indices = parse_indices(args.indices)
    reference_index = args.reference_index if args.reference_index is not None else indices[0]
    mesh = read_mesh(args.mesh_h5, args.geometry_path, args.topology_path, args.cell_tags_path)
    fluid_tag = None if args.all_cells else args.fluid_tag
    selected_topology, selected_tags = select_cells(mesh, fluid_tag)
    reader = make_reader(args.reader, args.mesh_h5, args.input)
    n_points = len(mesh.points)

    if args.combine_partials and args.frame != "ale":
        raise ValueError("Partial-file combination is supported only for --frame ale")
    if args.snapshot_chunk:
        if args.frame != "ale" or args.partial_output is None:
            raise ValueError("--snapshot-chunk requires --frame ale and --partial-output")
        rank, count = parse_chunk(args.snapshot_chunk)
        chunk_indices = indices[rank::count]
        if not chunk_indices:
            raise ValueError("This snapshot chunk is empty")
        velocity_sum = ale_sum(reader, chunk_indices, args.velocity_path, n_points)
        args.partial_output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            args.partial_output,
            velocity_sum=velocity_sum,
            count=np.array([len(chunk_indices)], dtype=np.int64),
            indices=np.asarray(chunk_indices, dtype=np.int64),
            mesh_h5=np.array(str(args.mesh_h5.resolve())),
            input_template=np.array(args.input),
            reader=np.array(args.reader),
            velocity_path=np.array(args.velocity_path),
        )
        print(f"wrote ALE partial: {args.partial_output}", flush=True)
        return

    reference_displacement = read_displacement(reader, reference_index, args.displacement_path, n_points)
    reference_points = mesh.points + reference_displacement
    output_points, output_topology, _, used_vertices = compact_tetra_data(reference_points, selected_topology)

    if args.combine_partials:
        total = np.zeros((n_points, 3), dtype=np.float64)
        total_count = 0
        combined_indices: list[int] = []
        for partial in args.combine_partials:
            with np.load(partial) as data:
                provenance = {
                    "mesh_h5": str(data["mesh_h5"].item()),
                    "input_template": str(data["input_template"].item()),
                    "reader": str(data["reader"].item()),
                    "velocity_path": str(data["velocity_path"].item()),
                }
                expected = {
                    "mesh_h5": str(args.mesh_h5.resolve()),
                    "input_template": args.input,
                    "reader": args.reader,
                    "velocity_path": args.velocity_path,
                }
                if provenance != expected:
                    raise ValueError(
                        f"Partial provenance mismatch in {partial}: {provenance} != {expected}"
                    )
                if data["velocity_sum"].shape != total.shape:
                    raise ValueError(f"Partial shape mismatch in {partial}")
                if int(data["count"][0]) != len(data["indices"]):
                    raise ValueError(f"Partial count/index mismatch in {partial}")
                total += data["velocity_sum"]
                total_count += int(data["count"][0])
                combined_indices.extend(int(value) for value in data["indices"])
        if sorted(combined_indices) != sorted(indices):
            raise ValueError("Combined partial indices do not exactly match --indices")
        mean_full = total / float(total_count)
        mean_output = mean_full[used_vertices]
        valid_fraction = np.ones(len(mean_output), dtype=np.float64)
    elif args.frame == "ale":
        mean_full = ale_sum(reader, indices, args.velocity_path, n_points) / float(len(indices))
        mean_output = mean_full[used_vertices]
        valid_fraction = np.ones(len(mean_output), dtype=np.float64)
    else:
        accumulator = np.zeros((len(output_points), 3), dtype=np.float64)
        valid_count = np.zeros(len(output_points), dtype=np.int64)
        for count, index in enumerate(indices, start=1):
            velocity = read_velocity(reader, index, args.velocity_path, n_points)
            displacement = read_displacement(reader, index, args.displacement_path, n_points)
            current_points = mesh.points + displacement
            compact_points, compact_topology, fields, current_used = compact_tetra_data(
                current_points, selected_topology, {"velocity": velocity}
            )
            if not np.array_equal(current_used, used_vertices):
                raise RuntimeError("Fluid vertex selection changed unexpectedly")
            grid = build_pyvista_grid(compact_points, compact_topology, fields)
            sampled, valid = sample_vector(grid, output_points, "velocity")
            accumulator[valid] += sampled[valid]
            valid_count[valid] += 1
            if count == 1 or count % 16 == 0 or count == len(indices):
                print(
                    f"Eulerian sample {count}/{len(indices)} (index {index}), "
                    f"valid={np.mean(valid):.2%}",
                    flush=True,
                )
        valid_fraction = valid_count / float(len(indices))
        mean_output = np.full_like(accumulator, np.nan)
        accepted = (valid_count > 0) & (valid_fraction >= args.minimum_valid_fraction)
        mean_output[accepted] = accumulator[accepted] / valid_count[accepted, None]

    field_name = args.field_name or (
        "eulerian_cycle_mean_velocity" if args.frame == "eulerian" else "ale_cycle_mean_velocity"
    )
    write_field_h5_xdmf(
        args.output_h5,
        args.output_xdmf,
        output_points,
        output_topology,
        mean_output,
        field_name,
        selected_tags,
    )
    metadata_path = args.metadata or args.output_xdmf.with_name(args.output_xdmf.stem + "_metadata.json")
    write_json(
        metadata_path,
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "definition": (
                "fixed-physical-point Eulerian arithmetic cycle mean"
                if args.frame == "eulerian"
                else "fixed-mesh-identity ALE arithmetic cycle mean"
            ),
            "frame": args.frame,
            "mesh_h5": str(args.mesh_h5),
            "input": args.input,
            "reader": args.reader,
            "velocity_path": args.velocity_path,
            "displacement_path": args.displacement_path,
            "indices": indices,
            "reference_index": reference_index,
            "fluid_tag": fluid_tag,
            "field_name": field_name,
            "valid_fraction_min": float(np.min(valid_fraction)),
            "valid_fraction_median": float(np.median(valid_fraction)),
            "output_vertices": int(len(mean_output)),
            "finite_output_vertices": int(np.count_nonzero(np.all(np.isfinite(mean_output), axis=1))),
            "output_h5": str(args.output_h5),
            "output_xdmf": str(args.output_xdmf),
        },
    )
    print(f"wrote {args.output_xdmf}", flush=True)
    print(f"wrote {metadata_path}", flush=True)


if __name__ == "__main__":
    main()
