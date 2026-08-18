#!/usr/bin/env python3
"""Integrate signed steady-streaming strength on zero-thickness axial slices."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

from steady_streaming_common import (
    MeshData,
    build_pyvista_grid,
    compact_tetra_data,
    h5_key,
    read_mesh,
    select_cells,
    slice_positions,
    strength_rows,
    write_csv,
    write_json,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Slice a tetrahedral steady-streaming field every 1 mm by default. Compute "
            "cranial integral(max(vz,0),dA), signed caudal integral(min(vz,0),dA), net, "
            "cross-sectional area, and area-normalized directional velocities."
        )
    )
    parser.add_argument("--field-h5", type=Path, required=True)
    parser.add_argument(
        "--mesh-h5",
        type=Path,
        default=None,
        help="Optional source of cell tags when the field HDF5 has none.",
    )
    parser.add_argument("--geometry-path", default="/Mesh/0/mesh/geometry")
    parser.add_argument("--topology-path", default="/Mesh/0/mesh/topology")
    parser.add_argument("--velocity-path", default="/VisualisationVector/0")
    parser.add_argument("--field-cell-tags-path", default="/Mesh/0/mesh/subdomain")
    parser.add_argument("--mesh-cell-tags-path", default="/domains/values")
    parser.add_argument("--fluid-tag", type=int, default=258514)
    parser.add_argument("--all-cells", action="store_true")
    parser.add_argument("--slice-spacing-mm", type=float, default=1.0)
    parser.add_argument(
        "--slice-chunk",
        default=None,
        help="Parallel slice partition K/N. Each job writes a separate CSV for later row concatenation.",
    )
    parser.add_argument(
        "--axial-sign",
        type=float,
        choices=(-1.0, 1.0),
        default=1.0,
        help="Use +1 when positive axial velocity is cranial; use -1 to reverse the convention.",
    )
    parser.add_argument(
        "--fm-z-mm",
        type=float,
        default=None,
        help="Optional z coordinate of the FM distance origin; default is the fluid field's z_max.",
    )
    parser.add_argument("--case", default="case")
    parser.add_argument("--label", default="Case")
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, default=None)
    return parser.parse_args()


def parse_chunk(spec: str | None) -> tuple[int, int]:
    if spec is None:
        return 0, 1
    try:
        rank_text, count_text = spec.split("/", maxsplit=1)
        rank, count = int(rank_text), int(count_text)
    except (ValueError, AttributeError) as exc:
        raise ValueError("--slice-chunk must have form K/N") from exc
    if count < 1 or rank < 0 or rank >= count:
        raise ValueError("--slice-chunk requires 0 <= K < N")
    return rank, count


def main() -> None:
    args = parse_args()
    field_mesh = read_mesh(
        args.field_h5,
        args.geometry_path,
        args.topology_path,
        args.field_cell_tags_path,
    )
    cell_tags = field_mesh.cell_tags
    if cell_tags is None and args.mesh_h5 is not None:
        external = read_mesh(args.mesh_h5, cell_tags_path=args.mesh_cell_tags_path)
        if external.topology.shape != field_mesh.topology.shape or not np.array_equal(
            external.topology, field_mesh.topology
        ):
            raise ValueError("External mesh topology does not exactly match the field topology")
        cell_tags = external.cell_tags
    mesh = MeshData(field_mesh.points, field_mesh.topology, cell_tags)
    selected_topology, selected_tags = select_cells(mesh, None if args.all_cells else args.fluid_tag)

    with h5py.File(args.field_h5, "r") as h5:
        key = h5_key(args.velocity_path)
        if key not in h5:
            raise KeyError(f"Velocity dataset {args.velocity_path} not found in {args.field_h5}")
        velocity = np.asarray(h5[key][:], dtype=np.float64)
    if velocity.shape != (len(mesh.points), 3):
        raise ValueError(f"Velocity shape {velocity.shape} does not match geometry ({len(mesh.points)}, 3)")
    velocity = velocity.copy()
    velocity[:, 2] *= args.axial_sign

    points, topology, fields, _ = compact_tetra_data(
        mesh.points, selected_topology, {"steady_streaming_velocity": velocity}
    )
    selected_velocity = fields["steady_streaming_velocity"]
    nonfinite = int(selected_velocity.size - np.count_nonzero(np.isfinite(selected_velocity)))
    if nonfinite:
        raise ValueError(
            f"Selected fluid field contains {nonfinite} nonfinite velocity components; "
            "restrict/remap the field before slice integration"
        )
    grid = build_pyvista_grid(
        points,
        topology,
        fields,
        {"subdomain": selected_tags},
    )
    rank, partitions = parse_chunk(args.slice_chunk)
    all_z_values = slice_positions(points, args.slice_spacing_mm)
    selected_z_values = all_z_values[rank::partitions]
    if len(selected_z_values) == 0:
        raise ValueError("Slice chunk is empty")
    z_upper = float(np.max(points[:, 2])) if args.fm_z_mm is None else args.fm_z_mm
    if z_upper < float(np.max(points[:, 2])) - 1.0e-9:
        raise ValueError("--fm-z-mm cannot be below the selected fluid field's z_max")
    rows = strength_rows(
        grid,
        "steady_streaming_velocity",
        args.slice_spacing_mm,
        args.case,
        args.label,
        2,
        z_values=selected_z_values,
        z_upper=z_upper,
    )
    if not rows:
        raise RuntimeError("No nondegenerate slice intersections were produced")
    write_csv(args.output_csv, rows)
    metadata_path = args.metadata or args.output_csv.with_name(args.output_csv.stem + "_metadata.json")
    write_json(
        metadata_path,
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "field_h5": str(args.field_h5),
            "mesh_h5_for_tags": str(args.mesh_h5) if args.mesh_h5 else None,
            "velocity_path": args.velocity_path,
            "fluid_tag": None if args.all_cells else args.fluid_tag,
            "slice_spacing_mm": args.slice_spacing_mm,
            "slice_chunk": args.slice_chunk,
            "slice_definition": "zero-thickness PyVista plane normal to z",
            "distance_from_fm_mm": "fm_z_mm minus z_slice",
            "fm_z_mm": z_upper,
            "cranial_strength": "integral(max(v_axial, 0), dA)",
            "caudal_strength": "integral(min(v_axial, 0), dA), stored signed negative",
            "mixed_sign_triangles": "exact integration of positive/negative parts of the linear nodal field",
            "axial_component": 2,
            "axial_sign": args.axial_sign,
            "case": args.case,
            "label": args.label,
            "output_csv": str(args.output_csv),
            "rows": len(rows),
        },
    )
    print(f"wrote {args.output_csv}", flush=True)
    print(f"wrote {metadata_path}", flush=True)


if __name__ == "__main__":
    main()
