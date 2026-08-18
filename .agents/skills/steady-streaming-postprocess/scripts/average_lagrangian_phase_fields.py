#!/usr/bin/env python3
"""Average independently computed Lagrangian starting-phase fields."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

from steady_streaming_common import (
    array_sha256,
    h5_key,
    parse_indices,
    read_mesh,
    write_field_h5_xdmf,
    write_json,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Average Lagrangian phase fields produced by independent --start-run jobs while "
            "verifying common geometry, topology, cell tags, and finite velocity values."
        )
    )
    parser.add_argument("--phase-dir", type=Path, required=True)
    parser.add_argument("--indices", required=True, help="Phase numbers, for example 0:20.")
    parser.add_argument("--template", default="lss_phase_{index:03d}.h5")
    parser.add_argument("--geometry-path", default="/Mesh/0/mesh/geometry")
    parser.add_argument("--topology-path", default="/Mesh/0/mesh/topology")
    parser.add_argument("--cell-tags-path", default="/Mesh/0/mesh/subdomain")
    parser.add_argument("--velocity-path", default="/VisualisationVector/0")
    parser.add_argument("--field-name", default="lagrangian_cycle_mean_velocity")
    parser.add_argument("--minimum-valid-fraction", type=float, default=0.0)
    parser.add_argument("--maximum-nearest-distance-mm", type=float, default=None)
    parser.add_argument(
        "--allow-partial-phase-average",
        action="store_true",
        help="Allow averaging a subset instead of all phase indices 0:start_count.",
    )
    parser.add_argument("--output-h5", type=Path, required=True)
    parser.add_argument("--output-xdmf", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, default=None)
    return parser.parse_args()


def fingerprint(h5: h5py.File, attribute: str, dataset: str) -> str:
    if attribute in h5.attrs:
        value = h5.attrs[attribute]
        return value.decode() if isinstance(value, bytes) else str(value)
    return array_sha256(np.asarray(h5[h5_key(dataset)][:]))


def attribute_text(h5: h5py.File, name: str) -> str:
    value = h5.attrs.get(name, "")
    return value.decode() if isinstance(value, bytes) else str(value)


def main() -> None:
    args = parse_args()
    if not 0.0 <= args.minimum_valid_fraction <= 1.0:
        raise ValueError("--minimum-valid-fraction must be in [0, 1]")
    if args.maximum_nearest_distance_mm is not None and args.maximum_nearest_distance_mm <= 0.0:
        raise ValueError("--maximum-nearest-distance-mm must be positive")
    indices = parse_indices(args.indices)
    paths = [args.phase_dir / args.template.format(index=index, i=index) for index in indices]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing phase fields: " + ", ".join(missing))

    first_mesh = read_mesh(paths[0], args.geometry_path, args.topology_path, args.cell_tags_path)
    if first_mesh.cell_tags is None:
        raise ValueError(f"Cell tags {args.cell_tags_path} are missing from {paths[0]}")
    expected_fingerprints: dict[str, str] = {}
    expected_phase_provenance: dict[str, object] = {}
    accumulator: np.ndarray | None = None
    sources: list[dict[str, object]] = []
    phase_diagnostics: list[dict[str, object]] = []
    for count, path in enumerate(paths, start=1):
        with h5py.File(path, "r") as h5:
            kind = str(h5.attrs.get("steady_streaming_kind", ""))
            phase_index = int(h5.attrs.get("phase_index", -1))
            if kind != "lagrangian_phase" or phase_index != indices[count - 1]:
                raise ValueError(
                    f"Phase provenance mismatch in {path}: kind={kind!r}, phase_index={phase_index}"
                )
            current_phase_provenance = {
                "method_version": int(h5.attrs.get("method_version", -1)),
                "field_name": attribute_text(h5, "field_name"),
                "start_count": int(h5.attrs.get("start_count", -1)),
                "period_s": float(h5.attrs.get("period_s", np.nan)),
                "mesh_h5": attribute_text(h5, "mesh_h5"),
                "input_template": attribute_text(h5, "input_template"),
                "reader": attribute_text(h5, "reader"),
                "snapshot_indices_json": attribute_text(h5, "snapshot_indices_json"),
                "velocity_path": attribute_text(h5, "velocity_path"),
                "displacement_path": attribute_text(h5, "displacement_path"),
                "integrator": attribute_text(h5, "integrator"),
                "idw_neighbors": int(h5.attrs.get("idw_neighbors", -1)),
            }
            if not expected_phase_provenance:
                expected_phase_provenance = current_phase_provenance
            elif current_phase_provenance != expected_phase_provenance:
                raise ValueError(f"Phase-run provenance mismatch in {path}")
            diagnostics = {
                "phase": phase_index,
                "valid_fraction": float(h5.attrs.get("valid_fraction", np.nan)),
                "nearest_distance_max_mm": float(
                    h5.attrs.get("nearest_distance_max_mm", np.nan)
                ),
                "nearest_distance_median_mm": float(
                    h5.attrs.get("nearest_distance_median_mm", np.nan)
                ),
                "heun_euler_fallback_updates": int(
                    h5.attrs.get("heun_euler_fallback_updates", -1)
                ),
            }
            if not np.isfinite(diagnostics["valid_fraction"]):
                raise ValueError(f"Missing particle-validity diagnostics in {path}")
            if diagnostics["valid_fraction"] < args.minimum_valid_fraction:
                raise ValueError(
                    f"Phase {phase_index} valid fraction {diagnostics['valid_fraction']:.2%} "
                    f"is below {args.minimum_valid_fraction:.2%}"
                )
            if not np.isfinite(diagnostics["nearest_distance_max_mm"]):
                raise ValueError(f"Missing IDW-distance diagnostics in {path}")
            if (
                args.maximum_nearest_distance_mm is not None
                and diagnostics["nearest_distance_max_mm"] > args.maximum_nearest_distance_mm
            ):
                raise ValueError(
                    f"Phase {phase_index} nearest remap distance "
                    f"{diagnostics['nearest_distance_max_mm']:.6g} mm exceeds "
                    f"{args.maximum_nearest_distance_mm:.6g} mm"
                )
            current_fingerprints = {
                "geometry": fingerprint(h5, "geometry_sha256", args.geometry_path),
                "topology": fingerprint(h5, "topology_sha256", args.topology_path),
                "cell_tags": fingerprint(h5, "cell_tags_sha256", args.cell_tags_path),
            }
            if not expected_fingerprints:
                expected_fingerprints = current_fingerprints
            elif current_fingerprints != expected_fingerprints:
                raise ValueError(f"Mesh fingerprint mismatch in {path}")
            key = h5_key(args.velocity_path)
            if key not in h5:
                raise KeyError(f"Velocity dataset {args.velocity_path} is missing from {path}")
            velocity = np.asarray(h5[key][:], dtype=np.float64)
        if velocity.shape != first_mesh.points.shape:
            raise ValueError(
                f"Velocity shape {velocity.shape} in {path} does not match geometry {first_mesh.points.shape}"
            )
        if not np.all(np.isfinite(velocity)):
            raise ValueError(f"Nonfinite velocity values found in {path}")
        if accumulator is None:
            accumulator = np.zeros_like(velocity)
        accumulator += velocity
        sources.append({"phase": indices[count - 1], "path": str(path)})
        phase_diagnostics.append(diagnostics)
        print(f"phase field {count}/{len(paths)}: {path}", flush=True)

    assert accumulator is not None
    start_count = int(expected_phase_provenance["start_count"])
    if not args.allow_partial_phase_average and indices != list(range(start_count)):
        raise ValueError(
            f"A complete phase average requires indices 0:{start_count}; got {indices}. "
            "Use --allow-partial-phase-average only intentionally."
        )
    mean_velocity = accumulator / float(len(paths))
    write_field_h5_xdmf(
        args.output_h5,
        args.output_xdmf,
        first_mesh.points,
        first_mesh.topology,
        mean_velocity,
        args.field_name,
        first_mesh.cell_tags,
        {
            "steady_streaming_kind": "lagrangian_phase_average",
            "method_version": int(expected_phase_provenance["method_version"]),
            "phase_count": len(indices),
            "start_count": start_count,
            "period_s": float(expected_phase_provenance["period_s"]),
            "mesh_h5": str(expected_phase_provenance["mesh_h5"]),
            "input_template": str(expected_phase_provenance["input_template"]),
            "reader": str(expected_phase_provenance["reader"]),
            "snapshot_indices_json": str(expected_phase_provenance["snapshot_indices_json"]),
            "velocity_path": str(expected_phase_provenance["velocity_path"]),
            "displacement_path": str(expected_phase_provenance["displacement_path"]),
            "integrator": str(expected_phase_provenance["integrator"]),
            "idw_neighbors": int(expected_phase_provenance["idw_neighbors"]),
        },
    )
    metadata_path = args.metadata or args.output_xdmf.with_name(args.output_xdmf.stem + "_metadata.json")
    write_json(
        metadata_path,
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "definition": "arithmetic average of common-mesh Lagrangian starting-phase fields",
            "phase_indices": indices,
            "sources": sources,
            "mesh_fingerprints": expected_fingerprints,
            "phase_provenance": expected_phase_provenance,
            "phase_diagnostics": phase_diagnostics,
            "minimum_valid_fraction": args.minimum_valid_fraction,
            "maximum_nearest_distance_mm": args.maximum_nearest_distance_mm,
            "allow_partial_phase_average": args.allow_partial_phase_average,
            "velocity_path": args.velocity_path,
            "field_name": args.field_name,
            "output_h5": str(args.output_h5),
            "output_xdmf": str(args.output_xdmf),
        },
    )
    print(f"wrote {args.output_xdmf}", flush=True)
    print(f"wrote {metadata_path}", flush=True)


if __name__ == "__main__":
    main()
