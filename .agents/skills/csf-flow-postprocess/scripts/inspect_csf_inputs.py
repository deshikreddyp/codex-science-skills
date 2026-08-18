#!/usr/bin/env python3
"""Inspect CSF XDMF or compact FEniCS HDF5 inputs without modifying them."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from slice_common import discover_xdmf_series, inspect_h5, numeric_h5_groups, parse_xdmf


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--results-dir", type=Path)
    source.add_argument("--snapshots-h5", type=Path)
    parser.add_argument("--mesh-h5", type=Path)
    parser.add_argument("--velocity-group", default="velocity")
    parser.add_argument("--pressure-group", default="pressure")
    parser.add_argument("--displacement-group", default="displacement")
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def inspect_xdmf(results_dir: Path) -> dict:
    paths = discover_xdmf_series(results_dir)
    if not paths:
        raise FileNotFoundError(f"No result*.xdmf files in {results_dir}")
    first = parse_xdmf(paths[0])
    last = parse_xdmf(paths[-1])
    return {
        "reader": "xdmf-series",
        "results_dir": str(results_dir),
        "file_count": len(paths),
        "first": first,
        "last": last,
        "times_monotone": first["time_s"] is None
        or last["time_s"] is None
        or last["time_s"] >= first["time_s"],
    }


def inspect_fenics(args: argparse.Namespace) -> dict:
    if args.mesh_h5 is None:
        raise ValueError("--mesh-h5 is required with --snapshots-h5")
    groups = {
        "velocity": numeric_h5_groups(args.snapshots_h5, args.velocity_group),
        "pressure": numeric_h5_groups(args.snapshots_h5, args.pressure_group),
        "displacement": numeric_h5_groups(args.snapshots_h5, args.displacement_group),
    }
    group_summary = {
        key: {"count": len(value), "first": value[:10], "last": value[-10:]}
        for key, value in groups.items()
    }
    with h5py.File(args.mesh_h5, "r") as h5:
        coordinates = h5["mesh/coordinates"]
        topology = h5["mesh/topology"]
        tags = np.unique(h5["domains/values"][:]) if "domains/values" in h5 else np.array([])
        mesh = {
            "coordinates_shape": list(coordinates.shape),
            "topology_shape": list(topology.shape),
            "domain_values": tags.tolist(),
        }
    return {
        "reader": "fenics-h5",
        "mesh_h5": str(args.mesh_h5),
        "snapshots_h5": str(args.snapshots_h5),
        "mesh": mesh,
        "numeric_groups": group_summary,
        "snapshot_file": inspect_h5(args.snapshots_h5),
    }


def main() -> None:
    args = parse_args()
    report = inspect_xdmf(args.results_dir) if args.results_dir else inspect_fenics(args)
    text = json.dumps(report, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
