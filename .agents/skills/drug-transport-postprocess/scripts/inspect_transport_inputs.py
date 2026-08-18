#!/usr/bin/env python3
"""Inspect drug-transport XDMF or compact FEniCS HDF5 inputs."""

from __future__ import annotations
import argparse, json
from pathlib import Path
import h5py, numpy as np
from slice_common import discover_xdmf_series, inspect_h5, numeric_h5_groups, parse_xdmf


def main():
    p = argparse.ArgumentParser(description=__doc__); src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--results-dir", type=Path); src.add_argument("--snapshots-h5", type=Path)
    p.add_argument("--mesh-h5", type=Path); p.add_argument("--concentration-group", default="concentration"); p.add_argument("--displacement-group", default="displacement"); p.add_argument("--output", type=Path); a = p.parse_args()
    if a.results_dir:
        paths = discover_xdmf_series(a.results_dir)
        if not paths: raise FileNotFoundError(f"No result*.xdmf files in {a.results_dir}")
        report = {"reader": "xdmf-series", "results_dir": str(a.results_dir), "file_count": len(paths), "first": parse_xdmf(paths[0]), "last": parse_xdmf(paths[-1])}
    else:
        if a.mesh_h5 is None: raise ValueError("--mesh-h5 is required with --snapshots-h5")
        with h5py.File(a.mesh_h5, "r") as h5:
            report_mesh = {"coordinates_shape": list(h5["mesh/coordinates"].shape), "topology_shape": list(h5["mesh/topology"].shape), "domain_values": np.unique(h5["domains/values"][:]).tolist() if "domains/values" in h5 else []}
        groups = {"concentration": numeric_h5_groups(a.snapshots_h5, a.concentration_group), "displacement": numeric_h5_groups(a.snapshots_h5, a.displacement_group)}
        summary = {key: {"count": len(value), "first": value[:10], "last": value[-10:]} for key, value in groups.items()}
        report = {"reader": "fenics-h5", "mesh_h5": str(a.mesh_h5), "snapshots_h5": str(a.snapshots_h5), "mesh": report_mesh, "numeric_groups": summary, "snapshot_file": inspect_h5(a.snapshots_h5)}
    text = json.dumps(report, indent=2)
    if a.output: a.output.parent.mkdir(parents=True, exist_ok=True); a.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__": main()
