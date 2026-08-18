#!/usr/bin/env python3
"""Exercise compact concentration/displacement HDF5 reconstruction through DOLFIN."""

from __future__ import annotations
import json, subprocess, sys, tempfile
from pathlib import Path
import numpy as np


def main():
    try:
        from dolfin import Constant, FunctionSpace, HDF5File, MeshFunction, UnitCubeMesh, VectorFunctionSpace, cells, interpolate
    except ImportError:
        print("SKIP: dolfin is unavailable")
        return
    script = Path(__file__).with_name("postprocess_concentration_slices.py")
    with tempfile.TemporaryDirectory() as raw:
        root = Path(raw); mesh_path = root / "mesh.h5"; snapshots = root / "snapshots.h5"; output = root / "output"
        mesh = UnitCubeMesh(2, 2, 2); domains = MeshFunction("size_t", mesh, mesh.topology().dim(), 258514)
        with HDF5File(mesh.mpi_comm(), str(mesh_path), "w") as hdf: hdf.write(mesh, "/mesh"); hdf.write(domains, "/domains")
        scalar = FunctionSpace(mesh, "CG", 1); vector = VectorFunctionSpace(mesh, "CG", 1); concentration = interpolate(Constant(2.25), scalar); displacement = interpolate(Constant((0.0, 0.0, 0.0)), vector)
        with HDF5File(mesh.mpi_comm(), str(snapshots), "w") as hdf:
            for index in (0, 1): hdf.write(concentration, f"/concentration/{index}"); hdf.write(displacement, f"/displacement/{index}")
        subprocess.run([sys.executable, str(script), "--snapshots-h5", str(snapshots), "--mesh-h5", str(mesh_path), "--dt-s", "0.1", "--output-dir", str(output), "--fm-z-mm", "0.5", "--z-bottom-mm", "0.5", "--no-plots"], check=True)
        data = np.load(output / "concentration_slice_metrics.npz"); valid = np.isfinite(data["area_mm2"]); assert np.allclose(data["c_net"][valid], 2.25 * data["area_mm2"][valid]); assert np.allclose(data["c_avg"][valid], 2.25)
        metadata = json.loads((output / "concentration_slice_metrics_metadata.json").read_text()); assert metadata["reader"] == "fenics-h5"
    print("Transport FEniCS-HDF5 end-to-end test passed")


if __name__ == "__main__": main()
