#!/usr/bin/env python3
"""Exercise compact FEniCS-HDF5 reconstruction through DOLFIN."""

from __future__ import annotations
import json, subprocess, sys, tempfile
from pathlib import Path
import numpy as np


def main():
    try:
        from dolfin import Constant, Function, FunctionSpace, HDF5File, MeshFunction, UnitCubeMesh, VectorFunctionSpace, cells, interpolate
    except ImportError:
        print("SKIP: dolfin is unavailable")
        return
    script = Path(__file__).with_name("postprocess_csf_slices.py")
    with tempfile.TemporaryDirectory() as raw:
        root = Path(raw); mesh_path = root / "mesh.h5"; snapshots = root / "snapshots.h5"; output = root / "output"
        mesh = UnitCubeMesh(2, 2, 2); domains = MeshFunction("size_t", mesh, mesh.topology().dim(), 258514)
        for cell in cells(mesh):
            if cell.midpoint().x() >= 0.5: domains[cell] = 258515
        with HDF5File(mesh.mpi_comm(), str(mesh_path), "w") as hdf: hdf.write(mesh, "/mesh"); hdf.write(domains, "/domains")
        vector = VectorFunctionSpace(mesh, "CG", 1); scalar = FunctionSpace(mesh, "CG", 1)
        velocity = interpolate(Constant((0.0, 0.0, 4.0)), vector); pressure = interpolate(Constant(7.5), scalar); displacement = interpolate(Constant((0.0, 0.0, 0.0)), vector)
        with HDF5File(mesh.mpi_comm(), str(snapshots), "w") as hdf:
            for index in (0, 1): hdf.write(velocity, f"/velocity/{index}"); hdf.write(pressure, f"/pressure/{index}"); hdf.write(displacement, f"/displacement/{index}")
        subprocess.run([sys.executable, str(script), "--snapshots-h5", str(snapshots), "--mesh-h5", str(mesh_path), "--dt-s", "0.1", "--output-dir", str(output), "--max-slices", "1", "--no-plots"], check=True)
        data = np.load(output / "csf_slice_metrics.npz"); assert data["area_mm2"].shape == (2,1); assert np.allclose(data["flow_rate_mm3_per_s"], 4 * data["area_mm2"]); assert np.allclose(data["mean_pressure"], 7.5); assert np.allclose(data["relative_area_deformation"], 0)
        metadata = json.loads((output / "csf_slice_metrics_metadata.json").read_text()); assert metadata["reader"] == "fenics-h5"
    print("CSF FEniCS-HDF5 end-to-end test passed")


if __name__ == "__main__": main()
