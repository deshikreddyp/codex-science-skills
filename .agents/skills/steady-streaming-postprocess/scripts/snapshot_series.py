#!/usr/bin/env python3
"""Readers for nodal vector fields stored as HDF5 arrays or FEniCS Functions."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np

from steady_streaming_common import format_snapshot, h5_key


class SnapshotReader:
    def read_vector(self, index: int, dataset_template: str) -> np.ndarray:
        raise NotImplementedError


class H5ArrayReader(SnapshotReader):
    def __init__(self, input_template: str):
        self.input_template = input_template

    def read_vector(self, index: int, dataset_template: str) -> np.ndarray:
        source = Path(format_snapshot(self.input_template, index))
        dataset = h5_key(format_snapshot(dataset_template, index))
        with h5py.File(source, "r") as h5:
            if dataset not in h5:
                raise KeyError(f"Dataset /{dataset} not found in {source}")
            values = np.asarray(h5[dataset][:], dtype=np.float64)
        if values.ndim != 2 or values.shape[1] != 3:
            raise ValueError(f"Expected vector dataset shape (N, 3), got {values.shape} at {source}:/{dataset}")
        return values


class FenicsFunctionReader(SnapshotReader):
    def __init__(self, mesh_h5: Path, input_template: str):
        try:
            from dolfin import Function, HDF5File, Mesh, VectorFunctionSpace, set_log_active, vertex_to_dof_map
        except ImportError as exc:
            raise RuntimeError("The fenics-function reader requires legacy dolfin") from exc

        set_log_active(False)
        self._Function = Function
        self._HDF5File = HDF5File
        self.input_template = input_template
        self.mesh = Mesh()
        with HDF5File(self.mesh.mpi_comm(), str(mesh_h5), "r") as hdf:
            hdf.read(self.mesh, "/mesh", False)
        self.space = VectorFunctionSpace(self.mesh, "CG", 1)
        self.vertex_dofs = vertex_to_dof_map(self.space).reshape(self.mesh.num_vertices(), 3)

    def read_vector(self, index: int, dataset_template: str) -> np.ndarray:
        source = format_snapshot(self.input_template, index)
        dataset = format_snapshot(dataset_template, index)
        function = self._Function(self.space)
        with self._HDF5File(self.mesh.mpi_comm(), source, "r") as hdf:
            hdf.read(function, dataset)
        local = function.vector().get_local()
        return np.asarray(local[self.vertex_dofs], dtype=np.float64)


def make_reader(kind: str, mesh_h5: Path, input_template: str) -> SnapshotReader:
    if kind == "h5-array":
        return H5ArrayReader(input_template)
    if kind == "fenics-function":
        return FenicsFunctionReader(mesh_h5, input_template)
    raise ValueError(f"Unknown reader kind: {kind}")
