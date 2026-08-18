#!/usr/bin/env python3
"""Readers and exact planar integration shared by drug-transport scripts."""

from __future__ import annotations

import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pyvista as pv
from vtkmodules.vtkIOXdmf2 import vtkXdmfReader

RESULT_RE = re.compile(r"result(\d+)\.xdmf$")


def result_index(path: Path) -> int:
    match = RESULT_RE.fullmatch(path.name)
    if not match:
        raise ValueError(f"Cannot parse result index from {path}")
    return int(match.group(1))


def discover_xdmf_series(results_dir: Path, first=None, last=None, stride=1) -> list[Path]:
    if stride <= 0:
        raise ValueError("stride must be positive")
    paths = sorted((p for p in results_dir.iterdir() if RESULT_RE.fullmatch(p.name)), key=result_index)
    paths = [p for p in paths if (first is None or result_index(p) >= first) and (last is None or result_index(p) <= last)]
    if paths and stride > 1:
        origin = result_index(paths[0]); paths = [p for p in paths if (result_index(p) - origin) % stride == 0]
    return paths


def parse_xdmf(path: Path) -> dict:
    root = ET.parse(path).getroot(); time_node = root.find(".//Time")
    attrs = {}
    for attr in root.findall(".//Attribute"):
        item = attr.find("DataItem")
        attrs[attr.attrib.get("Name", "")] = {"center": attr.attrib.get("Center"), "type": attr.attrib.get("AttributeType"), "dimensions": None if item is None else item.attrib.get("Dimensions"), "data_item": None if item is None or item.text is None else item.text.strip()}
    topology = root.find(".//Topology")
    return {"path": str(path), "result_index": result_index(path), "time_s": None if time_node is None else float(time_node.attrib["Value"]), "attributes": attrs, "topology_type": None if topology is None else topology.attrib.get("TopologyType")}


def read_xdmf_grid(path: Path):
    reader = vtkXdmfReader(); reader.SetFileName(str(path)); reader.Update(); grid = pv.wrap(reader.GetOutput())
    if grid is None or grid.n_points == 0:
        raise RuntimeError(f"XDMF reader returned an empty grid for {path}")
    return grid


def numeric_h5_groups(path: Path, group: str) -> list[int]:
    with h5py.File(path, "r") as h5:
        clean = group.strip("/")
        return [] if clean not in h5 else sorted(int(k) for k in h5[clean] if str(k).isdigit())


def inspect_h5(path: Path, max_datasets: int = 100) -> dict:
    datasets = {}; dataset_count = 0
    with h5py.File(path, "r") as h5:
        def visitor(name, obj):
            nonlocal dataset_count
            if isinstance(obj, h5py.Dataset): dataset_count += 1
            if isinstance(obj, h5py.Dataset) and len(datasets) < max_datasets: datasets[name] = {"shape": list(obj.shape), "dtype": str(obj.dtype)}
        h5.visititems(visitor)
    return {"path": str(path), "dataset_count": dataset_count, "datasets_shown": datasets, "datasets_truncated": dataset_count > len(datasets)}


def parse_indices(spec: str | None, available: Iterable[int]) -> list[int]:
    values = sorted(int(v) for v in available)
    if spec is None:
        return values
    selected = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            selected.append(int(token)); continue
        parts = token.split(":")
        if len(parts) not in (2, 3):
            raise ValueError(f"Invalid index token: {token}")
        start = int(parts[0]) if parts[0] else values[0]; stop = int(parts[1]) if parts[1] else values[-1] + 1
        step = int(parts[2]) if len(parts) == 3 and parts[2] else 1; selected.extend(range(start, stop, step))
    missing = sorted(set(selected) - set(values))
    if missing:
        raise ValueError(f"Requested unavailable indices: {missing[:20]}")
    return sorted(dict.fromkeys(selected))


def threshold_domain(grid, field_name: str, value: float, *, required: bool):
    if field_name in grid.cell_data:
        preference = "cell"
    elif field_name in grid.point_data:
        preference = "point"
    elif required:
        raise KeyError(f"Missing domain field {field_name!r}")
    else:
        return None
    selected = grid.threshold((value - 0.5, value + 0.5), scalars=field_name, preference=preference)
    if selected.n_cells == 0 and not required:
        return None
    if selected.n_cells == 0:
        raise RuntimeError(f"No cells selected for {field_name}={value}")
    return selected


def warp_if_available(grid, displacement_name: str):
    if displacement_name not in grid.point_data:
        return grid.copy(deep=False), False
    return grid.warp_by_vector(displacement_name, factor=1.0), True


def transport_slice_grid(fm_z: float, z_bottom: float, spacing_mm: float, max_slices=None):
    if spacing_mm <= 0:
        raise ValueError("slice spacing must be positive")
    count = int(math.floor((fm_z - z_bottom) / spacing_mm)) + 1
    if max_slices is not None:
        count = min(count, int(max_slices))
    if count <= 0:
        raise RuntimeError(f"Invalid z extent {z_bottom:g}..{fm_z:g}")
    distance = np.arange(count, dtype=float) * spacing_mm
    return distance, fm_z - distance


def integrate_surface(surface, scalar_name: str):
    if surface.n_cells == 0 or surface.n_points == 0:
        return math.nan, math.nan
    tri = surface.triangulate()
    if scalar_name not in tri.point_data:
        raise KeyError(f"Missing point scalar {scalar_name!r} on slice")
    faces = np.asarray(tri.faces).reshape(-1, 4)[:, 1:]
    points = np.asarray(tri.points, dtype=float)
    area_each = 0.5 * np.linalg.norm(np.cross(points[faces[:, 1]] - points[faces[:, 0]], points[faces[:, 2]] - points[faces[:, 0]]), axis=1)
    values = np.asarray(tri.point_data[scalar_name], dtype=float).reshape(-1)
    return float(np.sum(area_each)), float(np.sum(area_each * values[faces].mean(axis=1)))


def mesh_grid_from_h5(mesh_h5: Path):
    with h5py.File(mesh_h5, "r") as h5:
        coordinates = np.asarray(h5["mesh/coordinates"], dtype=float); topology = np.asarray(h5["mesh/topology"], dtype=np.int64)
        domains = np.asarray(h5["domains/values"]).reshape(-1) if "domains/values" in h5 else None
    cells = np.empty((topology.shape[0], 5), dtype=np.int64); cells[:, 0] = 4; cells[:, 1:] = topology
    grid = pv.UnstructuredGrid(cells.reshape(-1), np.full(topology.shape[0], pv.CellType.TETRA, dtype=np.uint8), coordinates)
    if domains is not None:
        grid.cell_data["Sub-domain"] = domains
    return grid, domains


class FenicsSeriesReader:
    """Reconstruct CG1 functions with DOLFIN; never assume vector_0 is nodal."""
    def __init__(self, mesh_h5: Path, snapshots_h5: Path):
        from dolfin import HDF5File, Function, FunctionSpace, Mesh, VectorFunctionSpace, vertex_to_dof_map
        self._HDF5File = HDF5File; self._Function = Function; self.snapshots_h5 = snapshots_h5; self.mesh = Mesh()
        with HDF5File(self.mesh.mpi_comm(), str(mesh_h5), "r") as hdf:
            hdf.read(self.mesh, "/mesh", False)
        self.scalar_space = FunctionSpace(self.mesh, "CG", 1); self.vector_space = VectorFunctionSpace(self.mesh, "CG", 1)
        self.scalar_dofs = vertex_to_dof_map(self.scalar_space); self.vector_dofs = vertex_to_dof_map(self.vector_space).reshape(self.mesh.num_vertices(), 3)
    def read_scalar(self, path: str):
        f = self._Function(self.scalar_space)
        with self._HDF5File(self.mesh.mpi_comm(), str(self.snapshots_h5), "r") as hdf: hdf.read(f, path)
        return np.asarray(f.vector().get_local()[self.scalar_dofs], dtype=float)
    def read_vector(self, path: str):
        f = self._Function(self.vector_space)
        with self._HDF5File(self.mesh.mpi_comm(), str(self.snapshots_h5), "r") as hdf: hdf.read(f, path)
        return np.asarray(f.vector().get_local()[self.vector_dofs], dtype=float)


def source_signature(paths: Iterable[Path], parameters: dict) -> str:
    digest = hashlib.sha256()
    expanded = set(Path(p).resolve() for p in paths)
    expanded.update(path.with_suffix(".h5") for path in list(expanded) if path.suffix == ".xdmf" and path.with_suffix(".h5").exists())
    for path in sorted(expanded):
        stat = path.stat(); digest.update(str(path).encode()); digest.update(str(stat.st_size).encode()); digest.update(str(stat.st_mtime_ns).encode())
    digest.update(json.dumps(parameters, sort_keys=True, default=str).encode()); return digest.hexdigest()


def parse_chunk(spec: str | None, count: int) -> np.ndarray:
    if spec is None:
        return np.arange(count, dtype=int)
    match = re.fullmatch(r"(\d+)/(\d+)", spec)
    if not match:
        raise ValueError("snapshot chunk must use K/N")
    index, total = map(int, match.groups())
    if total <= 0 or index < 0 or index >= total:
        raise ValueError("snapshot chunk requires 0 <= K < N")
    return np.arange(count, dtype=int)[index::total]


def ensure_output_dir(path: Path, overwrite: bool):
    names = ("net_concentration.csv", "avg_concentration.csv", "concentration_slice_metrics.csv", "concentration_slice_metrics.npz", "concentration_slice_metrics_metadata.json")
    if path.exists() and not overwrite and any((path / name).exists() for name in names):
        raise FileExistsError(f"Output deliverables already exist in {path}; pass --overwrite to replace them")
    path.mkdir(parents=True, exist_ok=True)
