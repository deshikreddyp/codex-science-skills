#!/usr/bin/env python3
"""Shared readers and exact planar-integration helpers for CSF post-processing."""

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
        raise ValueError(f"Cannot parse a result index from {path}")
    return int(match.group(1))


def discover_xdmf_series(
    results_dir: Path,
    first: int | None = None,
    last: int | None = None,
    stride: int = 1,
) -> list[Path]:
    if stride <= 0:
        raise ValueError("stride must be positive")
    paths = sorted(
        (path for path in results_dir.iterdir() if RESULT_RE.fullmatch(path.name)),
        key=result_index,
    )
    paths = [
        path
        for path in paths
        if (first is None or result_index(path) >= first)
        and (last is None or result_index(path) <= last)
    ]
    if paths and stride > 1:
        origin = result_index(paths[0])
        paths = [path for path in paths if (result_index(path) - origin) % stride == 0]
    return paths


def parse_indices(spec: str | None, available: Iterable[int]) -> list[int]:
    values = sorted(int(value) for value in available)
    if spec is None:
        return values
    selected: list[int] = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            selected.append(int(token))
            continue
        parts = token.split(":")
        if len(parts) not in (2, 3):
            raise ValueError(f"Invalid index token: {token}")
        start = int(parts[0]) if parts[0] else values[0]
        stop = int(parts[1]) if parts[1] else values[-1] + 1
        step = int(parts[2]) if len(parts) == 3 and parts[2] else 1
        selected.extend(range(start, stop, step))
    missing = sorted(set(selected) - set(values))
    if missing:
        raise ValueError(f"Requested unavailable indices: {missing[:20]}")
    return sorted(dict.fromkeys(selected))


def parse_xdmf(path: Path) -> dict:
    root = ET.parse(path).getroot()
    time_node = root.find(".//Time")
    time_value = None if time_node is None else float(time_node.attrib["Value"])
    attributes = {}
    for attr in root.findall(".//Attribute"):
        item = attr.find("DataItem")
        attributes[attr.attrib.get("Name", "")] = {
            "center": attr.attrib.get("Center"),
            "type": attr.attrib.get("AttributeType"),
            "dimensions": None if item is None else item.attrib.get("Dimensions"),
            "data_item": None if item is None or item.text is None else item.text.strip(),
        }
    topology = root.find(".//Topology")
    geometry = root.find(".//Geometry")
    return {
        "path": str(path),
        "result_index": result_index(path),
        "time_s": time_value,
        "attributes": attributes,
        "topology_type": None if topology is None else topology.attrib.get("TopologyType"),
        "number_of_cells": None
        if topology is None
        else int(topology.attrib.get("NumberOfElements", "0")),
        "geometry_type": None if geometry is None else geometry.attrib.get("GeometryType"),
    }


def read_xdmf_grid(path: Path) -> pv.UnstructuredGrid:
    reader = vtkXdmfReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = pv.wrap(reader.GetOutput())
    if grid is None or grid.n_points == 0:
        raise RuntimeError(f"XDMF reader returned an empty grid for {path}")
    return grid


def numeric_h5_groups(path: Path, group: str) -> list[int]:
    clean = group.strip("/")
    with h5py.File(path, "r") as h5:
        if clean not in h5:
            return []
        return sorted(int(key) for key in h5[clean].keys() if str(key).isdigit())


def inspect_h5(path: Path, max_datasets: int = 100) -> dict:
    groups: dict[str, dict] = {}
    dataset_count = 0
    with h5py.File(path, "r") as h5:
        def visitor(name: str, obj) -> None:
            nonlocal dataset_count
            if isinstance(obj, h5py.Dataset):
                dataset_count += 1
            if isinstance(obj, h5py.Dataset) and len(groups) < max_datasets:
                groups[name] = {"shape": list(obj.shape), "dtype": str(obj.dtype)}

        h5.visititems(visitor)
    return {"path": str(path), "dataset_count": dataset_count, "datasets_shown": groups, "datasets_truncated": dataset_count > len(groups)}


def threshold_domain(
    grid: pv.DataSet,
    field_name: str,
    value: float,
    band: float = 0.5,
    *,
    required: bool,
) -> pv.UnstructuredGrid | pv.PolyData | None:
    if field_name in grid.cell_data:
        preference = "cell"
    elif field_name in grid.point_data:
        preference = "point"
    elif required:
        raise KeyError(f"Missing domain field {field_name!r}")
    else:
        return None
    selected = grid.threshold(
        scalars=field_name,
        value=(value - band, value + band),
        preference=preference,
    )
    if selected.n_cells == 0 and not required:
        return None
    if selected.n_cells == 0:
        raise RuntimeError(f"No cells selected for {field_name}={value}")
    return selected


def warp_if_available(
    grid: pv.DataSet, displacement_name: str, *, required: bool = False
) -> tuple[pv.DataSet, bool]:
    if displacement_name not in grid.point_data:
        if required:
            raise KeyError(f"Missing point vector {displacement_name!r}")
        return grid.copy(deep=False), False
    return grid.warp_by_vector(displacement_name, factor=1.0), True


def flow_slice_grid(
    z_top: float,
    z_bottom: float,
    spacing_mm: float,
    max_slices: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if spacing_mm <= 0.0:
        raise ValueError("slice spacing must be positive")
    count = int(math.floor((z_top - z_bottom) / spacing_mm))
    if max_slices is not None:
        count = min(count, int(max_slices))
    if count <= 0:
        raise RuntimeError(f"Invalid z extent {z_bottom:g}..{z_top:g}")
    distance = (np.arange(count, dtype=float) + 0.5) * spacing_mm
    edges = np.arange(count + 1, dtype=float) * spacing_mm
    planes = z_top - distance
    return distance, edges, planes


def integrate_surface(
    surface: pv.DataSet,
    scalar_names: Iterable[str] = (),
    vector_name: str | None = None,
    axial_component: int = 2,
) -> dict[str, float]:
    result = {"area": math.nan}
    for name in scalar_names:
        result[f"integral:{name}"] = math.nan
    if vector_name is not None:
        result[f"integral:{vector_name}[{axial_component}]"] = math.nan
    if surface.n_cells == 0 or surface.n_points == 0:
        return result
    tri = surface.triangulate()
    if tri.n_cells == 0 or tri.n_points == 0:
        return result
    faces = np.asarray(tri.faces).reshape(-1, 4)[:, 1:]
    points = np.asarray(tri.points, dtype=float)
    p0, p1, p2 = points[faces[:, 0]], points[faces[:, 1]], points[faces[:, 2]]
    areas = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
    result["area"] = float(np.sum(areas))
    for name in scalar_names:
        if name not in tri.point_data:
            continue
        values = np.asarray(tri.point_data[name], dtype=float).reshape(-1)
        result[f"integral:{name}"] = float(np.sum(areas * values[faces].mean(axis=1)))
    if vector_name is not None and vector_name in tri.point_data:
        values = np.asarray(tri.point_data[vector_name], dtype=float)
        axial = values[:, axial_component]
        result[f"integral:{vector_name}[{axial_component}]"] = float(
            np.sum(areas * axial[faces].mean(axis=1))
        )
    return result


def max_vector_magnitude_on_surface(surface: pv.DataSet, vector_name: str) -> float:
    if surface.n_points == 0 or vector_name not in surface.point_data:
        return math.nan
    vectors = np.asarray(surface.point_data[vector_name], dtype=float)
    if vectors.ndim != 2:
        return math.nan
    return float(np.max(np.linalg.norm(vectors, axis=1)))


def mesh_grid_from_h5(mesh_h5: Path) -> tuple[pv.UnstructuredGrid, np.ndarray, np.ndarray, np.ndarray | None]:
    with h5py.File(mesh_h5, "r") as h5:
        coordinates = np.asarray(h5["mesh/coordinates"], dtype=float)
        topology = np.asarray(h5["mesh/topology"], dtype=np.int64)
        domains = np.asarray(h5["domains/values"]) if "domains/values" in h5 else None
    cells = np.empty((topology.shape[0], 5), dtype=np.int64)
    cells[:, 0] = 4
    cells[:, 1:] = topology
    celltypes = np.full(topology.shape[0], pv.CellType.TETRA, dtype=np.uint8)
    grid = pv.UnstructuredGrid(cells.reshape(-1), celltypes, coordinates)
    if domains is not None:
        grid.cell_data["Sub-domain"] = domains.reshape(-1)
    return grid, coordinates, topology, domains


class FenicsSeriesReader:
    """Read CG1 functions through DOLFIN and expose vertex-ordered NumPy arrays."""

    def __init__(self, mesh_h5: Path, snapshots_h5: Path):
        from dolfin import HDF5File, Function, FunctionSpace, Mesh, VectorFunctionSpace, vertex_to_dof_map

        self._HDF5File = HDF5File
        self._Function = Function
        self.mesh = Mesh()
        with HDF5File(self.mesh.mpi_comm(), str(mesh_h5), "r") as hdf:
            hdf.read(self.mesh, "/mesh", False)
        self.snapshots_h5 = snapshots_h5
        self.vector_space = VectorFunctionSpace(self.mesh, "CG", 1)
        self.scalar_space = FunctionSpace(self.mesh, "CG", 1)
        self.vector_dofs = vertex_to_dof_map(self.vector_space).reshape(self.mesh.num_vertices(), 3)
        self.scalar_dofs = vertex_to_dof_map(self.scalar_space)

    def read_vector(self, path: str) -> np.ndarray:
        function = self._Function(self.vector_space)
        with self._HDF5File(self.mesh.mpi_comm(), str(self.snapshots_h5), "r") as hdf:
            hdf.read(function, path)
        return np.asarray(function.vector().get_local()[self.vector_dofs], dtype=float)

    def read_scalar(self, path: str) -> np.ndarray:
        function = self._Function(self.scalar_space)
        with self._HDF5File(self.mesh.mpi_comm(), str(self.snapshots_h5), "r") as hdf:
            hdf.read(function, path)
        return np.asarray(function.vector().get_local()[self.scalar_dofs], dtype=float)


def source_signature(paths: Iterable[Path], parameters: dict) -> str:
    digest = hashlib.sha256()
    expanded = set(Path(item).resolve() for item in paths)
    expanded.update(path.with_suffix(".h5") for path in list(expanded) if path.suffix == ".xdmf" and path.with_suffix(".h5").exists())
    for path in sorted(expanded):
        stat = path.stat()
        digest.update(str(path).encode())
        digest.update(str(stat.st_size).encode())
        digest.update(str(stat.st_mtime_ns).encode())
    digest.update(json.dumps(parameters, sort_keys=True, default=str).encode())
    return digest.hexdigest()


def parse_chunk(spec: str | None, count: int) -> np.ndarray:
    if spec is None:
        return np.arange(count, dtype=int)
    match = re.fullmatch(r"(\d+)/(\d+)", spec)
    if not match:
        raise ValueError("snapshot chunk must use K/N")
    index, total = int(match.group(1)), int(match.group(2))
    if total <= 0 or index < 0 or index >= total:
        raise ValueError("snapshot chunk requires 0 <= K < N")
    return np.arange(count, dtype=int)[index::total]


def ensure_output_dir(path: Path, overwrite: bool) -> None:
    """Create an output directory and guard existing deliverables."""
    deliverables = (
        "csf_slice_metrics.csv",
        "csf_slice_metrics.npz",
        "csf_slice_metrics_metadata.json",
        "selected_snapshots.csv",
    )
    if path.exists() and not overwrite and any((path / name).exists() for name in deliverables):
        raise FileExistsError(
            f"Output deliverables already exist in {path}; pass --overwrite to replace them"
        )
    path.mkdir(parents=True, exist_ok=True)


def slice_at_z(grid: pv.DataSet, z_mm: float) -> pv.PolyData:
    return grid.slice(normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, float(z_mm)))
