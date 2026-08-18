#!/usr/bin/env python3
"""Shared mesh, HDF5, XDMF, and slice-integration utilities."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

import h5py
import numpy as np


@dataclass(frozen=True)
class MeshData:
    points: np.ndarray
    topology: np.ndarray
    cell_tags: np.ndarray | None


def h5_key(path: str) -> str:
    return path.lstrip("/")


def parse_indices(spec: str) -> list[int]:
    """Parse comma-separated integers and Python-style start:stop[:step] ranges."""
    result: list[int] = []
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            result.append(int(token))
            continue
        parts = token.split(":")
        if len(parts) not in (2, 3) or not parts[0] or not parts[1]:
            raise ValueError(f"Invalid index range: {token!r}")
        start, stop = int(parts[0]), int(parts[1])
        step = int(parts[2]) if len(parts) == 3 and parts[2] else 1
        if step == 0:
            raise ValueError("Index range step cannot be zero")
        result.extend(range(start, stop, step))
    if not result:
        raise ValueError("No snapshot indices were selected")
    if len(set(result)) != len(result):
        raise ValueError("Snapshot indices contain duplicates")
    return result


def format_snapshot(template: str, index: int) -> str:
    return template.format(index=index, i=index)


def first_existing_dataset(h5: h5py.File, candidates: Sequence[str]) -> str:
    for path in candidates:
        if h5_key(path) in h5:
            return path
    raise KeyError(f"None of these HDF5 datasets exists: {', '.join(candidates)}")


def read_mesh(
    mesh_h5: Path,
    geometry_path: str | None = None,
    topology_path: str | None = None,
    cell_tags_path: str | None = "/domains/values",
) -> MeshData:
    with h5py.File(mesh_h5, "r") as h5:
        geometry_path = geometry_path or first_existing_dataset(
            h5, ("/mesh/coordinates", "/mesh/geometry", "/Mesh/0/mesh/geometry")
        )
        topology_path = topology_path or first_existing_dataset(
            h5, ("/mesh/topology", "/Mesh/0/mesh/topology")
        )
        points = np.asarray(h5[h5_key(geometry_path)][:], dtype=np.float64)
        topology = np.asarray(h5[h5_key(topology_path)][:], dtype=np.int64)
        cell_tags = None
        if cell_tags_path and h5_key(cell_tags_path) in h5:
            cell_tags = np.asarray(h5[h5_key(cell_tags_path)][:])
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError(f"Expected geometry shape (N, 3), got {points.shape}")
    if topology.ndim != 2 or topology.shape[1] != 4:
        raise ValueError(f"Expected tetrahedral topology shape (M, 4), got {topology.shape}")
    if cell_tags is not None and len(cell_tags) != len(topology):
        raise ValueError("Cell-tag count does not match topology cell count")
    return MeshData(points=points, topology=topology, cell_tags=cell_tags)


def select_cells(mesh: MeshData, fluid_tag: int | None) -> tuple[np.ndarray, np.ndarray]:
    if fluid_tag is None:
        mask = np.ones(len(mesh.topology), dtype=bool)
    else:
        if mesh.cell_tags is None:
            raise ValueError("A fluid tag was requested, but the mesh has no cell-tag dataset")
        mask = mesh.cell_tags == fluid_tag
        if not np.any(mask):
            raise ValueError(f"No cells have fluid tag {fluid_tag}")
    tags = (
        np.asarray(mesh.cell_tags[mask])
        if mesh.cell_tags is not None
        else np.full(np.count_nonzero(mask), fluid_tag if fluid_tag is not None else 1, dtype=np.int64)
    )
    return mesh.topology[mask], tags


def compact_tetra_data(
    points: np.ndarray,
    topology_global: np.ndarray,
    point_fields: dict[str, np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], np.ndarray]:
    used = np.unique(topology_global.reshape(-1))
    global_to_local = np.full(len(points), -1, dtype=np.int64)
    global_to_local[used] = np.arange(len(used), dtype=np.int64)
    local_topology = global_to_local[topology_global]
    compact_fields = {name: np.asarray(values)[used] for name, values in (point_fields or {}).items()}
    return points[used], local_topology, compact_fields, used


def build_pyvista_grid(
    points: np.ndarray,
    topology: np.ndarray,
    point_fields: dict[str, np.ndarray] | None = None,
    cell_fields: dict[str, np.ndarray] | None = None,
):
    import pyvista as pv

    cells = np.empty((len(topology), 5), dtype=np.int64)
    cells[:, 0] = 4
    cells[:, 1:] = topology
    cell_types = np.full(len(topology), pv.CellType.TETRA, dtype=np.uint8)
    grid = pv.UnstructuredGrid(cells.reshape(-1), cell_types, points)
    for name, values in (point_fields or {}).items():
        grid.point_data[name] = values
    for name, values in (cell_fields or {}).items():
        grid.cell_data[name] = values
    return grid


def sample_vector(grid, locations: np.ndarray, field_name: str) -> tuple[np.ndarray, np.ndarray]:
    import pyvista as pv

    cloud = pv.PolyData(np.asarray(locations, dtype=np.float64))
    sampled = cloud.sample(grid)
    values = np.asarray(sampled.point_data[field_name], dtype=np.float64)
    mask_name = "vtkValidPointMask"
    valid = (
        np.asarray(sampled.point_data[mask_name], dtype=bool)
        if mask_name in sampled.point_data
        else np.all(np.isfinite(values), axis=1)
    )
    return values, valid


def positive_linear_triangle_integral(areas: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Integrate max(linear nodal field, 0) exactly on each triangle."""
    positive = np.zeros(values.shape[0], dtype=np.float64)
    n_positive = np.count_nonzero(values > 0.0, axis=1)

    mask = n_positive == 3
    positive[mask] = areas[mask] * values[mask].mean(axis=1)

    for i, j, k in ((0, 1, 2), (1, 0, 2), (2, 0, 1)):
        mask = (n_positive == 1) & (values[:, i] > 0.0)
        if np.any(mask):
            a = values[mask, i]
            b = values[mask, j]
            c = values[mask, k]
            t_b = a / (a - b)
            t_c = a / (a - c)
            positive[mask] = areas[mask] * t_b * t_c * a / 3.0

    for k, i, j in ((0, 1, 2), (1, 0, 2), (2, 0, 1)):
        mask = (n_positive == 2) & (values[:, k] <= 0.0)
        if np.any(mask):
            c = values[mask, k]
            a = values[mask, i]
            b = values[mask, j]
            total = areas[mask] * (a + b + c) / 3.0
            t_a = (-c) / (a - c)
            t_b = (-c) / (b - c)
            negative_corner = areas[mask] * t_a * t_b * c / 3.0
            positive[mask] = total - negative_corner
    return positive


def integrate_slice(grid, z_value: float, field_name: str, axial_component: int = 2) -> dict[str, float | int]:
    sliced = grid.slice(normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, float(z_value)))
    if sliced.n_cells == 0 or sliced.n_points == 0:
        return {"cranial": np.nan, "caudal": np.nan, "area": np.nan, "triangles": 0}
    tri = sliced.triangulate()
    if tri.n_cells == 0 or tri.n_points == 0:
        return {"cranial": np.nan, "caudal": np.nan, "area": np.nan, "triangles": 0}
    faces = np.asarray(tri.faces).reshape(-1, 4)
    triangles = faces[:, 1:]
    p0 = tri.points[triangles[:, 0]]
    p1 = tri.points[triangles[:, 1]]
    p2 = tri.points[triangles[:, 2]]
    areas = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
    axial = np.asarray(tri.point_data[field_name], dtype=np.float64)[:, axial_component]
    values = axial[triangles]
    cranial = float(np.sum(positive_linear_triangle_integral(areas, values)))
    caudal = float(-np.sum(positive_linear_triangle_integral(areas, -values)))
    return {
        "cranial": cranial,
        "caudal": caudal,
        "area": float(np.sum(areas)),
        "triangles": int(tri.n_cells),
    }


def slice_positions(points: np.ndarray, spacing_mm: float) -> np.ndarray:
    if spacing_mm <= 0.0:
        raise ValueError("Slice spacing must be positive")
    z_min = float(np.min(points[:, 2]))
    z_max = float(np.max(points[:, 2]))
    first = np.ceil(z_min / spacing_mm) * spacing_mm
    last = np.floor(z_max / spacing_mm) * spacing_mm
    return np.arange(first, last + 0.5 * spacing_mm, spacing_mm, dtype=np.float64)


def strength_rows(
    grid,
    field_name: str,
    spacing_mm: float,
    case: str,
    label: str,
    axial_component: int = 2,
    progress: Callable[[str], None] | None = print,
    z_values: np.ndarray | None = None,
    z_upper: float | None = None,
) -> list[dict[str, object]]:
    if z_values is None:
        z_values = slice_positions(np.asarray(grid.points), spacing_mm)
    if z_upper is None:
        z_upper = float(np.max(grid.points[:, 2]))
    rows: list[dict[str, object]] = []
    for count, z_value in enumerate(z_values, start=1):
        section = integrate_slice(grid, float(z_value), field_name, axial_component)
        if int(section["triangles"]) == 0:
            continue
        cranial = float(section["cranial"])
        caudal = float(section["caudal"])
        area = float(section["area"])
        net = cranial + caudal
        rows.append(
            {
                "case": case,
                "label": label,
                "z_mm": float(z_value),
                "distance_from_fm_mm": z_upper - float(z_value),
                "area_mm2": area,
                "cranial_strength_mm3_s": cranial,
                "caudal_strength_mm3_s": caudal,
                "net_strength_mm3_s": net,
                "cranial_strength_mL_s": cranial / 1000.0,
                "caudal_strength_mL_s": caudal / 1000.0,
                "net_strength_mL_s": net / 1000.0,
                "cranial_mean_velocity_mm_s": cranial / area if area > 0.0 else np.nan,
                "caudal_mean_velocity_mm_s": caudal / area if area > 0.0 else np.nan,
                "net_mean_velocity_mm_s": net / area if area > 0.0 else np.nan,
                "slice_triangles": int(section["triangles"]),
            }
        )
        if progress and (count == 1 or count % 50 == 0 or count == len(z_values)):
            progress(
                f"slice {count}/{len(z_values)}: d={z_upper-z_value:.1f} mm, "
                f"Q+={cranial/1000.0:.6g} mL/s, Q-={caudal/1000.0:.6g} mL/s"
            )
    return rows


STRENGTH_FIELDS = [
    "case",
    "label",
    "z_mm",
    "distance_from_fm_mm",
    "area_mm2",
    "cranial_strength_mm3_s",
    "caudal_strength_mm3_s",
    "net_strength_mm3_s",
    "cranial_strength_mL_s",
    "caudal_strength_mL_s",
    "net_strength_mL_s",
    "cranial_mean_velocity_mm_s",
    "caudal_mean_velocity_mm_s",
    "net_mean_velocity_mm_s",
    "slice_triangles",
]


def write_csv(path: Path, rows: Iterable[dict[str, object]], fields: Sequence[str] = STRENGTH_FIELDS) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def write_field_h5_xdmf(
    h5_path: Path,
    xdmf_path: Path,
    points: np.ndarray,
    topology: np.ndarray,
    velocity: np.ndarray,
    field_name: str,
    cell_tags: np.ndarray | None = None,
    h5_attributes: dict[str, object] | None = None,
) -> None:
    h5_path.parent.mkdir(parents=True, exist_ok=True)
    xdmf_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(h5_path, "w") as h5:
        h5.create_dataset("Mesh/0/mesh/geometry", data=np.asarray(points, dtype=np.float64))
        h5.create_dataset("Mesh/0/mesh/topology", data=np.asarray(topology, dtype=np.int64))
        h5.create_dataset("VisualisationVector/0", data=np.asarray(velocity, dtype=np.float64))
        if cell_tags is not None:
            h5.create_dataset("Mesh/0/mesh/subdomain", data=np.asarray(cell_tags, dtype=np.int64))
        h5.attrs["field_name"] = field_name
        h5.attrs["geometry_sha256"] = array_sha256(np.asarray(points, dtype=np.float64))
        h5.attrs["topology_sha256"] = array_sha256(np.asarray(topology, dtype=np.int64))
        if cell_tags is not None:
            h5.attrs["cell_tags_sha256"] = array_sha256(np.asarray(cell_tags, dtype=np.int64))
        for name, value in (h5_attributes or {}).items():
            h5.attrs[name] = value

    h5_reference = os.path.relpath(h5_path.resolve(), start=xdmf_path.parent.resolve())
    xdmf = ET.Element("Xdmf", Version="3.0")
    domain = ET.SubElement(xdmf, "Domain")
    grid = ET.SubElement(domain, "Grid", Name="steady_streaming", GridType="Uniform")
    topology_xml = ET.SubElement(
        grid, "Topology", TopologyType="Tetrahedron", NumberOfElements=str(len(topology))
    )
    topology_data = ET.SubElement(
        topology_xml,
        "DataItem",
        Dimensions=f"{len(topology)} 4",
        NumberType="Int",
        Format="HDF",
    )
    topology_data.text = f"{h5_reference}:/Mesh/0/mesh/topology"
    geometry_xml = ET.SubElement(grid, "Geometry", GeometryType="XYZ")
    geometry_data = ET.SubElement(
        geometry_xml,
        "DataItem",
        Dimensions=f"{len(points)} 3",
        NumberType="Float",
        Precision="8",
        Format="HDF",
    )
    geometry_data.text = f"{h5_reference}:/Mesh/0/mesh/geometry"
    attribute = ET.SubElement(
        grid, "Attribute", Name=field_name, AttributeType="Vector", Center="Node"
    )
    field_data = ET.SubElement(
        attribute,
        "DataItem",
        Dimensions=f"{len(velocity)} 3",
        NumberType="Float",
        Precision="8",
        Format="HDF",
    )
    field_data.text = f"{h5_reference}:/VisualisationVector/0"
    if cell_tags is not None:
        tag_attribute = ET.SubElement(
            grid, "Attribute", Name="subdomain", AttributeType="Scalar", Center="Cell"
        )
        tag_data = ET.SubElement(
            tag_attribute,
            "DataItem",
            Dimensions=str(len(cell_tags)),
            NumberType="Int",
            Format="HDF",
        )
        tag_data.text = f"{h5_reference}:/Mesh/0/mesh/subdomain"
    ET.indent(xdmf, space="  ")
    ET.ElementTree(xdmf).write(xdmf_path, encoding="utf-8", xml_declaration=True)


def array_sha256(values: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(values)
    return hashlib.sha256(memoryview(contiguous).cast("B")).hexdigest()
