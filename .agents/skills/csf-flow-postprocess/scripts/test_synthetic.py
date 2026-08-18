#!/usr/bin/env python3
"""Synthetic constant-field identities for CSF planar integration."""

from __future__ import annotations
import numpy as np
import pyvista as pv
from slice_common import integrate_surface, max_vector_magnitude_on_surface


def main():
    points = np.array([[0, 0, 0], [2, 0, 0], [2, 3, 0], [0, 3, 0]], dtype=float)
    surface = pv.PolyData(points, np.array([4, 0, 1, 2, 3]))
    surface.point_data["Velocity"] = np.tile([1.0, -2.0, 4.0], (4, 1)); surface.point_data["Pressure"] = np.full(4, 7.5); surface.point_data["Displacement"] = np.tile([0.3, 0.4, 0.0], (4, 1))
    result = integrate_surface(surface, ["Pressure"], "Velocity")
    assert np.isclose(result["area"], 6.0)
    assert np.isclose(result["integral:Velocity[2]"], 24.0), "Q must equal v_z A"
    assert np.isclose(result["integral:Pressure"] / result["area"], 7.5)
    assert np.isclose(max_vector_magnitude_on_surface(surface, "Displacement"), 0.5)
    areas = np.array([[6.0, 9.0], [6.0, 12.0]]); relative = (areas - areas[0]) / areas[0]
    assert np.allclose(relative[0], 0.0); assert np.allclose(relative[1], [0.0, 1 / 3])
    rigid = np.repeat(areas[:1], 2, axis=0); assert np.allclose((rigid - rigid[0]) / rigid[0], 0.0)
    print("CSF synthetic identities passed")


if __name__ == "__main__": main()
