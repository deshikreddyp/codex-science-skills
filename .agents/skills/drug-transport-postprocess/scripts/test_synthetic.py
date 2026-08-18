#!/usr/bin/env python3
"""Synthetic constant-field identities for transport planar integration."""

from __future__ import annotations
import numpy as np
import pyvista as pv
from slice_common import integrate_surface, transport_slice_grid


def main():
    points = np.array([[0, 0, 0], [2, 0, 0], [2, 3, 0], [0, 3, 0]], dtype=float)
    surface = pv.PolyData(points, np.array([4, 0, 1, 2, 3])); surface.point_data["Concentration"] = np.full(4, 2.25)
    area, cnet = integrate_surface(surface, "Concentration"); cavg = cnet / area
    assert np.isclose(area, 6.0); assert np.isclose(cnet, 13.5), "c_net must equal c A"; assert np.isclose(cavg, 2.25), "c_avg must equal c"
    distance, planes = transport_slice_grid(8.0, 5.0, 1.0); assert np.array_equal(distance, [0, 1, 2, 3]); assert np.array_equal(planes, [8, 7, 6, 5])
    print("Transport synthetic identities passed")


if __name__ == "__main__": main()
