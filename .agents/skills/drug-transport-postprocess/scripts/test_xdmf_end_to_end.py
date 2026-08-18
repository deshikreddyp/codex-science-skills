#!/usr/bin/env python3
"""Exercise the transport XDMF reader, exporters, and default plotter."""

from __future__ import annotations
import csv, json, os, subprocess, sys, tempfile
from pathlib import Path
import h5py, numpy as np


def write_snapshot(folder: Path, index: int, time_s: float):
    points = np.array([[0,0,-0.2],[2,0,-0.2],[0,2,-0.2],[0,0,2.2]], float); cells = np.array([[0,1,2,3]], np.uint32); concentration = np.full((4,1), 2.25); displacement = np.zeros((4,3))
    with h5py.File(folder / f"result{index}.h5", "w") as h5:
        h5["Mesh/0/mesh/topology"] = cells; h5["Mesh/0/mesh/geometry"] = points; h5["VisualisationVector/0"] = concentration; h5["VisualisationVector/1"] = displacement
    xml = f"""<?xml version="1.0"?><Xdmf Version="3.0"><Domain><Grid Name="mesh" GridType="Uniform"><Topology NumberOfElements="1" TopologyType="Tetrahedron"><DataItem Dimensions="1 4" NumberType="UInt" Format="HDF">result{index}.h5:/Mesh/0/mesh/topology</DataItem></Topology><Geometry GeometryType="XYZ"><DataItem Dimensions="4 3" Format="HDF">result{index}.h5:/Mesh/0/mesh/geometry</DataItem></Geometry><Time Value="{time_s}"/><Attribute Name="Concentration" AttributeType="Scalar" Center="Node"><DataItem Dimensions="4 1" Format="HDF">result{index}.h5:/VisualisationVector/0</DataItem></Attribute><Attribute Name="Displacement" AttributeType="Vector" Center="Node"><DataItem Dimensions="4 3" Format="HDF">result{index}.h5:/VisualisationVector/1</DataItem></Attribute></Grid></Domain></Xdmf>"""
    (folder / f"result{index}.xdmf").write_text(xml, encoding="utf-8")


def main():
    script = Path(__file__).with_name("postprocess_concentration_slices.py")
    merger = Path(__file__).with_name("merge_transport_partials.py")
    with tempfile.TemporaryDirectory() as raw:
        root = Path(raw); results = root / "results"; output = root / "output"; results.mkdir(); write_snapshot(results, 0, 0.0); write_snapshot(results, 1, 0.5)
        env = os.environ.copy(); env["MPLCONFIGDIR"] = str(root / "mpl")
        subprocess.run([sys.executable, str(script), "--results-dir", str(results), "--output-dir", str(output), "--max-slices", "3"], check=True, env=env)
        required = ["net_concentration.csv", "avg_concentration.csv", "concentration_slice_metrics.csv", "concentration_slice_metrics.npz", "selected_timesteps.csv", "concentration_slice_metrics_metadata.json", "c_avg_heatmap.png", "c_avg_heatmap.pdf", "c_net_heatmap.png", "c_net_heatmap.pdf"]
        assert all((output / name).is_file() for name in required)
        data = np.load(output / "concentration_slice_metrics.npz"); valid = np.isfinite(data["area_mm2"]); assert np.allclose(data["c_net"][valid], 2.25 * data["area_mm2"][valid]); assert np.allclose(data["c_avg"][valid], 2.25)
        metadata = json.loads((output / "concentration_slice_metrics_metadata.json").read_text()); assert "not global 3D mass" in metadata["definitions"]["c_net"]
        with (output / "concentration_slice_metrics.csv").open(newline="") as stream: assert len(list(csv.DictReader(stream))) == 6
        partials = root / "partials"; merged = root / "merged"
        for chunk in ("0/2", "1/2"):
            subprocess.run([sys.executable, str(script), "--results-dir", str(results), "--output-dir", str(partials), "--max-slices", "3", "--snapshot-chunk", chunk, "--no-plots"], check=True, env=env)
        subprocess.run([sys.executable, str(merger), str(partials), "--output-dir", str(merged), "--no-plots"], check=True, env=env)
        merged_data = np.load(merged / "concentration_slice_metrics.npz")
        for key in ("times_s", "distance_from_FM_mm", "area_mm2", "c_net", "c_avg"): assert np.allclose(merged_data[key], data[key], equal_nan=True)
    print("Transport XDMF end-to-end test passed")


if __name__ == "__main__": main()
