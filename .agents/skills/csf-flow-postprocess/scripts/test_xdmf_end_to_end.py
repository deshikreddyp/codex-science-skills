#!/usr/bin/env python3
"""Exercise the XDMF reader, FSI/domain logic, exporters, and plotter."""

from __future__ import annotations
import csv, json, os, subprocess, sys, tempfile
from pathlib import Path
import h5py, numpy as np


def write_snapshot(folder: Path, index: int, time_s: float):
    points = np.array([[0,0,0],[2,0,0],[0,2,0],[0,0,2],[4,0,0],[6,0,0],[4,2,0],[4,0,2]], float)
    cells = np.array([[0,1,2,3],[4,5,6,7]], np.uint32); velocity = np.tile([0.,0.,4.], (8,1)); pressure = np.full((8,1), 7.5); displacement = np.zeros((8,3)); displacement[4:] = [0.3,0.4,0.]
    h5_path = folder / f"result{index}.h5"
    with h5py.File(h5_path, "w") as h5:
        h5["Mesh/0/mesh/topology"] = cells; h5["Mesh/0/mesh/geometry"] = points; h5["VisualisationVector/0"] = velocity; h5["VisualisationVector/1"] = pressure; h5["VisualisationVector/2"] = displacement; h5["VisualisationVector/3"] = np.array([[0],[1]], np.int32)
    attributes = """<Attribute Name="Velocity" AttributeType="Vector" Center="Node"><DataItem Dimensions="8 3" Format="HDF">resultINDEX.h5:/VisualisationVector/0</DataItem></Attribute>
<Attribute Name="Pressure" AttributeType="Scalar" Center="Node"><DataItem Dimensions="8 1" Format="HDF">resultINDEX.h5:/VisualisationVector/1</DataItem></Attribute>
<Attribute Name="Displacement" AttributeType="Vector" Center="Node"><DataItem Dimensions="8 3" Format="HDF">resultINDEX.h5:/VisualisationVector/2</DataItem></Attribute>
<Attribute Name="Sub-domain" AttributeType="Scalar" Center="Cell"><DataItem Dimensions="2 1" Format="HDF">resultINDEX.h5:/VisualisationVector/3</DataItem></Attribute>""".replace("INDEX", str(index))
    xml = f"""<?xml version="1.0"?><Xdmf Version="3.0"><Domain><Grid Name="mesh" GridType="Uniform"><Topology NumberOfElements="2" TopologyType="Tetrahedron"><DataItem Dimensions="2 4" NumberType="UInt" Format="HDF">result{index}.h5:/Mesh/0/mesh/topology</DataItem></Topology><Geometry GeometryType="XYZ"><DataItem Dimensions="8 3" Format="HDF">result{index}.h5:/Mesh/0/mesh/geometry</DataItem></Geometry><Time Value="{time_s}"/>{attributes}</Grid></Domain></Xdmf>"""
    (folder / f"result{index}.xdmf").write_text(xml, encoding="utf-8")


def main():
    script = Path(__file__).with_name("postprocess_csf_slices.py")
    merger = Path(__file__).with_name("merge_csf_partials.py")
    with tempfile.TemporaryDirectory() as raw:
        root = Path(raw); results = root / "results"; output = root / "output"; results.mkdir(); write_snapshot(results, 0, 0.0); write_snapshot(results, 1, 0.25)
        env = os.environ.copy(); env["MPLCONFIGDIR"] = str(root / "mpl")
        subprocess.run([sys.executable, str(script), "--results-dir", str(results), "--output-dir", str(output), "--max-slices", "2"], check=True, env=env)
        required = ["csf_slice_metrics.csv", "csf_slice_metrics.npz", "selected_snapshots.csv", "csf_slice_metrics_metadata.json", "flow_rate_mm3_per_s_heatmap.png", "flow_rate_mm3_per_s_heatmap.pdf"]
        assert all((output / name).is_file() for name in required)
        data = np.load(output / "csf_slice_metrics.npz"); assert data["area_mm2"].shape == (2,2); assert np.allclose(data["flow_rate_mm3_per_s"], 4 * data["area_mm2"]); assert np.allclose(data["mean_pressure"], 7.5); assert np.allclose(data["relative_area_deformation"], 0); assert np.allclose(data["max_tissue_displacement_mm"], 0.5)
        metadata = json.loads((output / "csf_slice_metrics_metadata.json").read_text()); assert metadata["geometry"]["reference_snapshot_warped"]
        with (output / "csf_slice_metrics.csv").open(newline="") as stream: assert len(list(csv.DictReader(stream))) == 4
        partials = root / "partials"; merged = root / "merged"
        for chunk in ("0/2", "1/2"):
            subprocess.run([sys.executable, str(script), "--results-dir", str(results), "--output-dir", str(partials), "--max-slices", "2", "--snapshot-chunk", chunk, "--no-plots"], check=True, env=env)
        subprocess.run([sys.executable, str(merger), str(partials), "--output-dir", str(merged), "--no-plots"], check=True, env=env)
        merged_data = np.load(merged / "csf_slice_metrics.npz")
        for key in ("times_s", "distance_from_FM_mm", "area_mm2", "flow_rate_mm3_per_s", "relative_area_deformation", "mean_pressure", "max_tissue_displacement_mm"): assert np.allclose(merged_data[key], data[key], equal_nan=True)
    print("CSF XDMF end-to-end test passed")


if __name__ == "__main__": main()
