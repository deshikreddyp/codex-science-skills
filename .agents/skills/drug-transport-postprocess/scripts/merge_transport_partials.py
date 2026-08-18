#!/usr/bin/env python3
"""Verify and merge collision-safe drug-transport snapshot chunks."""

from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np


def wide(path, times, distances, values):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream); writer.writerow(["time_seconds", *[f"{x:.17g}" for x in times]])
        for j, d in enumerate(distances): writer.writerow([f"{d:.17g}", *[f"{x:.17g}" for x in values[:, j]]])


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("partials_root", type=Path); p.add_argument("--output-dir", type=Path, required=True); p.add_argument("--overwrite", action="store_true"); p.add_argument("--no-plots", action="store_true"); a = p.parse_args()
    folders = sorted(x.parent for x in a.partials_root.glob("partial_*/concentration_slice_metrics_metadata.json"))
    if not folders: raise FileNotFoundError(f"No transport partials under {a.partials_root}")
    metas = [json.loads((f / "concentration_slice_metrics_metadata.json").read_text()) for f in folders]; first = metas[0]
    for folder, meta in zip(folders[1:], metas[1:]):
        for key in ("source_signature", "field_paths", "slice_convention"):
            if meta.get(key) != first.get(key): raise ValueError(f"Partial mismatch for {key}: {folder}")
    chunks = [np.load(folder / "concentration_slice_metrics.npz") for folder in folders]; d = chunks[0]["distance_from_FM_mm"]; z = chunks[0]["z_plane_mm"]
    for folder, chunk in zip(folders[1:], chunks[1:]):
        if not np.array_equal(chunk["distance_from_FM_mm"], d) or not np.allclose(chunk["z_plane_mm"], z): raise ValueError(f"Slice-coordinate mismatch: {folder}")
    times = np.concatenate([c["times_s"] for c in chunks]); order = np.argsort(times)
    if len(np.unique(times)) != len(times): raise ValueError("Duplicate physical times across partials")
    arrays = {key: np.concatenate([c[key] for c in chunks], axis=0)[order] for key in ("area_mm2", "c_net", "c_avg")}
    if any(np.isinf(values).any() for values in arrays.values()): raise ValueError("Infinite metric value in partials")
    times = times[order]
    if a.output_dir.exists() and not a.overwrite and (a.output_dir / "concentration_slice_metrics.npz").exists(): raise FileExistsError("Output exists; pass --overwrite")
    a.output_dir.mkdir(parents=True, exist_ok=True); np.savez_compressed(a.output_dir / "concentration_slice_metrics.npz", times_s=times, distance_from_FM_mm=d, z_plane_mm=z, **arrays); wide(a.output_dir / "net_concentration.csv", times, d, arrays["c_net"]); wide(a.output_dir / "avg_concentration.csv", times, d, arrays["c_avg"])
    selected = []
    for folder in folders:
        with (folder / "selected_timesteps.csv").open(newline="", encoding="utf-8") as stream: selected.extend(csv.DictReader(stream))
    selected.sort(key=lambda row: float(row["time_s"]))
    with (a.output_dir / "selected_timesteps.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["snapshot_index", "time_s", "source_file"]); writer.writeheader(); writer.writerows(selected)
    with (a.output_dir / "concentration_slice_metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        fields = ["snapshot_index", "time_s", "distance_from_FM_mm", "z_plane_mm", "area_mm2", "c_net", "c_avg"]; writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for i, snapshot in enumerate(selected):
            for j, distance in enumerate(d): writer.writerow({"snapshot_index": snapshot["snapshot_index"], "time_s": times[i], "distance_from_FM_mm": distance, "z_plane_mm": z[j], **{key: values[i, j] for key, values in arrays.items()}})
    first["snapshot_chunk"] = None; first["merged_partials"] = [str(f) for f in folders]; first["merge_verification"] = {"source_signature": True, "field_paths": True, "slice_coordinates": True, "unique_times": True, "no_infinite_values": True, "nan_slices_allowed": True}
    (a.output_dir / "concentration_slice_metrics_metadata.json").write_text(json.dumps(first, indent=2) + "\n")
    if not a.no_plots:
        from plot_concentration_profiles import plot_all
        plot_all(a.output_dir)
    print(f"Merged {len(folders)} partials into {a.output_dir}")


if __name__ == "__main__": main()
