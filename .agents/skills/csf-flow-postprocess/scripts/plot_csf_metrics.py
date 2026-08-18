#!/usr/bin/env python3
"""Plot CSF metric heatmaps, optional waveforms, and optional spatial profiles."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PLOTS = {
    "flow_rate_mm3_per_s": ("Flow rate (mL/s)", 1e-3, "RdBu_r", True),
    "area_mm2": (r"Area (mm$^2$)", 1.0, "viridis", False),
    "relative_area_deformation": ("Relative area deformation (%)", 100.0, "RdBu_r", True),
    "mean_pressure": ("Area-weighted mean pressure", 1.0, "RdBu_r", True),
    "max_tissue_displacement_mm": ("Maximum tissue displacement (mm)", 1.0, "magma", False),
}


def style() -> None:
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Nimbus Roman", "Times New Roman", "DejaVu Serif"], "font.size": 10})


def parse_numbers(spec: str | None) -> list[float]:
    return [] if not spec else [float(v) for v in spec.split(",") if v.strip()]


def nearest(values: np.ndarray, requested: float) -> int:
    return int(np.nanargmin(np.abs(values - requested)))


def plot_all(output_dir: Path, waveform_distances: str | None = None, profile_times: str | None = None) -> None:
    style()
    data = np.load(output_dir / "csf_slice_metrics.npz")
    t, d = data["times_s"], data["distance_from_FM_mm"]
    for key, (label, scale, cmap, symmetric) in PLOTS.items():
        values = np.asarray(data[key]) * scale
        if not np.isfinite(values).any():
            continue
        fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
        kwargs = {}
        if symmetric:
            limit = float(np.nanmax(np.abs(values))); kwargs = {"vmin": -limit, "vmax": limit}
        mesh = ax.pcolormesh(t, d, values.T, shading="auto", cmap=cmap, **kwargs)
        ax.set_xlabel("Time (s)"); ax.set_ylabel("Distance from FM (mm)"); ax.invert_yaxis()
        fig.colorbar(mesh, ax=ax, label=label)
        for suffix in ("png", "pdf"):
            fig.savefig(output_dir / f"{key}_heatmap.{suffix}", dpi=300)
        plt.close(fig)
    requested_d = parse_numbers(waveform_distances)
    requested_t = parse_numbers(profile_times)
    for key, (label, scale, _, _) in PLOTS.items():
        values = np.asarray(data[key]) * scale
        if requested_d and np.isfinite(values).any():
            fig, ax = plt.subplots(figsize=(6.5, 4), constrained_layout=True)
            for value in requested_d:
                j = nearest(d, value); ax.plot(t, values[:, j], label=f"{d[j]:g} mm")
            ax.set(xlabel="Time (s)", ylabel=label); ax.legend(frameon=False)
            fig.savefig(output_dir / f"{key}_waveforms.png", dpi=300); fig.savefig(output_dir / f"{key}_waveforms.pdf"); plt.close(fig)
        if requested_t and np.isfinite(values).any():
            fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
            for value in requested_t:
                i = nearest(t, value); ax.plot(values[i], d, label=f"{t[i]:g} s")
            ax.set(xlabel=label, ylabel="Distance from FM (mm)"); ax.invert_yaxis(); ax.legend(frameon=False)
            fig.savefig(output_dir / f"{key}_profiles.png", dpi=300); fig.savefig(output_dir / f"{key}_profiles.pdf"); plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output_dir", type=Path)
    p.add_argument("--waveform-distances-mm")
    p.add_argument("--profile-times-s")
    a = p.parse_args(); plot_all(a.output_dir, a.waveform_distances_mm, a.profile_times_s)


if __name__ == "__main__":
    main()
