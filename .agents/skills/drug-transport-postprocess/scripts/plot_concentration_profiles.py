#!/usr/bin/env python3
"""Generate Paper2-style c_avg/c_net heatmaps and selected-time profiles."""

from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def style():
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Nimbus Roman", "Times New Roman", "DejaVu Serif"], "font.size": 10})


def parse_numbers(spec): return [] if not spec else [float(x) for x in spec.split(",") if x.strip()]


def plot_all(output_dir: Path, profile_times=None, cavg_limits=None, cnet_limits=None):
    style(); data = np.load(output_dir / "concentration_slice_metrics.npz"); t = data["times_s"]; d = data["distance_from_FM_mm"]
    for key, label, limits in (("c_avg", "Average concentration", cavg_limits), ("c_net", r"Cross-sectional $c_{net}$", cnet_limits)):
        values = np.asarray(data[key])
        if not np.isfinite(values).any(): continue
        fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True); kwargs = {} if limits is None else {"vmin": limits[0], "vmax": limits[1]}
        mesh = ax.pcolormesh(t, d, values.T, shading="auto", cmap="magma", **kwargs); ax.set(xlabel="Time (s)", ylabel="Distance from FM (mm)"); ax.invert_yaxis(); fig.colorbar(mesh, ax=ax, label=label)
        fig.savefig(output_dir / f"{key}_heatmap.png", dpi=300); fig.savefig(output_dir / f"{key}_heatmap.pdf"); plt.close(fig)
        requested = parse_numbers(profile_times)
        if requested:
            fig, ax = plt.subplots(figsize=(5, 5), constrained_layout=True)
            for value in requested:
                i = int(np.nanargmin(np.abs(t - value))); ax.plot(values[i], d, label=f"{t[i]:g} s")
            ax.set(xlabel=label, ylabel="Distance from FM (mm)"); ax.invert_yaxis(); ax.legend(frameon=False)
            fig.savefig(output_dir / f"{key}_profiles.png", dpi=300); fig.savefig(output_dir / f"{key}_profiles.pdf"); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("output_dir", type=Path); p.add_argument("--profile-times-s"); p.add_argument("--cavg-limits", nargs=2, type=float); p.add_argument("--cnet-limits", nargs=2, type=float); a = p.parse_args(); plot_all(a.output_dir, a.profile_times_s, a.cavg_limits, a.cnet_limits)


if __name__ == "__main__": main()
