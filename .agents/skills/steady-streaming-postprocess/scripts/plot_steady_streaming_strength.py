#!/usr/bin/env python3
"""Create publication-ready steady-streaming strength or mean-velocity profiles."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


PROFILE_COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
DIRECTION_COLORS = {"cranial": "#B2182B", "caudal": "#2166AC", "net": "#222222"}
LINESTYLES = ("-", (0, (5.0, 2.4)), (0, (1.3, 1.8)), (0, (7.0, 2.0, 1.4, 2.0)))


def parse_mapping(items: list[str], option: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"{option} value must have form LABEL=VALUE: {item!r}")
        label, value = item.split("=", maxsplit=1)
        if not label or not value:
            raise ValueError(f"{option} value must have form LABEL=VALUE: {item!r}")
        result[label] = value
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot standardized CSV profiles from compute_strength_by_slices.py as PNG, PDF, and SVG. "
            "No title is added."
        )
    )
    parser.add_argument(
        "--profile",
        action="append",
        required=True,
        help="Repeated LABEL=/path/to/profile.csv specification.",
    )
    parser.add_argument(
        "--case-filter",
        action="append",
        default=[],
        help="Optional repeated LABEL=case filter for multi-case CSV files.",
    )
    parser.add_argument("--profile-color", action="append", default=[], help="Optional LABEL=#RRGGBB.")
    parser.add_argument(
        "--profile-style",
        action="append",
        default=[],
        help="Optional LABEL=solid|dashed|dotted|dashdot.",
    )
    parser.add_argument("--mode", choices=("signed", "cranial", "caudal", "net"), default="signed")
    parser.add_argument("--quantity", choices=("strength", "mean-velocity"), default="strength")
    parser.add_argument("--color-by", choices=("direction", "profile"), default="direction")
    parser.add_argument("--absolute-caudal", action="store_true")
    parser.add_argument("--distance-unit", choices=("mm", "cm"), default="cm")
    parser.add_argument("--output-stem", type=Path, required=True)
    parser.add_argument("--font-family", default="Nimbus Roman")
    parser.add_argument("--font-size", type=float, default=22.0)
    parser.add_argument("--legend-font-size", type=float, default=18.0)
    parser.add_argument("--width-in", type=float, default=9.8)
    parser.add_argument("--height-in", type=float, default=5.35)
    parser.add_argument("--line-width", type=float, default=3.0)
    parser.add_argument("--dpi", type=int, default=600)
    return parser.parse_args()


def style_value(name: str):
    styles = {
        "solid": "-",
        "dashed": (0, (5.0, 2.4)),
        "dotted": (0, (1.3, 1.8)),
        "dashdot": (0, (7.0, 2.0, 1.4, 2.0)),
    }
    if name not in styles:
        raise ValueError(f"Unknown line style {name!r}; use {', '.join(styles)}")
    return styles[name]


def read_profile(path: Path, case_filter: str | None) -> dict[str, np.ndarray]:
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if case_filter is None or row.get("case") == case_filter:
                rows.append(row)
    if not rows:
        suffix = f" with case={case_filter}" if case_filter else ""
        raise ValueError(f"No rows found in {path}{suffix}")
    result: dict[str, np.ndarray] = {}
    for column in (
        "distance_from_fm_mm",
        "cranial_strength_mL_s",
        "caudal_strength_mL_s",
        "net_strength_mL_s",
        "cranial_mean_velocity_mm_s",
        "caudal_mean_velocity_mm_s",
        "net_mean_velocity_mm_s",
    ):
        if column not in rows[0]:
            raise KeyError(f"Required column {column!r} is missing from {path}")
        result[column] = np.asarray([float(row[column]) for row in rows], dtype=np.float64)
    order = np.argsort(result["distance_from_fm_mm"])
    return {name: values[order] for name, values in result.items()}


def main() -> None:
    args = parse_args()
    profiles = parse_mapping(args.profile, "--profile")
    case_filters = parse_mapping(args.case_filter, "--case-filter")
    custom_colors = parse_mapping(args.profile_color, "--profile-color")
    custom_styles = parse_mapping(args.profile_style, "--profile-style")
    unknown = set(case_filters) | set(custom_colors) | set(custom_styles)
    unknown -= set(profiles)
    if unknown:
        raise ValueError(f"Options refer to unknown profile labels: {sorted(unknown)}")

    plt.rcParams.update(
        {
            "font.family": args.font_family,
            "mathtext.fontset": "stix",
            "axes.labelsize": args.font_size,
            "xtick.labelsize": args.font_size,
            "ytick.labelsize": args.font_size,
            "legend.fontsize": args.legend_font_size,
            "axes.linewidth": 1.2,
            "xtick.major.width": 1.1,
            "ytick.major.width": 1.1,
            "xtick.major.size": 5.5,
            "ytick.major.size": 5.5,
            "savefig.dpi": args.dpi,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(args.width_in, args.height_in))
    fig.subplots_adjust(left=0.15, right=0.985, bottom=0.20, top=0.86)

    profile_handles: list[Line2D] = []
    directions = ("cranial", "caudal") if args.mode == "signed" else (args.mode,)
    for profile_index, (label, path_text) in enumerate(profiles.items()):
        data = read_profile(Path(path_text), case_filters.get(label))
        distance = data["distance_from_fm_mm"] / (10.0 if args.distance_unit == "cm" else 1.0)
        profile_color = custom_colors.get(label, PROFILE_COLORS[profile_index % len(PROFILE_COLORS)])
        profile_style = (
            style_value(custom_styles[label])
            if label in custom_styles
            else LINESTYLES[profile_index % len(LINESTYLES)]
        )
        for direction in directions:
            if args.quantity == "strength":
                column = f"{direction}_strength_mL_s"
            else:
                column = f"{direction}_mean_velocity_mm_s"
            values = data[column].copy()
            if direction == "caudal" and args.absolute_caudal:
                values = np.abs(values)
            if args.color_by == "direction":
                color = DIRECTION_COLORS[direction]
                linestyle = profile_style
            else:
                color = profile_color
                linestyle = "-" if direction in ("cranial", "net") else (0, (5.0, 2.4))
            plot_label = label if len(directions) == 1 else None
            ax.plot(
                distance,
                values,
                color=color,
                linestyle=linestyle,
                linewidth=args.line_width,
                label=plot_label,
            )
        if len(directions) > 1:
            handle_color = "0.15" if args.color_by == "direction" else profile_color
            handle_style = profile_style if args.color_by == "direction" else "-"
            profile_handles.append(
                Line2D([0], [0], color=handle_color, linestyle=handle_style, linewidth=args.line_width, label=label)
            )

    ax.axhline(0.0, color="0.22", linewidth=1.0)
    ax.set_xlabel(f"Distance from FM ({args.distance_unit})", labelpad=8)
    if args.quantity == "strength":
        ax.set_ylabel(r"$Q_{\mathrm{SS}}$ (mL s$^{-1}$)", labelpad=12)
    else:
        ax.set_ylabel(r"Mean steady-streaming velocity (mm s$^{-1}$)", labelpad=12)
    ax.grid(True, color="0.88", linewidth=1.0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.margins(x=0.01)

    if len(directions) == 1:
        ax.legend(frameon=False, loc="best", handlelength=2.6)
    else:
        first_legend = ax.legend(
            handles=profile_handles,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, 1.20),
            ncol=min(4, len(profile_handles)),
            handlelength=2.6,
        )
        ax.add_artist(first_legend)
        if args.color_by == "direction":
            direction_handles = [
                Line2D([0], [0], color=DIRECTION_COLORS["cranial"], linewidth=args.line_width, label="Cranial"),
                Line2D([0], [0], color=DIRECTION_COLORS["caudal"], linewidth=args.line_width, label="Caudal"),
            ]
        else:
            direction_handles = [
                Line2D([0], [0], color="0.15", linewidth=args.line_width, linestyle="-", label="Cranial"),
                Line2D(
                    [0],
                    [0],
                    color="0.15",
                    linewidth=args.line_width,
                    linestyle=(0, (5.0, 2.4)),
                    label="Caudal",
                ),
            ]
        ax.legend(handles=direction_handles, frameon=False, loc="best", handlelength=2.6)

    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    for extension in ("png", "pdf", "svg"):
        output = args.output_stem.with_suffix(f".{extension}")
        fig.savefig(output, bbox_inches="tight", pad_inches=0.08)
        print(output, flush=True)
    plt.close(fig)


if __name__ == "__main__":
    main()
