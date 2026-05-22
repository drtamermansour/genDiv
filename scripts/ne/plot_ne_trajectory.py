#!/usr/bin/env python3
"""Plot per-group Ne trajectory from an Ne_<tool>_summary.csv.

Default layout (`--broken-axis`, recommended for Standardbred GONE2 / SNeP
outputs): a two-panel broken-x figure with linear x. Left panel covers
recent generations (gen 1 - split-gen) with breed-history annotations;
right panel covers the ancestral tail to gen max-gen. The y-axis is log
because Ne spans two orders of magnitude across breed founding.

A single-panel log-x view (`--log-x`) is preserved for backward
compatibility.

Annotations:
  * Shaded band for gen 1-4 — GONE block-based estimator returns these as
    a single identical value (Novo et al. 2023).
  * Standardbred-specific historical reference lines:
      - 2011 (cap standardised to 140 mares/stallion both gaits)  -> gen 1.36
      - 2009 (USTA studbook cap introduced)                       -> gen 1.55
      - 1973 (closure of the Standardbred studbook)               -> gen 4.82
      - 1872 (breed founding, Hambletonian 10 era)                -> gen 14.0
    Disable with `--no-standardbred-events`.
  * Secondary x-axis (top) showing approximate calendar year at G = 11.0
    yr/gen (Waples 2026).

Usage:
    python plot_ne_trajectory.py \
        --summary-csv results_*/divStats/ne/gone2/Ne_gone2_summary.csv \
        --out results_*/divStats/ne/gone2/Ne_gone2_trajectory.png \
        --max-gen 50
"""
from __future__ import annotations

import argparse
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

GROUP_ORDER = ["wholePop", "Trotter", "Pacer"]
GROUP_COLOR = {
    "wholePop": "#2C3E50",
    "Trotter":  "#1F77B4",
    "Pacer":    "#D62728",
}
GROUP_MARKER = {
    "wholePop": "o",
    "Trotter":  "s",
    "Pacer":    "^",
}
# When multiple Methods are present, cycle through these linestyles so each
# (Group, Method) pair becomes a distinct line in the plot.
METHOD_LINESTYLES = ["-", "--", ":", "-."]

# Standardbred-specific historical reference events. Each entry is
# (year, label, line-style-hint, vertical-anchor 0-1 for label y).
# The 2009 USTA cap and the 2011 standardisation are 0.2 generations apart
# and both sit inside the GONE artifact zone (gen 1.36 and gen 1.55 at
# G=11.0); we collapse them into a single label drawn at the 2009 line.
SBRED_EVENTS = [
    (2011, None,                                            "in_artifact",            0.55),
    (2009, "USTA caps\n(2009 / 2011)",                      "in_artifact",            0.55),
    (1973, "Studbook\nclosure\n(1973)",                     "interpretable_boundary", 0.80),
    (1872, "Breed founding\n~1872",                         "breed_founding",         0.30),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary-csv", required=True,
                   help="Ne_<tool>_summary.csv with columns Group,Method,Generations_ago,Ne,...")
    p.add_argument("--out", required=True,
                   help="Output PNG path (a .pdf will be written alongside).")
    p.add_argument("--method", default=None,
                   help="Filter Method column to this value (default: use all rows).")
    p.add_argument("--max-gen", type=int, default=50,
                   help="Truncate the x-axis to this many generations (default 50).")
    p.add_argument("--min-gen", type=int, default=1,
                   help="Lower x-axis bound, generations (default 1).")
    p.add_argument("--split-gen", type=float, default=15.0,
                   help="Generation at which to split the broken x-axis "
                        "(default 15; below = recent panel, above = ancestral panel).")
    p.add_argument("--artifact-max-gen", type=int, default=4,
                   help="Highest generation in the GONE block-based artifact zone "
                        "to shade (default 4; Novo et al. 2023).")
    p.add_argument("--gen-time", type=float, default=11.0,
                   help="Years per generation for the top calendar axis "
                        "(default 11.0; Waples 2026 Standardbred mean).")
    p.add_argument("--current-year", type=int, default=2026,
                   help="Sampling year for converting generations -> calendar year.")
    p.add_argument("--no-standardbred-events", action="store_true",
                   help="Disable the four Standardbred-specific reference lines.")
    p.add_argument("--log-x", action="store_true",
                   help="Use single-panel log-x layout (legacy) instead of the "
                        "default broken-linear-x two-panel layout.")
    p.add_argument("--linear-y", action="store_true",
                   help="Use linear y-axis instead of log (default log).")
    p.add_argument("--title", default=None,
                   help="Plot title (default derived from method + group count).")
    return p.parse_args()


def _ordered_groups(present: list[str]) -> list[str]:
    out = [g for g in GROUP_ORDER if g in present]
    for g in present:
        if g not in out:
            out.append(g)
    return out


def _plot_lines(ax, df, groups, methods=None):
    """Plot one line per (group, method) pair. Color encodes group,
    linestyle encodes method. When only one method is present the
    linestyle defaults to solid; with multiple methods we cycle through
    METHOD_LINESTYLES so each method gets a distinct style."""
    if methods is None or len(methods) <= 1:
        method_linestyle = {(methods[0] if methods else None): METHOD_LINESTYLES[0]}
    else:
        method_linestyle = {m: METHOD_LINESTYLES[i % len(METHOD_LINESTYLES)]
                            for i, m in enumerate(methods)}
    for g in groups:
        for m, ls in method_linestyle.items():
            if m is None:
                sub = df[df["Group"] == g].sort_values("Generations_ago")
            else:
                sub = df[(df["Group"] == g) & (df["Method"] == m)].sort_values("Generations_ago")
            if sub.empty:
                continue
            color = GROUP_COLOR.get(g, "grey")
            marker = GROUP_MARKER.get(g, "o")
            label = g if (m is None or len(method_linestyle) == 1) else f"{g} ({m})"
            ax.plot(sub["Generations_ago"], sub["Ne"],
                    color=color, linestyle=ls, marker=marker,
                    linewidth=1.6, markersize=4.0, alpha=0.95,
                    label=label)


def _draw_artifact_band(ax, min_gen, artifact_max_gen, label_y_pos=None,
                        annotate=True):
    """Shade the GONE block-based-estimator artifact zone (gen 1-4 by default)."""
    if artifact_max_gen < min_gen:
        return
    ax.axvspan(min_gen - 0.5, artifact_max_gen + 0.5,
               color="lightgrey", alpha=0.40, zorder=0)
    if not annotate:
        return
    y0, y1 = ax.get_ylim()
    if ax.get_yscale() == "log":
        y_label = (y1 ** 0.88) * (y0 ** 0.12)
    else:
        y_label = y0 + 0.88 * (y1 - y0)
    mid = (min_gen + artifact_max_gen) / 2.0
    ax.text(mid, y_label,
            "GONE artifact\nzone (Novo 2023)",
            ha="center", va="center",
            fontsize=7.5, color="dimgrey", fontstyle="italic")


def _draw_events(ax, events, gen_time, current_year, panel_gen_range):
    """Draw historical reference lines that fall inside `panel_gen_range`.

    Each event tuple is (year, label_or_None, kind, y_anchor) where
    y_anchor in [0,1] is the relative vertical position of the label
    inside the panel (used to stagger labels so they don't overlap).
    A None label still draws the line but no text.
    """
    lo, hi = panel_gen_range
    y0, y1 = ax.get_ylim()
    for year, label, kind, y_anchor in events:
        gen = (current_year - year) / gen_time
        if not (lo <= gen <= hi):
            continue
        # Visual style depends on which "zone" the event sits in
        if kind == "in_artifact":
            style = dict(linestyle=":", color="#888", linewidth=0.9, alpha=0.7)
            label_color = "#555"
            label_fontsize = 7.5
        elif kind == "interpretable_boundary":
            style = dict(linestyle="--", color="#333", linewidth=1.2, alpha=0.9)
            label_color = "#222"
            label_fontsize = 8.0
        elif kind == "breed_founding":
            style = dict(linestyle="--", color="grey", linewidth=1.2, alpha=0.85)
            label_color = "dimgrey"
            label_fontsize = 8.0
        else:
            style = dict(linestyle="--", color="grey", linewidth=1.0, alpha=0.6)
            label_color = "grey"
            label_fontsize = 8.0
        ax.axvline(gen, **style)
        if label is None:
            continue
        if ax.get_yscale() == "log":
            y_label = (y1 ** y_anchor) * (y0 ** (1 - y_anchor))
        else:
            y_label = y0 + y_anchor * (y1 - y0)
        # Place label slightly to the right of the event line, unless the
        # event is in the right half of the panel where right-of-line would
        # overflow — in that case place it to the left.
        x_frac = (gen - lo) / max(1e-9, hi - lo)
        if x_frac < 0.7:
            ha, x_off = "left", gen + (hi - lo) * 0.02
        else:
            ha, x_off = "right", gen - (hi - lo) * 0.02
        ax.text(x_off, y_label, label,
                fontsize=label_fontsize, color=label_color,
                ha=ha, va="center",
                bbox=dict(boxstyle="round,pad=0.20",
                          facecolor="white",
                          edgecolor="lightgrey",
                          alpha=0.9))


def _add_year_axis(ax, gen_time, current_year, gen_ticks, label=True):
    """Attach a secondary x-axis with calendar-year tick labels at the
    given generation positions."""
    secax = ax.twiny()
    secax.set_xscale(ax.get_xscale())
    secax.set_xlim(ax.get_xlim())
    secax.set_xticks(gen_ticks)
    secax.set_xticklabels([f"{int(round(current_year - g * gen_time))}" for g in gen_ticks],
                          fontsize=8)
    secax.minorticks_off()
    if label:
        secax.set_xlabel(f"Approximate calendar year (G = {gen_time:.1f} yr/gen, Waples 2026)",
                         fontsize=9)


def _broken_axis_marks(ax_left, ax_right, d=0.018):
    """Draw the diagonal break marks on the inner edges of two adjacent
    axes that share an x-coordinate range with a discontinuity."""
    kwargs = dict(transform=ax_left.transAxes, color="k", clip_on=False, lw=1)
    ax_left.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    ax_left.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
    kwargs.update(transform=ax_right.transAxes)
    ax_right.plot((-d, +d), (-d, +d), **kwargs)
    ax_right.plot((-d, +d), (1 - d, 1 + d), **kwargs)


def plot_log_x(args, df, methods, groups):
    """Single-panel log-x layout (legacy)."""
    G, yr = args.gen_time, args.current_year
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    _plot_lines(ax, df, groups, methods=list(methods))
    ax.set_xscale("log")
    if not args.linear_y:
        ax.set_yscale("log")
    ax.set_xlim(args.min_gen, args.max_gen)
    ax.set_xlabel("Generations ago (log scale)")
    ax.set_ylabel("Effective population size, $N_e$"
                  + ("" if args.linear_y else " (log scale)"))
    _draw_artifact_band(ax, args.min_gen, args.artifact_max_gen)
    if not args.no_standardbred_events:
        _draw_events(ax, SBRED_EVENTS, G, yr,
                     (args.min_gen, args.max_gen))
    tick_gens = [g for g in (1, 5, 10, 14, 20, 30, 50)
                 if args.min_gen <= g <= args.max_gen]
    _add_year_axis(ax, G, yr, tick_gens)
    return fig, (ax,)


def plot_broken(args, df, methods, groups):
    """Two-panel broken-linear-x layout (default)."""
    G, yr = args.gen_time, args.current_year
    split = args.split_gen
    width_ratio_left = max(1.5, split - args.min_gen)
    width_ratio_right = max(0.6, args.max_gen - split)
    fig, (ax_l, ax_r) = plt.subplots(
        1, 2, figsize=(10.5, 5.5),
        gridspec_kw={"width_ratios": [width_ratio_left, width_ratio_right],
                     "wspace": 0.05},
    )
    for ax in (ax_l, ax_r):
        _plot_lines(ax, df, groups, methods=list(methods))
        if not args.linear_y:
            ax.set_yscale("log")
        ax.grid(True, which="both", linestyle=":", alpha=0.30)

    ax_l.set_xlim(args.min_gen - 0.5, split + 0.5)
    ax_r.set_xlim(split - 0.5, args.max_gen + 0.5)

    # Style the inner spines / ticks for the broken-axis effect
    ax_l.spines["right"].set_visible(False)
    ax_r.spines["left"].set_visible(False)
    ax_r.tick_params(left=False, labelleft=False)
    _broken_axis_marks(ax_l, ax_r)

    # x-tick locations: data ticks shown every odd generation on the
    # recent panel; year-axis ticks use a sparser subset to avoid label
    # collisions.
    left_data_ticks  = [g for g in range(int(args.min_gen),
                                         int(split) + 1) if g % 2 == 1 or g in (args.min_gen,)]
    right_data_ticks = [g for g in range(int(split), int(args.max_gen) + 1, 5)
                        if g >= split]
    ax_l.set_xticks(left_data_ticks)
    ax_r.set_xticks(right_data_ticks)

    # Year-axis ticks: pick ~4-6 evenly-spaced positions per panel so the
    # calendar labels don't collide. Works for any data range (e.g. GONE2
    # starts at gen 1, SNeP starts at gen ~12).
    def _spaced_year_ticks(lo, hi, target=5):
        span = hi - lo
        if span <= 0:
            return []
        # Round to a "nice" step that gives ~target ticks
        for step in (1, 2, 3, 5, 10, 20, 50, 100):
            if span / step <= target + 1:
                break
        first = (int(lo) // step) * step
        if first < lo:
            first += step
        return list(range(int(first), int(hi) + 1, step))

    left_year_ticks  = _spaced_year_ticks(args.min_gen, split, target=4)
    right_year_ticks = _spaced_year_ticks(split, args.max_gen, target=5)
    if not left_year_ticks:
        left_year_ticks = left_data_ticks
    if not right_year_ticks:
        right_year_ticks = right_data_ticks

    # Annotations on the recent panel only
    _draw_artifact_band(ax_l, args.min_gen, args.artifact_max_gen)
    if not args.no_standardbred_events:
        _draw_events(ax_l, SBRED_EVENTS, G, yr,
                     (args.min_gen - 0.5, split + 0.5))
        _draw_events(ax_r, SBRED_EVENTS, G, yr,
                     (split - 0.5, args.max_gen + 0.5))

    # Calendar-year secondary axes per panel (no per-axis labels — one
    # shared figure-level label sits above them).
    _add_year_axis(ax_l, G, yr, left_year_ticks, label=False)
    _add_year_axis(ax_r, G, yr, right_year_ticks, label=False)

    ax_l.set_ylabel("Effective population size, $N_e$"
                    + ("" if args.linear_y else " (log scale)"))
    # Hide x-axis labels at the bottom of each subplot since we use a shared figure label
    ax_l.set_xlabel("")
    ax_r.set_xlabel("")

    return fig, (ax_l, ax_r)


def main() -> int:
    args = parse_args()

    df = pd.read_csv(args.summary_csv)
    required = {"Group", "Method", "Generations_ago", "Ne"}
    if not required.issubset(df.columns):
        sys.exit(f"[plot_ne_trajectory] CSV missing required columns: "
                 f"{required - set(df.columns)}")

    if args.method:
        df = df[df["Method"] == args.method]
    df = df[(df["Generations_ago"] >= args.min_gen) &
            (df["Generations_ago"] <= args.max_gen)].copy()
    df["Ne"] = pd.to_numeric(df["Ne"], errors="coerce")
    df = df.dropna(subset=["Ne"])

    if df.empty:
        sys.exit(f"[plot_ne_trajectory] no rows after filtering "
                 f"(method={args.method}, gen in [{args.min_gen},{args.max_gen}])")

    methods = df["Method"].unique()
    groups = _ordered_groups(list(df["Group"].unique()))

    if args.log_x:
        fig, axes = plot_log_x(args, df, methods, groups)
    else:
        fig, axes = plot_broken(args, df, methods, groups)

    # Title
    if args.title is None:
        method_label = methods[0] if len(methods) == 1 else "+".join(methods)
        args.title = (f"{method_label} per-generation $N_e$ trajectory "
                      f"(gen {args.min_gen}-{args.max_gen})")

    if not args.log_x:
        # Two-panel: leave clear separation between title, year axis,
        # year-axis label, and the plot area.
        G = args.gen_time
        fig.subplots_adjust(top=0.82, bottom=0.13, left=0.085, right=0.97)
        # Title at the very top
        fig.text(0.5, 0.96, args.title, ha="center", va="center",
                 fontsize=11, fontweight="bold")
        # Year-axis caption sits between the title and the year tick labels
        fig.text(0.5, 0.91,
                 f"Approximate calendar year (G = {G:.1f} yr/gen, Waples 2026)",
                 ha="center", va="center", fontsize=9, color="dimgrey")
        # Shared bottom axis label
        fig.text(0.5, 0.04, "Generations ago", ha="center", va="center", fontsize=10)
        # Legend on the right panel where there are no annotation labels
        axes[-1].legend(loc="lower right", frameon=False, fontsize=9)
    else:
        fig.suptitle(args.title, y=0.995, fontsize=11)
        axes[0].legend(loc="upper left", frameon=False, fontsize=9)
        fig.tight_layout()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    pdf = os.path.splitext(args.out)[0] + ".pdf"
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot_ne_trajectory] wrote {args.out}")
    print(f"[plot_ne_trajectory] wrote {pdf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
