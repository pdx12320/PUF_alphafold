#!/usr/bin/env python3
"""Regenerate the Wiki figures from recorded repository data (no retraining).

Run from any directory:
    python wiki/scripts/make_figures.py

Optional extended layout QA:
    python wiki/scripts/make_figures.py --qa-tools-dir /path/to/figure-qa/scripts

The optional directory provides audit_panel_alignment.py. The standalone run
still measures all plot areas and blocks unequal comparable panel geometry.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "wiki/figures"
DATA = ROOT / "wiki/figure_data"
QA = OUT / "qa"
MM = 1 / 25.4
BLUE = "#287C9D"
ORANGE = "#C47742"
DARK = "#243440"
GREY = "#D9E1E5"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "sans-serif"],
    "font.size": 7,
    "axes.titlesize": 8.5,
    "axes.labelsize": 7,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 7,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "text.color": DARK,
    "axes.labelcolor": DARK,
    "xtick.color": DARK,
    "ytick.color": DARK,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.facecolor": "white",
})

SOURCES = {
    "scaffold_metrics": "research/results_20260922/scaffold_RF/results/metrics.csv",
    "scaffold_predictions": "research/results_20260922/scaffold_RF/results/predictions_indexed.csv",
    "scaffold_labels": "research/results_20260922/scaffold_RF/results/construct_ids.csv",
    "five_construct_predictions": "research/results_20260922/five_construct_holdout/predictions.csv",
    "mpnn_site_scores": "research/aice_mpnn_20260922/results/aice_single_ranked.csv",
}
SOURCE_FRAMES = {}
require_matplotlib_panel_alignment = None


def read_source(key):
    frame = pd.read_csv(ROOT / SOURCES[key])
    SOURCE_FRAMES[key] = frame
    return frame


def panel_label(ax, label):
    ax.annotate(label, (0, 1), xycoords="axes fraction", xytext=(-21, 11),
                textcoords="offset points", weight="bold", fontsize=8.5,
                ha="left", va="bottom", annotation_clip=False)


def alignment_gate(fig, stem, axes):
    """Measure final physical plot areas; use the extended auditor if supplied."""
    fig.canvas.draw()
    if require_matplotlib_panel_alignment is not None:
        require_matplotlib_panel_alignment(fig, axes=axes,
                       json_out=QA / f"{stem}.alignment.json",
                       overlay_svg=QA / f"{stem}.alignment.svg",
                       tolerance_pt=1.5, gutter_tolerance_pt=1.5,
                       require_panel_labels=len(axes) > 1, strict=True)
        return
    size = fig.get_size_inches() * 72
    rectangles = [np.array([ax.get_position().x0, ax.get_position().y0,
                            ax.get_position().x1, ax.get_position().y1])
                  * np.tile(size, 2) for ax in axes]
    checks = {}
    if len(rectangles) > 1:
        rect = np.array(rectangles)
        checks = {"bottom_range_pt": float(np.ptp(rect[:, 1])),
                  "top_range_pt": float(np.ptp(rect[:, 3])),
                  "width_range_pt": float(np.ptp(rect[:, 2] - rect[:, 0])),
                  "height_range_pt": float(np.ptp(rect[:, 3] - rect[:, 1]))}
        if len(rectangles) > 2:
            checks["gutter_range_pt"] = float(np.ptp(rect[1:, 0] - rect[:-1, 2]))
        assert all(v <= 1.5 for v in checks.values()), checks
    report = {"backend": "python-matplotlib", "tolerance_pt": 1.5,
              "verdict": "PASS" if checks else "NOT APPLICABLE",
              "plot_area_rectangles_pt": [x.tolist() for x in rectangles],
              "checks": checks}
    (QA / f"{stem}.alignment.json").write_text(json.dumps(report, indent=2) + "\n")


def save_figure(fig, stem, axes):
    alignment_gate(fig, stem, axes)
    # Fixed dimensions preserve the measured physical layout; text is editable.
    active = {"fig5_model3_nontrm_consensus", "fixed4_scaffold_test", "training20_s12", "fixed5_endpoint_test"}
    destination = OUT if stem in active else ROOT / "research/archive/wiki_figures"
    destination.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination / f"{stem}.svg")
    fig.savefig(destination / f"{stem}.pdf")
    fig.savefig(destination / f"{stem}.png", dpi=600)
    plt.close(fig)


def matrix_counts(actual, predicted, labels):
    lookup = {v: i for i, v in enumerate(labels)}
    counts = np.zeros((len(labels), len(labels)), dtype=int)
    for a, p in zip(actual, predicted):
        counts[lookup[a], lookup[p]] += 1
    assert counts.sum() == len(actual)
    return counts


def confusion_panel(ax, counts, display, title, vmax):
    cmap = plt.get_cmap("Blues")
    norm = Normalize(0, vmax)
    ax.imshow(counts, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    for row in range(len(display)):
        for col in range(len(display)):
            val = counts[row, col]
            ax.text(col, row, str(val), ha="center", va="center", fontsize=11,
                    color="white" if norm(val) > .55 else DARK,
                    weight="bold" if row == col else "normal")
    ax.set_xticks(range(len(display)), display)
    ax.set_yticks(range(len(display)), display)
    ax.set_xlabel("Predicted class", labelpad=7)
    ax.set_ylabel("Measured class", labelpad=7)
    ax.set_title(title, pad=12)
    ax.tick_params(length=0, pad=5)
    for spine in ax.spines.values():
        spine.set_visible(False)


def model1_figures():
    metrics = read_source("scaffold_metrics")
    predictions = read_source("scaffold_predictions")
    labels = read_source("scaffold_labels")
    joined = predictions.merge(labels, on="id", validate="many_to_one")
    for validation, number, title, expected in [
        ("fixed4", 1, "Four-construct assessment", (4, 3)),
        ("LOCO", 2, "Scaffold cross-validation", (24, 22)),
    ]:
        frame = joined.loc[(joined.dataset == "PUF12_24") &
                           (joined.model == "RF_CP_structure9") &
                           (joined.validation == validation)].copy()
        counts = matrix_counts(frame.true_work, frame.prediction, [0, 1])
        assert len(frame) == expected[0] and np.trace(counts) == expected[1]
        metric = metrics.loc[(metrics.dataset == "PUF12_24") &
                             (metrics.model == "RF_CP_structure9") &
                             (metrics.validation == validation)].iloc[0]
        assert np.array_equal(counts, [[metric.TN, metric.FP], [metric.FN, metric.TP]])
        suffix = "four_construct" if number == 1 else "loco"
        stem = f"fig{number}_model1_{suffix}_confusion"
        frame.to_csv(DATA / f"{stem}_predictions.csv", index=False)
        pd.DataFrame(counts, index=["Non-work", "Work"], columns=["Non-work", "Work"]).to_csv(
            DATA / f"{stem}_counts.csv", index_label="Measured class")
        fig, ax = plt.subplots(figsize=(89 * MM, 78 * MM))
        fig.subplots_adjust(left=.26, right=.95, bottom=.24, top=.79)
        confusion_panel(ax, counts, ["Non-work", "Work"], title, vmax=int(counts.max()))
        fig.text(.60, .91, "Random forest · CP + structure", ha="center", fontsize=7)
        fig.text(.60, .06, f"{expected[1]}/{expected[0]} correct · n = {expected[0]} constructs",
                 ha="center", fontsize=7)
        save_figure(fig, stem, [ax])


def model2_figures():
    frame = read_source("five_construct_predictions")
    assert len(frame) == 15 and frame.target.nunique() == 5
    assert (frame.status == "ok").all() and not frame.isna().any().any()
    frame["matched_control_percent"] = frame.experimental_percent - frame.delta_pp
    frame["correct"] = frame.true_class == frame.predicted_class
    assert frame.correct.sum() == 13
    frame.to_csv(DATA / "model2_selected_five_complete.csv", index=False)
    definitions = [
        ("C295", ["decrease", "not_decrease"], ["Decrease", "No decrease"]),
        ("C388", ["decrease", "within_selected_interval", "increase"],
         ["Decrease", "Within\ninterval", "Increase"]),
        ("C871", ["decrease", "not_decrease"], ["Decrease", "No decrease"]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(183 * MM, 84 * MM))
    fig.subplots_adjust(left=.14, right=.985, bottom=.29, top=.75, wspace=.72)
    counts_long = []
    for i, (endpoint, labels, display) in enumerate(definitions):
        sub = frame.loc[frame.endpoint == endpoint]
        counts = matrix_counts(sub.true_class, sub.predicted_class, labels)
        assert len(sub) == 5
        confusion_panel(axes[i], counts, display,
                        f"{endpoint}  ·  {np.trace(counts)}/5 correct", vmax=5)
        panel_label(axes[i], chr(97+i))
        for a, actual in enumerate(labels):
            for p, prediction in enumerate(labels):
                counts_long.append({"endpoint": endpoint, "measured_class": actual,
                                    "predicted_class": prediction, "count": counts[a, p]})
    fig.text(.5, .93, "Selected five constructs: endpoint agreement", ha="center", fontsize=9)
    fig.text(.5, .07, "Each construct withheld individually · Fold-specific class thresholds",
             ha="center", fontsize=7)
    fig.align_xlabels(axes)
    pd.DataFrame(counts_long).to_csv(DATA / "fig3_model2_five_construct_confusion_counts.csv", index=False)
    save_figure(fig, "fig3_model2_five_construct_confusion", list(axes))

    order = frame.target.drop_duplicates().tolist()
    display_names = {"P4-R5-SNE+P7-R5-SNE": "P4/P7-SNE"}
    row_labels = [display_names.get(x, x) for x in order]
    columns = ["C295", "C388", "C871"]
    absolute = frame.pivot(index="target", columns="endpoint", values="experimental_percent").loc[order, columns]
    delta = frame.pivot(index="target", columns="endpoint", values="delta_pp").loc[order, columns]
    absolute.to_csv(DATA / "fig4_model2_measured_percent.csv")
    delta.to_csv(DATA / "fig4_model2_matched_delta_pp.csv")
    fig, axes = plt.subplots(1, 2, figsize=(183 * MM, 94 * MM))
    fig.subplots_adjust(left=.13, right=.985, bottom=.31, top=.79, wspace=.53)
    for idx, (values, title, cmap, norm, fmt, ticks, unit) in enumerate([
        (absolute.values, "Measured editing", "YlGnBu", Normalize(0, 100), ".1f", [0, 50, 100], "Editing (%)"),
        (delta.values, "Change from matched control", "RdBu_r", TwoSlopeNorm(vmin=-60, vcenter=0, vmax=60), "+.1f", [-60, 0, 60], "Change (percentage points)"),
    ]):
        ax = axes[idx]
        im = ax.imshow(values, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
        for row, col in np.ndindex(values.shape):
            value = values[row, col]
            rgba = im.cmap(norm(value))
            luminance = .2126*rgba[0] + .7152*rgba[1] + .0722*rgba[2]
            ax.text(col, row, format(value, fmt), ha="center", va="center", fontsize=7.5,
                    color="white" if luminance < .52 else DARK)
        ax.set_xticks(range(3), ["C295 (byst.)", "C388 (on)", "C871 (byst.)"])
        ax.set_yticks(range(5), row_labels)
        ax.set_title(title, pad=12)
        ax.tick_params(length=0, pad=5)
        for spine in ax.spines.values():
            spine.set_visible(False)
        panel_label(ax, chr(97 + idx))
        position = ax.get_position()
        cax = fig.add_axes([position.x0, .165, position.width, .025])
        cb = fig.colorbar(im, cax=cax, orientation="horizontal", ticks=ticks)
        cb.set_label(unit, labelpad=3)
        cb.outline.set_visible(False)
        cb.ax.tick_params(length=2)
    fig.text(.5, .95, "Measured editing across the five selected constructs", ha="center", fontsize=9)
    fig.text(.5, .025, "P4/P7-SNE: P4-R5-SNE + P7-R5-SNE · Construct-level values", ha="center", fontsize=6.5)
    save_figure(fig, "fig4_model2_measured_and_control_delta", list(axes))


def model3_figure():
    frame = read_source("mpnn_site_scores")
    selected = frame.loc[frame.consensus].copy()
    assert len(selected) == 17
    assert not selected.is_specificity_pos.any()
    assert (selected.ligandmpnn_mut == selected.proteinmpnn_mut).all()
    selected["mutation"] = selected.wt + selected.pos_1based.astype(str) + selected.ligandmpnn_mut
    selected["generated_sequences_per_model"] = 10000
    selected.to_csv(DATA / "fig5_model3_nontrm_consensus.csv", index=False)
    fig, ax = plt.subplots(figsize=(120 * MM, 126 * MM))
    fig.subplots_adjust(left=.18, right=.965, bottom=.14, top=.85)
    y = np.arange(len(selected))
    ligand = selected.ligandmpnn_f_max.to_numpy()
    protein = selected.proteinmpnn_f_max.to_numpy()
    ax.hlines(y, np.minimum(ligand, protein), np.maximum(ligand, protein), color=GREY, linewidth=1.3)
    ax.scatter(ligand, y, s=17, color=BLUE, marker="o", label="LigandMPNN", zorder=3)
    ax.scatter(protein, y, s=17, color=ORANGE, marker="s", label="ProteinMPNN", zorder=3)
    ax.set_yticks(y, selected.mutation)
    ax.set_ylim(len(selected)-.35, -.65)
    ax.set_xlim(0, 1.035)
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_xlabel("Frequency among generated sequences", labelpad=7)
    ax.tick_params(axis="y", length=0, pad=5)
    ax.grid(axis="x", color="#EAEFF2", linewidth=.5)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.legend(loc="lower center", bbox_to_anchor=(.5, 1.015), ncol=2, frameon=False,
              handlelength=1, columnspacing=1.7)
    save_figure(fig, "fig5_model3_nontrm_consensus", [ax])


def main():
    global require_matplotlib_panel_alignment
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qa-tools-dir", type=Path)
    args = parser.parse_args()
    if args.qa_tools_dir:
        sys.path.insert(0, str(args.qa_tools_dir))
        from audit_panel_alignment import require_matplotlib_panel_alignment
    for directory in (OUT, DATA, QA):
        directory.mkdir(parents=True, exist_ok=True)
    model1_figures()
    model2_figures()
    model3_figure()
    provenance = []
    for key, frame in SOURCE_FRAMES.items():
        path = ROOT / SOURCES[key]
        provenance.append({"key": key, "repository_path": SOURCES[key],
                           "rows": len(frame), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    (DATA / "source_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print("Wrote 5 figures (SVG, PDF, 600-dpi PNG), source CSVs and alignment reports.")


if __name__ == "__main__":
    main()
