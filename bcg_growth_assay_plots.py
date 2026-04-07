#!/usr/bin/env python3
"""
BCG Growth Assay: Bar charts of fold change (7dpi / 4hpi)
in integrated fluorescence (BCG-dsRed Count × Geometric Mean),
stratified by carprofen dose.

Grid layout: rows = MOI (B0, B1, B10, B100),
             columns = compartment (EC, IC).
X-axis: carprofen dose (0, 1, 10, 100 µg/mL).
Y-axis: Fold Change (Integrated Fluorescence, 7d/4h).
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = os.path.dirname(os.path.abspath(__file__))

# ── Constants ──────────────────────────────────────────────────────────────
MOIS = [0, 1, 10, 100]
MOI_LABELS = ["B0", "B1", "B10", "B100"]
CARPROFEN_DOSES = [0, 1, 10, 100]
CARPROFEN_LABELS = ["C0", "C1", "C10", "C100"]
COMPARTMENTS = ["EC", "IC"]
COMP_LABELS = {"EC": "Extracellular", "IC": "Intracellular"}
CARPROFEN_COLOURS = {0: "#2166AC", 1: "#4393C3", 10: "#92C5DE", 100: "#D1E5F0"}
DONOR_FILES = [f"AIF007-{i}_Table.xls" for i in range(1, 7)]


# ── Parse helpers ──────────────────────────────────────────────────────────
def parse_pct(val):
    if isinstance(val, str):
        return float(val.replace("%", "").strip())
    return float(val)


def parse_filename(fname):
    m = re.match(r"(\w+)_(EC|IC)_C(\d+)B(\d+)\.fcs", str(fname))
    if not m:
        return None
    return {
        "timepoint": m.group(1),
        "compartment": m.group(2),
        "carprofen": int(m.group(3)),
        "moi": int(m.group(4)),
    }


# ── Load all donors ───────────────────────────────────────────────────────
all_rows = []
for donor_idx, fname in enumerate(DONOR_FILES, start=1):
    fpath = os.path.join(BASE, fname)
    if not os.path.exists(fpath):
        print(f"  Warning: {fname} not found, skipping.")
        continue
    df = pd.read_excel(fpath)
    # Exclude Mean/SD summary rows
    df = df[~df.iloc[:, 0].isin(["Mean", "SD"])]

    for _, row in df.iterrows():
        parsed = parse_filename(row.iloc[0])
        if parsed is None:
            continue
        count = row["-Beads/BCG-dsRed | Count"]
        gmean = row["-Beads/BCG-dsRed | Geometric Mean (PE-A)"]
        try:
            count = float(count)
            gmean = float(gmean)
        except (ValueError, TypeError):
            continue
        integrated = count * gmean
        parsed["donor"] = f"D{donor_idx}"
        parsed["integrated_fluor"] = integrated
        all_rows.append(parsed)

data = pd.DataFrame(all_rows)
print(f"Loaded {len(data)} observations from {len(DONOR_FILES)} donors.")

# ── Compute FC = 7dpi / 4hpi for each (donor, compartment, carprofen, moi) ─
fc_rows = []
for donor in data["donor"].unique():
    for comp in COMPARTMENTS:
        for carp in CARPROFEN_DOSES:
            for moi in MOIS:
                mask_base = (
                    (data["donor"] == donor)
                    & (data["compartment"] == comp)
                    & (data["carprofen"] == carp)
                    & (data["moi"] == moi)
                )
                v4h = data[mask_base & (data["timepoint"] == "4hpi")]
                v7d = data[mask_base & (data["timepoint"] == "7dpi")]
                if len(v4h) == 1 and len(v7d) == 1:
                    val_4h = v4h["integrated_fluor"].values[0]
                    val_7d = v7d["integrated_fluor"].values[0]
                    if val_4h > 0:
                        fc = val_7d / val_4h
                    else:
                        fc = np.nan
                    fc_rows.append({
                        "donor": donor,
                        "compartment": comp,
                        "carprofen": carp,
                        "moi": moi,
                        "fc": fc,
                    })

fc_data = pd.DataFrame(fc_rows)
fc_data = fc_data.dropna(subset=["fc"])
print(f"Computed {len(fc_data)} fold-change values.")


# ═══════════════════════════════════════════════════════════════════════════
# Plot A: Facet grid — rows = MOI, columns = Compartment (bar charts)
# ═══════════════════════════════════════════════════════════════════════════

def make_facet_grid_bars():
    fig, axes = plt.subplots(
        len(MOIS), len(COMPARTMENTS),
        figsize=(8, 14),
        sharey="row", sharex=True,
    )

    bar_width = 0.6
    x_positions = np.arange(len(CARPROFEN_DOSES))

    for row_idx, moi in enumerate(MOIS):
        for col_idx, comp in enumerate(COMPARTMENTS):
            ax = axes[row_idx, col_idx]
            subset = fc_data[(fc_data["moi"] == moi) & (fc_data["compartment"] == comp)]

            means = []
            sems = []
            for carp in CARPROFEN_DOSES:
                grp = subset[subset["carprofen"] == carp]["fc"]
                means.append(grp.mean() if len(grp) > 0 else 0)
                sems.append(grp.sem() if len(grp) > 1 else 0)

            # Bars
            bars = ax.bar(
                x_positions, means, bar_width,
                yerr=sems, capsize=4,
                color=[CARPROFEN_COLOURS[d] for d in CARPROFEN_DOSES],
                edgecolor="black", linewidth=0.6,
                error_kw=dict(elinewidth=1.2, capthick=1.2, color="black"),
                zorder=2,
            )

            # Individual donor points
            for i, carp in enumerate(CARPROFEN_DOSES):
                grp = subset[subset["carprofen"] == carp]["fc"]
                if len(grp) > 0:
                    jitter = np.random.default_rng(42).uniform(-0.12, 0.12, size=len(grp))
                    ax.scatter(
                        np.full(len(grp), x_positions[i]) + jitter,
                        grp.values,
                        color="black", s=20, alpha=0.6, zorder=3,
                    )

            ax.set_xticks(x_positions)
            ax.set_xticklabels(CARPROFEN_LABELS, fontsize=10)
            ax.axhline(y=1, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=10)

            # Column titles (top row only)
            if row_idx == 0:
                ax.set_title(COMP_LABELS[comp], fontsize=13, fontweight="bold")
            # Row labels (left column only)
            if col_idx == 0:
                ax.set_ylabel(
                    f"FC (7d/4h)\n{MOI_LABELS[row_idx]}",
                    fontsize=11, fontweight="bold",
                )
            # X-axis label (bottom row only)
            if row_idx == len(MOIS) - 1:
                ax.set_xlabel("Carprofen dose (µg/mL)", fontsize=11)

    n_donors = fc_data["donor"].nunique()
    fig.suptitle(
        f"BCG Growth Assay — Fold Change (Integrated Fluorescence, 7d/4h)\n"
        f"Mean ± SEM | n={n_donors} donors",
        fontsize=14, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fname = "BCG_growth_FC_facet_grid.png"
    fig.savefig(os.path.join(BASE, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {fname}")


# ═══════════════════════════════════════════════════════════════════════════
# Plot B: Same grid but with continuous x-axis (carprofen as numeric)
# ═══════════════════════════════════════════════════════════════════════════

def make_facet_grid_continuous():
    fig, axes = plt.subplots(
        len(MOIS), len(COMPARTMENTS),
        figsize=(8, 14),
        sharey="row", sharex=True,
    )

    # Evenly-spaced x for carprofen doses
    CARP_POS = {0: 0, 1: 1, 10: 2, 100: 3}
    JITTER_SCALE = 0.08

    for row_idx, moi in enumerate(MOIS):
        for col_idx, comp in enumerate(COMPARTMENTS):
            ax = axes[row_idx, col_idx]
            subset = fc_data[(fc_data["moi"] == moi) & (fc_data["compartment"] == comp)]

            # Summary stats
            summary = subset.groupby("carprofen")["fc"].agg(["mean", "sem"])
            summary = summary.reindex(CARPROFEN_DOSES)
            valid = summary.dropna(subset=["mean"])

            if not valid.empty:
                x_vals = np.array([CARP_POS[d] for d in valid.index])
                ax.errorbar(
                    x_vals, valid["mean"], yerr=valid["sem"],
                    color="#2166AC", marker="o",
                    markersize=8, linewidth=2, capsize=5, capthick=1.5,
                    zorder=4,
                )

            # Individual donor points
            for carp in CARPROFEN_DOSES:
                grp = subset[subset["carprofen"] == carp]["fc"]
                if len(grp) > 0:
                    jitter = np.random.default_rng(42).uniform(
                        -JITTER_SCALE, JITTER_SCALE, size=len(grp)
                    )
                    ax.scatter(
                        np.full(len(grp), CARP_POS[carp]) + jitter,
                        grp.values,
                        color="#2166AC", s=25, alpha=0.4, zorder=3,
                    )

            ax.set_xticks(list(CARP_POS.values()))
            ax.set_xticklabels(CARPROFEN_LABELS, fontsize=10)
            ax.axhline(y=1, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=10)

            if row_idx == 0:
                ax.set_title(COMP_LABELS[comp], fontsize=13, fontweight="bold")
            if col_idx == 0:
                ax.set_ylabel(
                    f"FC (7d/4h)\n{MOI_LABELS[row_idx]}",
                    fontsize=11, fontweight="bold",
                )
            if row_idx == len(MOIS) - 1:
                ax.set_xlabel("Carprofen dose (µg/mL)", fontsize=11)

    n_donors = fc_data["donor"].nunique()
    fig.suptitle(
        f"BCG Growth Assay — Fold Change (Integrated Fluorescence, 7d/4h)\n"
        f"Mean ± SEM | n={n_donors} donors",
        fontsize=14, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fname = "BCG_growth_FC_continuous.png"
    fig.savefig(os.path.join(BASE, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {fname}")


# ═══════════════════════════════════════════════════════════════════════════
# Plot C: Individual panels per MOI (EC and IC side by side)
# ═══════════════════════════════════════════════════════════════════════════

def make_individual_moi_panels():
    for moi_idx, moi in enumerate(MOIS):
        fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
        bar_width = 0.6
        x_positions = np.arange(len(CARPROFEN_DOSES))

        for col_idx, comp in enumerate(COMPARTMENTS):
            ax = axes[col_idx]
            subset = fc_data[(fc_data["moi"] == moi) & (fc_data["compartment"] == comp)]

            means = []
            sems = []
            for carp in CARPROFEN_DOSES:
                grp = subset[subset["carprofen"] == carp]["fc"]
                means.append(grp.mean() if len(grp) > 0 else 0)
                sems.append(grp.sem() if len(grp) > 1 else 0)

            bars = ax.bar(
                x_positions, means, bar_width,
                yerr=sems, capsize=4,
                color=[CARPROFEN_COLOURS[d] for d in CARPROFEN_DOSES],
                edgecolor="black", linewidth=0.6,
                error_kw=dict(elinewidth=1.2, capthick=1.2, color="black"),
                zorder=2,
            )

            for i, carp in enumerate(CARPROFEN_DOSES):
                grp = subset[subset["carprofen"] == carp]["fc"]
                if len(grp) > 0:
                    jitter = np.random.default_rng(42).uniform(-0.12, 0.12, size=len(grp))
                    ax.scatter(
                        np.full(len(grp), x_positions[i]) + jitter,
                        grp.values,
                        color="black", s=25, alpha=0.6, zorder=3,
                    )

            ax.set_xticks(x_positions)
            ax.set_xticklabels(CARPROFEN_LABELS, fontsize=11)
            ax.set_xlabel("Carprofen dose (µg/mL)", fontsize=12)
            ax.set_title(COMP_LABELS[comp], fontsize=13, fontweight="bold")
            ax.axhline(y=1, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=11)

            if col_idx == 0:
                ax.set_ylabel("FC (Integrated Fluorescence, 7d/4h)",
                              fontsize=11, fontweight="bold")

        n_donors = fc_data["donor"].nunique()
        fig.suptitle(
            f"BCG Growth Assay — {MOI_LABELS[moi_idx]} | Mean ± SEM | n={n_donors}",
            fontsize=14, fontweight="bold",
        )
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        fname = f"BCG_growth_FC_{MOI_LABELS[moi_idx]}.png"
        fig.savefig(os.path.join(BASE, fname), dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {fname}")


# ── Run ───────────────────────────────────────────────────────────────────
print("\nGenerating plots...")
make_facet_grid_bars()
make_facet_grid_continuous()
make_individual_moi_panels()
print("\nDone.")
