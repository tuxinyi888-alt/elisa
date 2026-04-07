#!/usr/bin/env python3
"""
BCG Growth Assay plots, split into Exp 1–3 and Exp 4–6 (n=3 each).

Two plot types:
  1. 4h baseline: Absolute integrated fluorescence at 4hpi for ALL MOIs
     (B0, B1, B10, B100) — confirms infection establishment.
     Grid: rows = MOI, columns = compartment (EC, IC).

  2. 7d/4h Fold Change: Only for INFECTED cultures (B1, B10, B100).
     B0 excluded — no biological meaning for FC in uninfected wells.
     Grid: rows = MOI (B1, B10, B100), columns = compartment (EC, IC).
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))

# ── Constants ──────────────────────────────────────────────────────────────
ALL_MOIS = [0, 1, 10, 100]
ALL_MOI_LABELS = {0: "B0", 1: "B1", 10: "B10", 100: "B100"}
INFECTED_MOIS = [1, 10, 100]  # exclude B0 for FC
CARPROFEN_DOSES = [0, 1, 10, 100]
CARPROFEN_LABELS = ["C0", "C1", "C10", "C100"]
COMPARTMENTS = ["EC", "IC"]
COMP_LABELS = {"EC": "Extracellular", "IC": "Intracellular"}
CARPROFEN_COLOURS = {0: "#2166AC", 1: "#4393C3", 10: "#92C5DE", 100: "#D1E5F0"}

GROUPS = {
    "Exp1-3": {"files": [1, 2, 3], "label": "Exp 1\u20133"},
    "Exp4-6": {"files": [4, 5, 6], "label": "Exp 4\u20136"},
    "All":    {"files": [1, 2, 3, 4, 5, 6], "label": "All Experiments"},
}


# ── Parse helpers ──────────────────────────────────────────────────────────
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
for donor_idx in range(1, 7):
    fpath = os.path.join(BASE, f"AIF007-{donor_idx}_Table.xls")
    if not os.path.exists(fpath):
        print(f"  Warning: AIF007-{donor_idx}_Table.xls not found, skipping.")
        continue
    df = pd.read_excel(fpath)
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
        parsed["donor"] = f"D{donor_idx}"
        parsed["exp_num"] = donor_idx
        parsed["integrated_fluor"] = count * gmean
        all_rows.append(parsed)

data = pd.DataFrame(all_rows)
print(f"Loaded {len(data)} observations from {data['donor'].nunique()} donors.")

# ── Compute FC = 7dpi / 4hpi (all MOIs, filter later) ────────────────────
fc_rows = []
for donor in data["donor"].unique():
    exp_num = data[data["donor"] == donor]["exp_num"].iloc[0]
    for comp in COMPARTMENTS:
        for carp in CARPROFEN_DOSES:
            for moi in ALL_MOIS:
                mask = (
                    (data["donor"] == donor)
                    & (data["compartment"] == comp)
                    & (data["carprofen"] == carp)
                    & (data["moi"] == moi)
                )
                v4h = data[mask & (data["timepoint"] == "4hpi")]
                v7d = data[mask & (data["timepoint"] == "7dpi")]
                if len(v4h) == 1 and len(v7d) == 1:
                    val_4h = v4h["integrated_fluor"].values[0]
                    val_7d = v7d["integrated_fluor"].values[0]
                    fc = val_7d / val_4h if val_4h > 0 else np.nan
                    fc_rows.append({
                        "donor": donor,
                        "exp_num": exp_num,
                        "compartment": comp,
                        "carprofen": carp,
                        "moi": moi,
                        "fc": fc,
                    })

fc_all = pd.DataFrame(fc_rows).dropna(subset=["fc"])
print(f"Computed {len(fc_all)} fold-change values total.")


# ═══════════════════════════════════════════════════════════════════════════
# PLOT 1: 4h Baseline — Absolute integrated fluorescence at 4hpi
#         All MOIs (B0, B1, B10, B100), rows=MOI, cols=compartment
# ═══════════════════════════════════════════════════════════════════════════

def make_baseline_4h(raw_data, group_label, file_suffix):
    """
    Bar chart of absolute integrated fluorescence at 4hpi.
    Rows = MOI (B0, B1, B10, B100), Columns = Compartment (EC, IC).
    """
    baseline = raw_data[raw_data["timepoint"] == "4hpi"].copy()
    n_donors = baseline["donor"].nunique()
    moi_list = ALL_MOIS
    moi_labels = [ALL_MOI_LABELS[m] for m in moi_list]

    fig, axes = plt.subplots(
        len(moi_list), len(COMPARTMENTS),
        figsize=(8, 14),
        sharey="row", sharex=True,
    )

    bar_width = 0.6
    x_positions = np.arange(len(CARPROFEN_DOSES))

    for row_idx, moi in enumerate(moi_list):
        for col_idx, comp in enumerate(COMPARTMENTS):
            ax = axes[row_idx, col_idx]
            subset = baseline[
                (baseline["moi"] == moi) & (baseline["compartment"] == comp)
            ]

            means, sems = [], []
            for carp in CARPROFEN_DOSES:
                grp = subset[subset["carprofen"] == carp]["integrated_fluor"]
                means.append(grp.mean() if len(grp) > 0 else 0)
                sems.append(grp.sem() if len(grp) > 1 else 0)

            ax.bar(
                x_positions, means, bar_width,
                yerr=sems, capsize=4,
                color=[CARPROFEN_COLOURS[d] for d in CARPROFEN_DOSES],
                edgecolor="black", linewidth=0.6,
                error_kw=dict(elinewidth=1.2, capthick=1.2, color="black"),
                zorder=2,
            )

            # Individual donor points
            for i, carp in enumerate(CARPROFEN_DOSES):
                grp = subset[subset["carprofen"] == carp]["integrated_fluor"]
                if len(grp) > 0:
                    jitter = np.random.default_rng(42).uniform(-0.12, 0.12, size=len(grp))
                    ax.scatter(
                        np.full(len(grp), x_positions[i]) + jitter,
                        grp.values,
                        color="black", s=20, alpha=0.6, zorder=3,
                    )

            ax.set_xticks(x_positions)
            ax.set_xticklabels(CARPROFEN_LABELS, fontsize=10)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=10)

            if row_idx == 0:
                ax.set_title(COMP_LABELS[comp], fontsize=13, fontweight="bold")
            if col_idx == 0:
                ax.set_ylabel(
                    f"Integrated Fluorescence\n{moi_labels[row_idx]}",
                    fontsize=10, fontweight="bold",
                )
            if row_idx == len(moi_list) - 1:
                ax.set_xlabel("Carprofen dose (\u00b5g/mL)", fontsize=11)

    fig.suptitle(
        f"BCG Growth Assay \u2014 {group_label}\n"
        f"4h Baseline (Integrated Fluorescence) | Mean \u00b1 SEM | n={n_donors}",
        fontsize=14, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fname = f"BCG_baseline_4h_{file_suffix}.png"
    fig.savefig(os.path.join(BASE, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {fname}")


# ═══════════════════════════════════════════════════════════════════════════
# PLOT 2: 7d/4h Fold Change — INFECTED ONLY (B1, B10, B100)
#         Rows = MOI (B1, B10, B100), Columns = Compartment (EC, IC)
# ═══════════════════════════════════════════════════════════════════════════

def make_fc_infected_bars(fc_data, group_label, file_suffix):
    """3×2 facet grid bar charts for infected MOIs only."""
    moi_list = INFECTED_MOIS
    moi_labels = [ALL_MOI_LABELS[m] for m in moi_list]
    n_donors = fc_data["donor"].nunique()

    fig, axes = plt.subplots(
        len(moi_list), len(COMPARTMENTS),
        figsize=(8, 11),
        sharey="row", sharex=True,
    )

    bar_width = 0.6
    x_positions = np.arange(len(CARPROFEN_DOSES))

    for row_idx, moi in enumerate(moi_list):
        for col_idx, comp in enumerate(COMPARTMENTS):
            ax = axes[row_idx, col_idx]
            subset = fc_data[
                (fc_data["moi"] == moi) & (fc_data["compartment"] == comp)
            ]

            means, sems = [], []
            for carp in CARPROFEN_DOSES:
                grp = subset[subset["carprofen"] == carp]["fc"]
                means.append(grp.mean() if len(grp) > 0 else 0)
                sems.append(grp.sem() if len(grp) > 1 else 0)

            ax.bar(
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
                        color="black", s=20, alpha=0.6, zorder=3,
                    )

            ax.set_xticks(x_positions)
            ax.set_xticklabels(CARPROFEN_LABELS, fontsize=10)
            ax.axhline(y=1, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=10)

            if row_idx == 0:
                ax.set_title(COMP_LABELS[comp], fontsize=13, fontweight="bold")
            if col_idx == 0:
                ax.set_ylabel(
                    f"FC (7d/4h)\n{moi_labels[row_idx]}",
                    fontsize=11, fontweight="bold",
                )
            if row_idx == len(moi_list) - 1:
                ax.set_xlabel("Carprofen dose (\u00b5g/mL)", fontsize=11)

    fig.suptitle(
        f"BCG Growth Assay \u2014 {group_label}\n"
        f"Fold Change (Integrated Fluorescence, 7d/4h) | Mean \u00b1 SEM | n={n_donors}\n"
        f"Infected cultures only (B0 excluded)",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fname = f"BCG_growth_FC_infected_{file_suffix}.png"
    fig.savefig(os.path.join(BASE, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {fname}")


def make_fc_infected_continuous(fc_data, group_label, file_suffix):
    """3×2 facet grid with connected mean ± SEM for infected MOIs only."""
    moi_list = INFECTED_MOIS
    moi_labels = [ALL_MOI_LABELS[m] for m in moi_list]
    n_donors = fc_data["donor"].nunique()

    fig, axes = plt.subplots(
        len(moi_list), len(COMPARTMENTS),
        figsize=(8, 11),
        sharey="row", sharex=True,
    )

    CARP_POS = {0: 0, 1: 1, 10: 2, 100: 3}
    JITTER_SCALE = 0.08

    for row_idx, moi in enumerate(moi_list):
        for col_idx, comp in enumerate(COMPARTMENTS):
            ax = axes[row_idx, col_idx]
            subset = fc_data[
                (fc_data["moi"] == moi) & (fc_data["compartment"] == comp)
            ]

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
                    f"FC (7d/4h)\n{moi_labels[row_idx]}",
                    fontsize=11, fontweight="bold",
                )
            if row_idx == len(moi_list) - 1:
                ax.set_xlabel("Carprofen dose (\u00b5g/mL)", fontsize=11)

    fig.suptitle(
        f"BCG Growth Assay \u2014 {group_label}\n"
        f"Fold Change (Integrated Fluorescence, 7d/4h) | Mean \u00b1 SEM | n={n_donors}\n"
        f"Infected cultures only (B0 excluded)",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fname = f"BCG_growth_FC_infected_continuous_{file_suffix}.png"
    fig.savefig(os.path.join(BASE, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {fname}")


# ═══════════════════════════════════════════════════════════════════════════
# Run for each group
# ═══════════════════════════════════════════════════════════════════════════
for grp_key, grp_cfg in GROUPS.items():
    exp_nums = grp_cfg["files"]
    label = grp_cfg["label"]

    # Raw data subset for this group
    raw_subset = data[data["exp_num"].isin(exp_nums)]
    # FC subset — infected only
    fc_infected = fc_all[
        (fc_all["exp_num"].isin(exp_nums)) & (fc_all["moi"].isin(INFECTED_MOIS))
    ]

    donors_in = sorted(raw_subset["donor"].unique())
    n = len(donors_in)
    print(f"\n{'='*60}")
    print(f"  {label} — donors: {', '.join(donors_in)} (n={n})")
    print(f"  4h baseline obs: {len(raw_subset[raw_subset['timepoint']=='4hpi'])}")
    print(f"  FC values (infected): {len(fc_infected)}")
    print(f"{'='*60}")

    print(f"\nGenerating plots for {label}...")
    # 1. 4h baseline (all MOIs including B0)
    make_baseline_4h(raw_subset, label, grp_key)
    # 2. FC plots (infected only: B1, B10, B100)
    make_fc_infected_bars(fc_infected, label, grp_key)
    make_fc_infected_continuous(fc_infected, label, grp_key)

print("\nDone.")
