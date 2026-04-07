#!/usr/bin/env python3
"""
Dose-response plots for PGE2, TNFα, IL-1β, and IL-6
stratified by carprofen dose.

X-axis: BCG dose (0, 1, 10, 100)
Y-axis: cytokine concentration (pg/mL)
Colour/line: carprofen dose (0, 1, 10, 100 µg/mL)

Each panel shows individual donor points, mean ± SEM, and
connecting lines through the group means.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit

# ── Paths ───────────────────────────────────────────────────────────────────
BASE = os.path.dirname(os.path.abspath(__file__))
ELISA_XLSX = os.path.join(BASE, "ELISA.xlsx")

PGE2_PLATE2_RAW = os.path.join(
    BASE, "20260401_AIF006_ELISA_PGE2_Plate2(Plate 1 - Sheet1).csv"
)
PGE2_PLATE3_RAW = os.path.join(
    BASE, "20260401_AIF006_ELISA_PGE2_Plate3(Plate 2 - Sheet1).csv"
)
PGE2_LAYOUT = os.path.join(BASE, "AIF006_ELISA_PGE2_Analysis(Layout).csv")


# ── Helper: parse coded factor levels into numeric ──────────────────────────
def parse_numeric(coded):
    """C0 → 0, C100 → 100, L10 → 10, B1 → 1, etc."""
    m = re.search(r"(\d+)", str(coded))
    return int(m.group(1)) if m else np.nan


# ═══════════════════════════════════════════════════════════════════════════
# 1. Load IL-1β, IL-6, TNFα from ELISA.xlsx
# ═══════════════════════════════════════════════════════════════════════════
def load_cytokine(sheet, value_col, display_name):
    """Read one sheet, return tidy df with numeric Carprofen/Stim_Conc."""
    df = pd.read_excel(ELISA_XLSX, sheet_name=sheet)
    df["Carprofen_num"] = df["Carprofen"].apply(parse_numeric)
    df["Stim_Conc_num"] = df["Stim_Conc"].apply(parse_numeric)
    df = df.rename(columns={value_col: "value"})
    df["cytokine"] = display_name
    return df[["Donor", "Carprofen_num", "Stimulus", "Stim_Conc_num", "value", "cytokine"]]


il1b = load_cytokine("IL1b", "IL1b_pgml", "IL-1β")
il6 = load_cytokine("IL6", "IL6_pgml", "IL-6")
tnfa = load_cytokine("TNFa", "TNFa_pgml", "TNFα")


# ═══════════════════════════════════════════════════════════════════════════
# 2. Process PGE2 from raw plate-reader OD (competitive ELISA)
# ═══════════════════════════════════════════════════════════════════════════

def read_plate_od(csv_path):
    """
    Parse Synergy H1 plate-reader CSV.  Return dict  {(row, col): OD_405}
    where row ∈ 'A'..'H' and col ∈ 1..12.
    We use 405 nm readings only.
    """
    rows_data = {}
    with open(csv_path) as fh:
        lines = fh.readlines()
    for line in lines:
        parts = line.strip().split(",")
        if len(parts) < 14:
            continue
        row_letter = parts[1].strip()
        wavelength = parts[-1].strip()
        if row_letter in "ABCDEFGH" and wavelength == "405":
            for col_idx in range(2, 14):
                col_num = col_idx - 1
                try:
                    od = float(parts[col_idx])
                except (ValueError, IndexError):
                    od = np.nan
                rows_data[(row_letter, col_num)] = od
    return rows_data


def parse_layout_block(layout_df, start_row):
    """
    Parse one plate block from the layout CSV (8 rows of data, A-H).
    """
    std_types = {}
    std_concs = {}
    sample_ids = {}
    for i in range(8):
        row = layout_df.iloc[start_row + i]
        row_letter = str(row.iloc[0]).strip()
        std_type = str(row.iloc[1]).strip()
        std_types[row_letter] = std_type
        for c in [2, 3]:
            try:
                conc = float(row.iloc[c])
                std_concs[(row_letter, c)] = conc
            except (ValueError, TypeError):
                pass
        for c in range(4, 13):
            val = str(row.iloc[c]).strip()
            if val and val != "nan" and val != "NaN":
                sample_ids[(row_letter, c)] = val
    return std_types, std_concs, sample_ids


def logistic4pl(x, bottom, top, ec50, slope):
    """4-parameter logistic function."""
    return bottom + (top - bottom) / (1 + (x / ec50) ** slope)


def fit_std_curve_and_predict(od_data, std_types, std_concs, sample_ids):
    """
    Competitive ELISA: high analyte → low OD.
    """
    nsb_ods = [od_data[(r, 1)] for r in "ABCDEFGH"
               if std_types.get(r) == "NSB" and (r, 1) in od_data]
    b0_ods = [od_data[(r, 1)] for r in "ABCDEFGH"
              if std_types.get(r) == "B0" and (r, 1) in od_data]

    nsb = np.mean(nsb_ods) if nsb_ods else 0
    b0 = np.mean(b0_ods)

    if b0 - nsb <= 0:
        return {}

    std_points = []
    for (row, col), conc in std_concs.items():
        if (row, col) in od_data and conc > 0:
            od = od_data[(row, col)]
            pct = (od - nsb) / (b0 - nsb) * 100
            std_points.append((conc, pct))

    if len(std_points) < 4:
        return {}

    std_points.sort()
    concs = np.array([p[0] for p in std_points])
    pcts = np.array([p[1] for p in std_points])

    try:
        popt, _ = curve_fit(
            logistic4pl, concs, pcts,
            p0=[0, 100, np.median(concs), -1],
            maxfev=10000,
        )
    except RuntimeError:
        return {}

    def predict_conc(pct_val):
        bottom, top, ec50, slope = popt
        if slope == 0:
            return np.nan
        ratio = (top - bottom) / (pct_val - bottom) - 1
        if ratio <= 0:
            return np.nan
        try:
            return ec50 * (ratio ** (1.0 / slope))
        except (ValueError, ZeroDivisionError):
            return np.nan

    results = {}
    for (row, col), sid in sample_ids.items():
        plate_col = col
        if (row, plate_col) in od_data:
            od = od_data[(row, plate_col)]
            pct = (od - nsb) / (b0 - nsb) * 100
            predicted = predict_conc(pct)
            results[sid] = predicted

    return results


def parse_sample_id(sid):
    """
    Parse 'Exp1 C0 B10 (1/1)' → (donor, carprofen, stimulus, stim_conc).
    """
    m = re.match(r"Exp(\d+)\s+C(\d+)\s+([BL])(\d+)", sid)
    if not m:
        return None
    exp_num = int(m.group(1))
    donor_map = {1: "D1", 2: "D2", 3: "D3", 4: "D4"}
    donor = donor_map.get(exp_num, f"D{exp_num}")
    carprofen = int(m.group(2))
    stim_letter = m.group(3)
    stimulus = "LPS" if stim_letter == "L" else "BCG"
    stim_conc = int(m.group(4))
    return donor, carprofen, stimulus, stim_conc


def process_pge2():
    """Process PGE2 from raw plate reader data using Plates 2 and 3."""
    layout = pd.read_csv(PGE2_LAYOUT, header=None)

    plate_configs = [
        {"raw_file": PGE2_PLATE2_RAW, "layout_start": 21, "name": "Plate2"},
        {"raw_file": PGE2_PLATE3_RAW, "layout_start": 31, "name": "Plate3"},
    ]

    all_results = []
    for cfg in plate_configs:
        od_data = read_plate_od(cfg["raw_file"])
        std_types, std_concs, sample_ids = parse_layout_block(
            layout, cfg["layout_start"]
        )

        predictions = fit_std_curve_and_predict(
            od_data, std_types, std_concs, sample_ids
        )

        for sid, conc in predictions.items():
            parsed = parse_sample_id(sid)
            if parsed is None:
                continue
            donor, carprofen, stimulus, stim_conc = parsed
            all_results.append({
                "Donor": donor,
                "Carprofen_num": carprofen,
                "Stimulus": stimulus,
                "Stim_Conc_num": stim_conc,
                "value": conc,
                "cytokine": "PGE2",
            })

    return pd.DataFrame(all_results)


pge2 = process_pge2()

# ═══════════════════════════════════════════════════════════════════════════
# 3. Combine all cytokines, filter to BCG only
# ═══════════════════════════════════════════════════════════════════════════
all_data = pd.concat([pge2, tnfa, il1b, il6], ignore_index=True)
bcg = all_data[all_data["Stimulus"] == "BCG"].copy()
bcg = bcg.dropna(subset=["value"])

# ═══════════════════════════════════════════════════════════════════════════
# 4. Create the plots
# ═══════════════════════════════════════════════════════════════════════════

CYTOKINE_ORDER = ["PGE2", "TNFα", "IL-1β", "IL-6"]
CARPROFEN_DOSES = [0, 1, 10, 100]
BCG_DOSES = [0, 1, 10, 100]
# Map real BCG dose to evenly-spaced x position
BCG_POS = {0: 0, 1: 1, 10: 2, 100: 3}
BCG_LABELS = ["0", "1", "10", "100"]
COLOURS = {0: "#2166AC", 1: "#4DAC26", 10: "#F4A582", 100: "#D6604D"}
MARKERS = {0: "o", 1: "s", 10: "^", 100: "D"}

# Small jitter offsets for each carprofen dose so points don't overlap
JITTER = {0: -0.12, 1: -0.04, 10: 0.04, 100: 0.12}


def make_plot_single_panel():
    """
    Option A: 4 one-panel figures, one per cytokine.
    Each panel has BCG dose on x, cytokine on y, stratified by carprofen.
    """
    for cyto in CYTOKINE_ORDER:
        subset = bcg[bcg["cytokine"] == cyto]
        if subset.empty:
            print(f"  ⚠ No data for {cyto} (BCG), skipping.")
            continue

        fig, ax = plt.subplots(figsize=(6, 5))

        for cdose in CARPROFEN_DOSES:
            grp = subset[subset["Carprofen_num"] == cdose]
            if grp.empty:
                continue

            # Individual donor points with jitter (use evenly-spaced x)
            x_pos = grp["Stim_Conc_num"].map(BCG_POS).values + JITTER[cdose]
            ax.scatter(
                x_pos, grp["value"],
                color=COLOURS[cdose], marker=MARKERS[cdose],
                s=40, alpha=0.5, edgecolors="none", zorder=3,
            )

            # Mean ± SEM
            summary = grp.groupby("Stim_Conc_num")["value"].agg(["mean", "sem", "count"])
            summary = summary.reindex(BCG_DOSES)
            valid = summary.dropna(subset=["mean"])
            if valid.empty:
                continue

            x_mean = np.array([BCG_POS[d] for d in valid.index]) + JITTER[cdose]
            ax.errorbar(
                x_mean, valid["mean"], yerr=valid["sem"],
                color=COLOURS[cdose], marker=MARKERS[cdose],
                markersize=8, linewidth=1.8, capsize=4, capthick=1.5,
                label=f"Carprofen {cdose} µg/mL",
                zorder=4,
            )

        ax.set_xticks(list(BCG_POS.values()))
        ax.set_xticklabels(BCG_LABELS)
        ax.set_xlabel("BCG dose", fontsize=13, fontweight="bold")

        unit = "pg/mL"
        ax.set_ylabel(f"{cyto} ({unit})", fontsize=13, fontweight="bold")
        ax.set_title(
            f"{cyto} dose-response by carprofen dose",
            fontsize=14, fontweight="bold",
        )
        ax.legend(fontsize=9, frameon=True, loc="best")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=11)

        fig.tight_layout()
        fname = f"BCG_dose_response_{cyto.replace('-', '').replace('α', 'a').replace('β', 'b')}_by_carprofen.png"
        fig.savefig(os.path.join(BASE, fname), dpi=300)
        plt.close(fig)
        print(f"  Saved {fname}")


def make_plot_grid():
    """
    Option B: 2×2 grid with all four cytokines in one figure.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for idx, cyto in enumerate(CYTOKINE_ORDER):
        ax = axes[idx]
        subset = bcg[bcg["cytokine"] == cyto]

        for cdose in CARPROFEN_DOSES:
            grp = subset[subset["Carprofen_num"] == cdose]
            if grp.empty:
                continue

            # Individual donor points (evenly-spaced x)
            x_pos = grp["Stim_Conc_num"].map(BCG_POS).values + JITTER[cdose]
            ax.scatter(
                x_pos, grp["value"],
                color=COLOURS[cdose], marker=MARKERS[cdose],
                s=30, alpha=0.45, edgecolors="none", zorder=3,
            )

            # Mean ± SEM with connecting line
            summary = grp.groupby("Stim_Conc_num")["value"].agg(["mean", "sem"])
            summary = summary.reindex(BCG_DOSES)
            valid = summary.dropna(subset=["mean"])
            if valid.empty:
                continue

            x_mean = np.array([BCG_POS[d] for d in valid.index]) + JITTER[cdose]
            ax.errorbar(
                x_mean, valid["mean"], yerr=valid["sem"],
                color=COLOURS[cdose], marker=MARKERS[cdose],
                markersize=7, linewidth=1.6, capsize=3, capthick=1.2,
                label=f"{cdose}" if idx == 0 else None,
                zorder=4,
            )

        ax.set_xticks(list(BCG_POS.values()))
        ax.set_xticklabels(BCG_LABELS)
        ax.set_xlabel("BCG dose", fontsize=11)
        ax.set_ylabel(f"{cyto} (pg/mL)", fontsize=11, fontweight="bold")
        ax.set_title(cyto, fontsize=13, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=10)

    # Shared legend
    legend_handles = [
        Line2D([0], [0], color=COLOURS[d], marker=MARKERS[d], markersize=7,
               linewidth=1.6, label=f"Carprofen {d} µg/mL")
        for d in CARPROFEN_DOSES
    ]
    fig.legend(
        handles=legend_handles, loc="lower center",
        ncol=4, fontsize=10, frameon=True,
        bbox_to_anchor=(0.5, -0.02),
    )

    fig.suptitle(
        "Cytokine dose-response to BCG, stratified by carprofen dose",
        fontsize=15, fontweight="bold", y=1.01,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.99])
    fname = "BCG_dose_response_all_cytokines_by_carprofen.png"
    fig.savefig(os.path.join(BASE, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {fname}")


def make_plot_facet_grid():
    """
    Option C: Grid of panels — one row per cytokine, one column per
    carprofen dose. This gives the clearest separation.
    """
    fig, axes = plt.subplots(
        4, 4, figsize=(16, 14),
        sharey="row", sharex=True,
    )

    for row_idx, cyto in enumerate(CYTOKINE_ORDER):
        subset = bcg[bcg["cytokine"] == cyto]
        for col_idx, cdose in enumerate(CARPROFEN_DOSES):
            ax = axes[row_idx, col_idx]
            grp = subset[subset["Carprofen_num"] == cdose]

            if not grp.empty:
                # Individual donor points (evenly-spaced x)
                for donor in grp["Donor"].unique():
                    ddata = grp[grp["Donor"] == donor]
                    x_pos = ddata["Stim_Conc_num"].map(BCG_POS).values
                    ax.scatter(
                        x_pos, ddata["value"],
                        color=COLOURS[cdose], s=35, alpha=0.5,
                        edgecolors="none", zorder=3,
                    )

                # Mean ± SEM
                summary = grp.groupby("Stim_Conc_num")["value"].agg(["mean", "sem"])
                summary = summary.reindex(BCG_DOSES)
                valid = summary.dropna(subset=["mean"])
                if not valid.empty:
                    x_mean = np.array([BCG_POS[d] for d in valid.index])
                    ax.errorbar(
                        x_mean, valid["mean"], yerr=valid["sem"],
                        color=COLOURS[cdose], marker="o",
                        markersize=7, linewidth=2, capsize=4, capthick=1.5,
                        zorder=4,
                    )

            ax.set_xticks(list(BCG_POS.values()))
            ax.set_xticklabels(BCG_LABELS, fontsize=9)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            if row_idx == 0:
                ax.set_title(
                    f"Carprofen {cdose} µg/mL",
                    fontsize=12, fontweight="bold", color=COLOURS[cdose],
                )
            if row_idx == 3:
                ax.set_xlabel("BCG dose", fontsize=11)
            if col_idx == 0:
                ax.set_ylabel(f"{cyto} (pg/mL)", fontsize=11, fontweight="bold")

    fig.suptitle(
        "Cytokine dose-response to BCG, stratified by carprofen dose",
        fontsize=16, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fname = "BCG_dose_response_facet_grid.png"
    fig.savefig(os.path.join(BASE, fname), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {fname}")


# ── Run all three plot variants ─────────────────────────────────────────
print("\nData summary (BCG only, non-NA):")
for cyto in CYTOKINE_ORDER:
    n = bcg[bcg["cytokine"] == cyto].shape[0]
    print(f"  {cyto}: {n} observations")

print("\nGenerating plots...")
make_plot_single_panel()
make_plot_grid()
make_plot_facet_grid()
print("\nDone.")
