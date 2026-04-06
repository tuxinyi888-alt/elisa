###############################################
# IL-1β ELISA Figures – LPS & BCG
# Python version matching R script output
###############################################

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import FancyBboxPatch

# ── Config ──────────────────────────────────
base_dir = os.path.dirname(os.path.abspath(__file__))
export_dir = os.path.join(base_dir, "AIF006_IL1b_Figures")
os.makedirs(export_dir, exist_ok=True)

# ── Import data ─────────────────────────────
lps = pd.read_csv(os.path.join(base_dir, "AIF006_ELISA_IL1b(LPS).csv"))
bcg = pd.read_csv(os.path.join(base_dir, "AIF006_ELISA_IL1b(BCG).csv"))
all_data = pd.concat([lps, bcg], ignore_index=True)

# ── Colour palette matching R script ────────
carprofen_levels = [0, 1, 10, 100]
colors = {0: "#4DAF4A", 1: "#377EB8", 10: "#FF7F00", 100: "#E41A1C"}

stim_levels = [1, 10, 100]  # drop Stim_Conc == 0 (all NA)

# ── Summary statistics (mean ± SEM) ────────
def compute_summary(df):
    rows = []
    for carp in carprofen_levels:
        for stim in stim_levels:
            sub = df[(df["Carprofen"] == carp) & (df["Stim_Conc"] == stim)]["IL1b_pgml"].dropna()
            n = len(sub)
            if n > 0:
                mean_val = sub.mean()
                sem_val = sub.std(ddof=1) / np.sqrt(n) if n > 1 else 0
            else:
                mean_val = 0
                sem_val = 0
            rows.append({"Carprofen": carp, "Stim_Conc": stim,
                         "mean": mean_val, "sem": sem_val, "n": n})
    return pd.DataFrame(rows)


# ── Shared theme helper ─────────────────────
def apply_theme(ax, title, xlabel):
    """Apply the R theme_minimal style to an axes."""
    # Title and labels
    ax.set_title(title, fontsize=18, fontweight="bold", pad=10)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel("IL-1β (pg/mL)", fontsize=18)

    # Axis tick labels – bold, size 14
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(14)
        label.set_fontweight("bold")

    # Remove grid (theme_minimal + panel.grid.major/minor = element_blank)
    ax.grid(False)

    # Black border (panel.border)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_color("black")

    # Tick marks inward, length 2 mm ≈ 5.67 pt
    ax.tick_params(axis="both", which="both", direction="out",
                   length=5.67, width=1, color="black")

    # y-axis starts at 0
    ax.set_ylim(bottom=0)


# ── Plotting function ──────────────────────
def make_il1b_plot(stimulus_name, ax=None, standalone=True):
    df = all_data[all_data["Stimulus"] == stimulus_name].copy()
    summary = compute_summary(df)

    if ax is None:
        fig, ax = plt.subplots(figsize=(18 / 2.54, 14 / 2.54))  # cm → inches
    else:
        fig = ax.figure

    n_groups = len(stim_levels)       # 3 stim conc groups
    n_bars = len(carprofen_levels)    # 4 carprofen levels
    group_width = 0.8
    bar_width = group_width / n_bars

    x_positions = np.arange(n_groups)

    for i, carp in enumerate(carprofen_levels):
        sub = summary[summary["Carprofen"] == carp]
        offsets = x_positions - group_width / 2 + bar_width * (i + 0.5)

        means = sub["mean"].values
        sems = sub["sem"].values

        # Bar
        ax.bar(offsets, means, width=bar_width * 0.875,
               color=colors[carp], edgecolor="black", linewidth=0.8,
               label=str(carp), zorder=2)

        # Error bars
        ax.errorbar(offsets, means, yerr=sems,
                     fmt="none", ecolor="black", elinewidth=1,
                     capsize=3, capthick=1, zorder=3)

        # Individual donor points
        for j, stim in enumerate(stim_levels):
            pts = df[(df["Carprofen"] == carp) & (df["Stim_Conc"] == stim)]["IL1b_pgml"].dropna()
            if len(pts) > 0:
                jitter = np.zeros(len(pts))  # no jitter needed, centred on bar
                ax.scatter(np.full(len(pts), offsets[j]) + jitter, pts.values,
                           s=25, facecolors="white", edgecolors="black",
                           linewidths=0.8, zorder=4)

    ax.set_xticks(x_positions)
    ax.set_xticklabels([str(s) for s in stim_levels])

    title = f"IL-1\u03B2 ({stimulus_name})"
    xlabel = f"{stimulus_name} concentration (ng/mL)"
    apply_theme(ax, title, xlabel)

    # Legend
    legend = ax.legend(title="Carprofen (\u03BCM)", fontsize=14,
                       title_fontsize=14, frameon=False,
                       loc="best")

    if standalone:
        fig.tight_layout()

    return fig, ax


# ── Generate individual figures ─────────────
fig_lps, _ = make_il1b_plot("LPS")
fig_lps.savefig(os.path.join(export_dir, "AIF006_IL1b_LPS.png"), dpi=300, bbox_inches="tight")
fig_lps.savefig(os.path.join(export_dir, "AIF006_IL1b_LPS.pdf"), bbox_inches="tight")
print("Saved LPS figure")

fig_bcg, _ = make_il1b_plot("BCG")
fig_bcg.savefig(os.path.join(export_dir, "AIF006_IL1b_BCG.png"), dpi=300, bbox_inches="tight")
fig_bcg.savefig(os.path.join(export_dir, "AIF006_IL1b_BCG.pdf"), bbox_inches="tight")
print("Saved BCG figure")

# ── Combined panel (A / B) ──────────────────
fig_comb, axes = plt.subplots(1, 2, figsize=(34 / 2.54, 16 / 2.54))

make_il1b_plot("LPS", ax=axes[0], standalone=False)
make_il1b_plot("BCG", ax=axes[1], standalone=False)

# Panel labels
axes[0].text(-0.12, 1.08, "A", transform=axes[0].transAxes,
             fontsize=20, fontweight="bold", va="top")
axes[1].text(-0.12, 1.08, "B", transform=axes[1].transAxes,
             fontsize=20, fontweight="bold", va="top")

fig_comb.tight_layout()
fig_comb.savefig(os.path.join(export_dir, "AIF006_IL1b_Combined.png"), dpi=300, bbox_inches="tight")
fig_comb.savefig(os.path.join(export_dir, "AIF006_IL1b_Combined.pdf"), bbox_inches="tight")
print("Saved Combined figure")

print(f"\nAll figures exported to: {export_dir}")
