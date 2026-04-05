"""
Generate Figures 2, 3, and 4 for carprofen/MDM report.

Figure 2: ELISA — LPS group (IL-1β, IL-6, TNF-α)
Figure 3: ELISA — BCG group (IL-1β, IL-6, TNF-α)
Figure 4: BCG growth assay — EC and IC, media-only vs +carprofen
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator
import xlrd

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'xtick.labelsize': 7.5,
    'ytick.labelsize': 7.5,
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'legend.frameon': False,
    'legend.fontsize': 7.5,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# Carprofen colours: 4 levels from light→dark teal/blue
CARP_COLORS  = ['#AECDE8', '#5BA4CF', '#2171B5', '#084594']
CARP_LABELS  = ['0 µg/mL', '1 µg/mL', '10 µg/mL', '100 µg/mL']
CARP_LEVELS  = ['C0', 'C1', 'C10', 'C100']

# ── Load ELISA data ──────────────────────────────────────────────────────────
def load_elisa():
    dfs = {}
    for sheet, col in [('IL1b', 'IL1b_pgml'), ('IL6', 'IL6_pgml'), ('TNFa', 'TNFa_pgml')]:
        import openpyxl
        wb = openpyxl.load_workbook('/home/user/elisa/ELISA.xlsx')
        ws = wb[sheet]
        rows = list(ws.iter_rows(values_only=True))
        header = rows[0]
        records = []
        for row in rows[1:]:
            records.append(dict(zip(header[:5], row[:5])))
        df = pd.DataFrame(records)
        df.columns = ['Donor', 'Carprofen', 'Stimulus', 'Stim_Conc', col]
        df[col] = pd.to_numeric(df[col], errors='coerce')
        dfs[sheet] = df
    return dfs

elisa = load_elisa()

# ── ELISA helper: mean±SEM per (carprofen, stim_conc) ───────────────────────
def elisa_stats(df, cytokine_col, stimulus, stim_levels):
    sub = df[(df['Stimulus'] == stimulus) & (df['Stim_Conc'].isin(stim_levels))].copy()
    stats = (sub.groupby(['Carprofen', 'Stim_Conc'])[cytokine_col]
               .agg(['mean', 'sem', lambda x: x.dropna().tolist()])
               .rename(columns={'<lambda_0>': 'vals'})
               .reset_index())
    return stats

# ── Generic ELISA figure (3-panel, line plot) ────────────────────────────────
def make_elisa_figure(stimulus, stim_levels, stim_xlabel, stim_ticklabels,
                      ylabels, titles, panel_letters, filename):
    """
    Line plot: x = stimulus concentration, one line per carprofen level.
    Solid circles = mean, error bars = SEM, small dots = individual donors.
    """
    cytokines = [('IL1b', 'IL1b_pgml'), ('IL6', 'IL6_pgml'), ('TNFa', 'TNFa_pgml')]
    n_stim    = len(stim_levels)
    x_pos     = np.arange(n_stim)          # evenly-spaced categorical positions

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8))
    fig.subplots_adjust(wspace=0.45)

    for ax, (sheet, col), ylabel, title, letter in zip(
            axes, cytokines, ylabels, titles, panel_letters):

        stats = elisa_stats(elisa[sheet], col, stimulus, stim_levels)

        for carp, color, clabel in zip(CARP_LEVELS, CARP_COLORS, CARP_LABELS):
            means, sems, all_vals = [], [], []
            for sl in stim_levels:
                row = stats[(stats['Carprofen'] == carp) & (stats['Stim_Conc'] == sl)]
                if len(row):
                    m = row['mean'].values[0]
                    s = row['sem'].values[0]
                    v = row['vals'].values[0]
                else:
                    m, s, v = np.nan, np.nan, []
                means.append(m if (m is not None and not (isinstance(m, float) and np.isnan(m))) else np.nan)
                sems.append(s if (s is not None and not (isinstance(s, float) and np.isnan(s))) else 0)
                all_vals.append(v)

            means = np.array(means, dtype=float)
            sems  = np.array(sems,  dtype=float)
            sems[np.isnan(sems)] = 0

            # Draw line + error bars (only where mean is not NaN)
            valid = ~np.isnan(means)
            if valid.any():
                ax.plot(x_pos[valid], means[valid], color=color, linewidth=1.4,
                        marker='o', markersize=5, markeredgewidth=0,
                        label=clabel, zorder=3)
                ax.errorbar(x_pos[valid], means[valid], yerr=sems[valid],
                            fmt='none', ecolor=color, elinewidth=1.0,
                            capsize=3, capthick=1.0, zorder=2)

            # Individual donor dots (small, slightly transparent)
            jitter_rng = np.random.default_rng(42)
            for xi, vals in zip(x_pos, all_vals):
                if vals:
                    jitter = jitter_rng.uniform(-0.12, 0.12, len(vals))
                    ax.scatter(xi + jitter, vals, s=9, color=color,
                               zorder=1, linewidths=0, alpha=0.55)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(stim_ticklabels, fontsize=7.5)
        ax.set_xlabel(stim_xlabel, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_title(title, fontsize=9, pad=4)
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=3))
        ax.set_xlim(-0.4, n_stim - 0.6)
        ax.set_ylim(bottom=0)

        ax.text(-0.18, 1.05, letter, transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='top', ha='left')

    # Single legend on last panel
    axes[2].legend(title='Carprofen', title_fontsize=7.5,
                   loc='upper left', fontsize=7.5,
                   handlelength=1.4, handleheight=0.9,
                   markerscale=0.9)

    fig.savefig(filename)
    plt.close(fig)
    print(f'Saved {filename}')


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: ELISA — LPS
# ─────────────────────────────────────────────────────────────────────────────
make_elisa_figure(
    stimulus='LPS',
    stim_levels=['L0', 'L1', 'L10', 'L100'],
    stim_xlabel='LPS concentration (ng/mL)',
    stim_ticklabels=['0', '1', '10', '100'],
    ylabels=['IL-1β (pg/mL)', 'IL-6 (pg/mL)', 'TNF-α (pg/mL)'],
    titles=['IL-1β', 'IL-6', 'TNF-α'],
    panel_letters=['A', 'B', 'C'],
    filename='/home/user/elisa/Figure2_ELISA_LPS.pdf'
)

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3: ELISA — BCG
# ─────────────────────────────────────────────────────────────────────────────
make_elisa_figure(
    stimulus='BCG',
    stim_levels=['B0', 'B1', 'B10', 'B100'],
    stim_xlabel='BCG (MOI)',
    stim_ticklabels=['0', '1', '10', '100'],
    ylabels=['IL-1β (pg/mL)', 'IL-6 (pg/mL)', 'TNF-α (pg/mL)'],
    titles=['IL-1β', 'IL-6', 'TNF-α'],
    panel_letters=['A', 'B', 'C'],
    filename='/home/user/elisa/Figure3_ELISA_BCG.pdf'
)

print('Figures 2 and 3 done.')

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4: BCG growth assay
# ─────────────────────────────────────────────────────────────────────────────

def load_bcg_data():
    """
    Returns a DataFrame with columns:
    condition, donor, timepoint, fraction, carprofen, moi, bcg_count, bead_count
    condition: 'media' (files 1-3) or 'carprofen' (files 4-6)
    """
    records = []
    file_map = {
        1: ('media', 'D1'), 2: ('media', 'D2'), 3: ('media', 'D3'),
        4: ('carprofen', 'D1'), 5: ('carprofen', 'D2'), 6: ('carprofen', 'D3'),
    }
    for fnum, (condition, donor) in file_map.items():
        fname = f'/home/user/elisa/AIF007-{fnum}_Table.xls'
        wb = xlrd.open_workbook(fname)
        sh = wb.sheet_by_name('Sheet0')
        for r in range(1, sh.nrows):
            row = sh.row_values(r)
            sample = str(row[0]).replace('.fcs', '')
            if not sample or 'Mean' in sample or 'SD' in sample:
                continue
            bead_count = row[2]
            bcg_count  = row[6]
            if not bead_count:
                continue
            # Parse sample name: e.g. 4hpi_EC_C0B10
            parts = sample.split('_')
            if len(parts) < 3:
                continue
            tp_str   = parts[0]       # 4hpi or 7dpi
            frac     = parts[1]       # EC or IC
            carp_moi = parts[2]       # C0B10 etc.
            # Extract carprofen and moi
            import re
            m = re.match(r'(C\d+)(B\d+)', carp_moi)
            if not m:
                continue
            carp = m.group(1)
            moi  = m.group(2)
            tp   = 4 if '4hpi' in tp_str else 7
            records.append(dict(condition=condition, donor=donor,
                                timepoint=tp, fraction=frac,
                                carprofen=carp, moi=moi,
                                bcg_count=float(bcg_count),
                                bead_count=float(bead_count)))
    return pd.DataFrame(records)

bcg_df = load_bcg_data()

def compute_fold_change(df):
    """Compute fold change = (7dpi_BCG/bead) / (4hpi_BCG/bead) per group."""
    df = df.copy()
    df['norm'] = df['bcg_count'] / df['bead_count']
    t4 = df[df['timepoint'] == 4].set_index(
        ['condition', 'donor', 'fraction', 'carprofen', 'moi'])['norm']
    t7 = df[df['timepoint'] == 7].set_index(
        ['condition', 'donor', 'fraction', 'carprofen', 'moi'])['norm']
    fc = (t7 / t4).reset_index()
    fc.columns = ['condition', 'donor', 'fraction', 'carprofen', 'moi', 'fold_change']
    return fc

fc_df = compute_fold_change(bcg_df)

def fc_stats(fc_df, fraction, condition, moi_levels, carp_levels):
    sub = fc_df[(fc_df['fraction'] == fraction) &
                (fc_df['condition'] == condition) &
                (fc_df['moi'].isin(moi_levels)) &
                (fc_df['carprofen'].isin(carp_levels))]
    stats = (sub.groupby(['carprofen', 'moi'])['fold_change']
               .agg(['mean', 'sem', list])
               .rename(columns={'list': 'vals'})
               .reset_index())
    return stats


# Figure 4 layout: 2 rows × 4 columns
# Row 0: EC  (panels A-D)
# Row 1: IC  (panels E-H)
# Columns: B0, B1, B10, B100

MOI_LEVELS  = ['B0', 'B1', 'B10', 'B100']
MOI_LABELS  = ['0', '1', '10', '100']

COND_COLORS = {'media': '#5B9BD5', 'carprofen': '#E06060'}
COND_LABELS = {'media': 'Media only', 'carprofen': '+Carprofen'}

# For Figure 4: x-axis = carprofen level, bars grouped by condition
n_carp   = len(CARP_LEVELS)
bar_w    = 0.35
x_centers = np.arange(n_carp)

def make_figure4(filename):
    fig, axes = plt.subplots(2, 4, figsize=(7.5, 4.8))
    fig.subplots_adjust(wspace=0.38, hspace=0.55)

    row_labels = ['EC', 'IC']
    fractions  = ['EC', 'IC']
    panel_letters_row = [['A', 'B', 'C', 'D'], ['E', 'F', 'G', 'H']]

    # Collect per-panel y limits (independent scaling per panel)
    panel_ymax = {}
    for ri, frac in enumerate(fractions):
        for ci, moi in enumerate(MOI_LEVELS):
            ymax = 0
            for condition in ['media', 'carprofen']:
                st = fc_stats(fc_df, frac, condition, [moi], CARP_LEVELS)
                for _, row in st.iterrows():
                    m = row['mean']
                    s = row['sem'] if not np.isnan(row['sem']) else 0
                    if not np.isnan(m):
                        ymax = max(ymax, m + s)
                        for v in row['vals']:
                            ymax = max(ymax, v)
            panel_ymax[(ri, ci)] = ymax * 1.18

    # Draw
    for ri, (frac, rlabel) in enumerate(zip(fractions, row_labels)):
        for ci, (moi, mlabel) in enumerate(zip(MOI_LEVELS, MOI_LABELS)):
            ax = axes[ri, ci]

            for cond_i, condition in enumerate(['media', 'carprofen']):
                color = COND_COLORS[condition]
                st = fc_stats(fc_df, frac, condition, [moi], CARP_LEVELS)
                offsets = x_centers + (cond_i - 0.5) * bar_w

                for xi, carp in zip(offsets, CARP_LEVELS):
                    row = st[st['carprofen'] == carp]
                    if len(row) == 0:
                        continue
                    m = row['mean'].values[0]
                    s = row['sem'].values[0] if not np.isnan(row['sem'].values[0]) else 0
                    vals = row['vals'].values[0]
                    if np.isnan(m):
                        continue

                    ax.bar(xi, m, width=bar_w * 0.85, color=color,
                           yerr=s,
                           error_kw={'elinewidth': 0.8, 'capsize': 2, 'capthick': 0.8},
                           zorder=2)

                    if vals:
                        jitter = np.random.default_rng(99).uniform(
                            -bar_w * 0.15, bar_w * 0.15, len(vals))
                        ax.scatter([xi + j for j in jitter], vals,
                                   s=10, color='black', zorder=3, linewidths=0, alpha=0.85)

            # Reference line at y=1
            ax.axhline(1, color='gray', linewidth=0.7, linestyle='--', zorder=1)

            ax.set_xticks(x_centers)
            ax.set_xticklabels(['0', '1', '10', '100'], fontsize=7)
            ax.set_ylim(0, panel_ymax[(ri, ci)])
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4, min_n_ticks=3))

            # Column titles (MOI) only on top row
            if ri == 0:
                ax.set_title(f'MOI {mlabel}', fontsize=8.5, pad=3)

            # Row label on the left
            if ci == 0:
                ax.set_ylabel(f'{rlabel}\nFold change\n(normalised to 4 hpi)', fontsize=7.5)
            else:
                ax.set_ylabel('')

            # x-axis label only on bottom row
            if ri == 1:
                ax.set_xlabel('Carprofen (µg/mL)', fontsize=7.5)

            # Panel letter
            letter = panel_letters_row[ri][ci]
            ax.text(-0.28, 1.08, letter, transform=ax.transAxes,
                    fontsize=10, fontweight='bold', va='top', ha='left')

    # Legend
    patches = [mpatches.Patch(color=COND_COLORS[c], label=COND_LABELS[c])
               for c in ['media', 'carprofen']]
    fig.legend(handles=patches, loc='lower center', ncol=2,
               bbox_to_anchor=(0.5, -0.02), fontsize=8, frameon=False)

    fig.savefig(filename)
    plt.close(fig)
    print(f'Saved {filename}')

make_figure4('/home/user/elisa/Figure4_BCG_growth_assay.pdf')

print('\nAll figures generated successfully.')
