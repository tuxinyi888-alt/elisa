"""
Generate Figures 2, 3, 4 — publication-quality, R/ggplot-equivalent style.
Uses Python / matplotlib + scipy / statsmodels.
"""

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.ticker import LogLocator, NullFormatter
import openpyxl, xlrd, re, warnings
warnings.filterwarnings('ignore')

# ── rcParams ─────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':          'sans-serif',
    'font.sans-serif':      ['Arial','DejaVu Sans'],
    'font.size':            7,
    'axes.titlesize':       7.5,
    'axes.labelsize':       7,
    'xtick.labelsize':      6.5,
    'ytick.labelsize':      6.5,
    'axes.linewidth':       0.5,
    'xtick.major.width':    0.5,
    'ytick.major.width':    0.5,
    'xtick.major.size':     2.5,
    'ytick.major.size':     2.5,
    'axes.spines.top':      False,
    'axes.spines.right':    False,
    'legend.frameon':       False,
    'legend.fontsize':      6.5,
    'figure.dpi':           150,
    'savefig.dpi':          300,
    'savefig.bbox':         'tight',
})

# ── Colour palettes ───────────────────────────────────────────────────────────
LPS_PAL = {0: '#D9D9D9', 1: '#FDBB84', 10: '#E34A33', 100: '#7F0000'}
BCG_PAL = {0: '#D9D9D9', 1: '#9ECAE1', 10: '#3182BD', 100: '#08306B'}

# ─────────────────────────────────────────────────────────────────────────────
# ELISA data
# ─────────────────────────────────────────────────────────────────────────────
def load_elisa():
    records = []
    for sheet, col in [('IL1b','IL1b_pgml'),('IL6','IL6_pgml'),('TNFa','TNFa_pgml')]:
        wb = openpyxl.load_workbook('/home/user/elisa/ELISA.xlsx')
        ws = wb[sheet]
        rows = list(ws.iter_rows(values_only=True))
        for row in rows[1:]:
            donor, carp, stim, sc, val = row[:5]
            if donor is None: continue
            cn = int(re.sub(r'[^\d]','', str(carp)))
            sn = int(re.sub(r'[^\d]','', str(sc)))
            records.append({'Donor': donor, 'Carprofen': carp, 'carp_num': cn,
                            'Stimulus': stim, 'Stim_Conc': sc, 'stim_num': sn,
                            'cytokine': col,
                            'carp_x': 0.5 if cn == 0 else float(cn),
                            'value': float(val) if val is not None else np.nan})
    return pd.DataFrame(records)

elisa = load_elisa()

CYTOKINES = [
    ('IL1b_pgml', 'IL-1β (pg/mL)'),
    ('IL6_pgml',  'IL-6 (pg/mL)'),
    ('TNFa_pgml', 'TNF-α (pg/mL)'),
]

# ─────────────────────────────────────────────────────────────────────────────
# Log-log linear fit + 95% CI
# ─────────────────────────────────────────────────────────────────────────────
def log_log_fit(xs, ys, x_range):
    """Fit y ~ x in log10-log10 space. Returns (y_fit, y_lo, y_hi) or None."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    mask = (xs > 0) & (ys > 0) & np.isfinite(ys)
    if mask.sum() < 3:
        return None
    lx, ly = np.log10(xs[mask]), np.log10(ys[mask])
    n = len(lx)
    sl, ic, *_ = stats.linregress(lx, ly)
    ly_hat   = ic + sl * lx
    se_res   = np.sqrt(np.sum((ly - ly_hat)**2) / (n - 2))
    lx_pred  = np.log10(x_range)
    ly_pred  = ic + sl * lx_pred
    lx_mean  = lx.mean()
    ss_x     = np.sum((lx - lx_mean)**2)
    t95      = stats.t.ppf(0.975, n - 2)
    se_mean  = se_res * np.sqrt(1/n + (lx_pred - lx_mean)**2 / ss_x)
    return 10**ly_pred, 10**(ly_pred - t95*se_mean), 10**(ly_pred + t95*se_mean)

# ─────────────────────────────────────────────────────────────────────────────
# 2-way ANOVA
# ─────────────────────────────────────────────────────────────────────────────
def run_anova(df_sub):
    """df_sub already filtered to one cytokine & stimulus; has 'value' column."""
    d = df_sub.dropna(subset=['value'])
    d = d[d['value'] > 0].copy()
    if len(d) < 8:
        return ''
    def fp(p):
        if np.isnan(p): return 'n.s.'
        if p < 0.001: return 'p < 0.001'
        return f'p = {p:.3f}' + ('' if p < 0.05 else ' (n.s.)')
    try:
        m  = smf.ols('value ~ C(carp_num) * C(stim_num)', data=d).fit()
        tb = anova_lm(m, typ=2)
        ri = tb.index.tolist()
        return (f"2-way ANOVA\n"
                f"Carprofen: {fp(tb.loc[ri[0],'PR(>F)'])}\n"
                f"Stimulus:  {fp(tb.loc[ri[1],'PR(>F)'])}\n"
                f"Interaction: {fp(tb.loc[ri[2],'PR(>F)'])}")
    except:
        return ''

# ─────────────────────────────────────────────────────────────────────────────
# Draw one ELISA panel onto an existing Axes
# ─────────────────────────────────────────────────────────────────────────────
X_SMOOTH = np.logspace(np.log10(0.45), np.log10(125), 200)
X_BREAKS = [0.5, 1, 10, 100]
X_LABELS = ['0', '1', '10', '100']

def semilog_fit(xs, ys, x_range):
    """Fit y ~ log10(x) in log-x / linear-y space. Returns (y_fit, y_lo, y_hi)."""
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    mask = (xs > 0) & np.isfinite(ys) & (ys >= 0)
    if mask.sum() < 4:
        return None
    lx, ly = np.log10(xs[mask]), ys[mask]
    n = len(lx)
    sl, ic, *_ = stats.linregress(lx, ly)
    ly_hat  = ic + sl * lx
    se_res  = np.sqrt(np.sum((ly - ly_hat)**2) / (n - 2))
    lx_pred = np.log10(x_range)
    y_pred  = ic + sl * lx_pred
    lx_mean = lx.mean()
    ss_x    = np.sum((lx - lx_mean)**2)
    t95     = stats.t.ppf(0.975, n - 2)
    se_mean = se_res * np.sqrt(1/n + (lx_pred - lx_mean)**2 / ss_x)
    return y_pred, y_pred - t95*se_mean, y_pred + t95*se_mean

def draw_elisa_panel(ax, df_panel, cytokine_col, cyt_label, palette,
                     dose_levels, anova_txt='', log_y=True,
                     show_xlabel=True, show_ylabel=True, title=''):
    """Draw dose-response lines (log-x; log-y or linear-y depending on log_y)."""
    all_vals = []

    for dose in dose_levels:
        d   = df_panel[df_panel['stim_num'] == dose].copy()
        col = palette[dose]

        good = d.dropna(subset=['value'])
        good = good[good['value'] > (0 if log_y else -np.inf)]
        if log_y:
            good = good[good['value'] > 0]
        else:
            good = good[good['value'] >= 0]

        if len(good):
            ax.scatter(good['carp_x'], good['value'],
                       color=col, s=8, alpha=0.55, zorder=3, linewidths=0)
            all_vals.extend(good['value'].tolist())

        # Smooth + CI
        if log_y:
            fit = log_log_fit(good['carp_x'].values, good['value'].values,
                              X_SMOOTH) if len(good) >= 4 else None
        else:
            fit = semilog_fit(good['carp_x'].values, good['value'].values,
                              X_SMOOTH) if len(good) >= 4 else None

        if fit:
            yf, yl, yu = fit
            if not log_y:
                yl = np.maximum(yl, 0)   # CI can't go below 0 on linear scale
            ax.plot(X_SMOOTH, yf, color=col, linewidth=0.9, zorder=4,
                    label=str(dose))
            ax.fill_between(X_SMOOTH, yl, yu, color=col, alpha=0.12, zorder=2)
        elif len(good) >= 2:
            pts = good.sort_values('carp_x')
            ax.plot(pts['carp_x'], pts['value'], color=col, linewidth=0.7,
                    linestyle='--', zorder=4, label=str(dose), alpha=0.7)

    ax.set_xscale('log')
    ax.set_xticks(X_BREAKS)
    ax.set_xticklabels(X_LABELS)
    ax.set_xlim(0.4, 135)
    ax.xaxis.set_minor_formatter(NullFormatter())

    if log_y:
        ax.set_yscale('log')
        ax.yaxis.set_minor_formatter(NullFormatter())
        if all_vals:
            ax.set_ylim(max(0.05, np.nanmin(all_vals) * 0.3),
                        np.nanmax(all_vals) * 5)
    else:
        ax.set_yscale('linear')
        if all_vals:
            ax.set_ylim(0, np.nanmax(all_vals) * 1.35)

    if show_xlabel:
        ax.set_xlabel('Carprofen (µg/mL)')
    if show_ylabel:
        ax.set_ylabel(cyt_label)
    ax.set_title(title, pad=3)

    if anova_txt:
        ax.text(0.03, 0.97, anova_txt, transform=ax.transAxes,
                fontsize=5.5, color='#444444', va='top', ha='left',
                linespacing=1.4)

    ax.yaxis.grid(True, which='major', color='#E8E8E8', linewidth=0.3, zorder=0)
    ax.set_axisbelow(True)

# ─────────────────────────────────────────────────────────────────────────────
# Build full ELISA figure: 3 rows (cytokines) × 5 cols (D1-D4 + Average)
# ─────────────────────────────────────────────────────────────────────────────
def make_elisa_figure(stimulus, palette, dose_title, dose_levels,
                      fig_label, filename, log_y=True):
    """
    Produces TWO separate files:
      {filename}_donors.pdf   — 3 rows (cytokines) × 4 cols (donors)
      {filename}_average.pdf  — 1 row × 3 cols (cytokines), average + ANOVA
    """
    df = elisa[elisa['Stimulus'] == stimulus].copy()
    donors      = ['D1', 'D2', 'D3', 'D4']
    donor_titles= ['Donor 1', 'Donor 2', 'Donor 3', 'Donor 4']
    pan_let_d   = 'ABCDEFGHIJKL'   # 12 letters for 3×4
    pan_let_a   = 'ABC'             # 3 letters for average row

    # ── pre-compute ANOVA for each cytokine ──────────────────────────────────
    anova_txts = {}
    for col_id, _ in CYTOKINES:
        df_cyt = df[df['cytokine'] == col_id]
        anova_txts[col_id] = run_anova(df_cyt)

    # ─────────────────────────────────────────────────────────────────────────
    # Figure A: 3 rows × 4 cols  (individual donors)
    # ─────────────────────────────────────────────────────────────────────────
    fig_d, axes_d = plt.subplots(3, 4, figsize=(10.5, 7.5))
    fig_d.subplots_adjust(wspace=0.35, hspace=0.52)

    for ri, (col_id, cyt_label) in enumerate(CYTOKINES):
        df_cyt = df[df['cytokine'] == col_id]
        for ci, (donor, dtitle) in enumerate(zip(donors, donor_titles)):
            ax = axes_d[ri, ci]
            draw_elisa_panel(
                ax           = ax,
                df_panel     = df_cyt[df_cyt['Donor'] == donor],
                cytokine_col = col_id,
                cyt_label    = cyt_label,
                palette      = palette,
                dose_levels  = dose_levels,
                log_y        = log_y,
                show_xlabel  = (ri == 2),
                show_ylabel  = (ci == 0),
                title        = dtitle,
            )
            if ci == 0:
                ax.text(-0.20, 1.07, pan_let_d[ri],
                        transform=ax.transAxes,
                        fontsize=9, fontweight='bold', va='top')

    # Row labels on the right edge
    for ri, (_, cyt_label) in enumerate(CYTOKINES):
        axes_d[ri, -1].annotate(
            cyt_label, xy=(1.20, 0.5), xycoords='axes fraction',
            rotation=270, va='center', ha='left', fontsize=7, color='#333333')

    # Legend below figure
    handles = [plt.Line2D([0], [0], color=palette[d], linewidth=1.4,
                          label=str(d)) for d in dose_levels]
    fig_d.legend(handles=handles, title=dose_title,
                 title_fontsize=6.5, loc='lower center', fontsize=6.5,
                 handlelength=1.4, ncol=len(dose_levels),
                 bbox_to_anchor=(0.5, -0.02))

    fig_d.suptitle(fig_label + ' — individual donors',
                   fontsize=9, y=1.005, fontweight='normal')
    fig_d.savefig(filename + '_donors.pdf', bbox_inches='tight')
    fig_d.savefig(filename + '_donors.png', bbox_inches='tight')
    plt.close(fig_d)
    print(f'Saved {filename}_donors')

    # ─────────────────────────────────────────────────────────────────────────
    # Figure B: 1 row × 3 cols  (average + ANOVA)
    # ─────────────────────────────────────────────────────────────────────────
    fig_a, axes_a = plt.subplots(1, 3, figsize=(8.5, 3.0))
    fig_a.subplots_adjust(wspace=0.42)

    for ci, (col_id, cyt_label) in enumerate(CYTOKINES):
        ax = axes_a[ci]
        df_cyt = df[df['cytokine'] == col_id]
        draw_elisa_panel(
            ax           = ax,
            df_panel     = df_cyt,
            cytokine_col = col_id,
            cyt_label    = cyt_label,
            palette      = palette,
            dose_levels  = dose_levels,
            anova_txt    = anova_txts[col_id],
            log_y        = log_y,
            show_xlabel  = True,
            show_ylabel  = True,
            title        = cyt_label.split(' ')[0],   # just the cytokine name
        )
        ax.text(-0.20, 1.07, pan_let_a[ci],
                transform=ax.transAxes,
                fontsize=10, fontweight='bold', va='top')

    # Legend below figure
    handles = [plt.Line2D([0], [0], color=palette[d], linewidth=1.4,
                          label=str(d)) for d in dose_levels]
    fig_a.legend(handles=handles, title=dose_title,
                 title_fontsize=7, loc='lower center', fontsize=7,
                 handlelength=1.4, ncol=len(dose_levels),
                 bbox_to_anchor=(0.5, -0.06))

    fig_a.suptitle(fig_label + ' — average (n=4)',
                   fontsize=9, y=1.04, fontweight='normal')
    fig_a.savefig(filename + '_average.pdf', bbox_inches='tight')
    fig_a.savefig(filename + '_average.png', bbox_inches='tight')
    plt.close(fig_a)
    print(f'Saved {filename}_average')

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — LPS
# ─────────────────────────────────────────────────────────────────────────────
make_elisa_figure(
    stimulus    = 'LPS',
    palette     = LPS_PAL,
    dose_title  = 'LPS (ng/mL)',
    dose_levels = [0, 1, 10, 100],
    log_y       = True,
    fig_label   = 'Figure 2 | ELISA: Effect of carprofen on LPS-induced cytokine production in human MDMs',
    filename    = '/home/user/elisa/Figure2_ELISA_LPS',
)

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — BCG  (log x, linear y — data too sparse for log-y)
# ─────────────────────────────────────────────────────────────────────────────
make_elisa_figure(
    stimulus    = 'BCG',
    palette     = BCG_PAL,
    dose_title  = 'BCG (MOI)',
    dose_levels = [0, 1, 10, 100],
    log_y       = False,
    fig_label   = 'Figure 3 | ELISA: Effect of carprofen on BCG-induced cytokine production in human MDMs',
    filename    = '/home/user/elisa/Figure3_ELISA_BCG',
)

# ─────────────────────────────────────────────────────────────────────────────
# BCG flow cytometry — load & compute fold change
# ─────────────────────────────────────────────────────────────────────────────
def load_bcg():
    records = []
    file_map = {1:('Media only','D1'), 2:('Media only','D2'), 3:('Media only','D3'),
                4:('+Carprofen','D1'), 5:('+Carprofen','D2'), 6:('+Carprofen','D3')}
    for fnum,(cond,donor) in file_map.items():
        wb  = xlrd.open_workbook(f'/home/user/elisa/AIF007-{fnum}_Table.xls')
        sh  = wb.sheet_by_name('Sheet0')
        for r in range(1, sh.nrows):
            row = sh.row_values(r)
            s   = str(row[0]).replace('.fcs','')
            if not s or 'Mean' in s or 'SD' in s: continue
            m = re.match(r'(\d+)[a-z]+_(\w+)_C(\d+)B(\d+)', s)
            if not m: continue
            tp, frac, carp, moi = int(m[1]),m[2],int(m[3]),int(m[4])
            beads, bcg = float(row[2] or 0), float(row[6] or 0)
            if beads == 0: continue
            records.append(dict(condition=cond, donor=donor, timepoint=tp,
                                fraction=frac, carprofen=carp, moi=moi,
                                norm=bcg/beads))
    df  = pd.DataFrame(records)
    KEY = ['condition','donor','fraction','carprofen','moi']
    t4  = df[df['timepoint']==4][KEY+['norm']].rename(columns={'norm':'n4'})
    t7  = df[df['timepoint']==7][KEY+['norm']].rename(columns={'norm':'n7'})
    fc  = t4.merge(t7, on=KEY)
    fc['fold_change'] = fc['n7'] / fc['n4']
    fc['carp_x'] = fc['carprofen'].apply(lambda v: 0.5 if v==0 else float(v))
    return fc.dropna(subset=['fold_change'])

fc_df = load_bcg()

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 — BCG growth assay: bar chart + SEM + donor dots
#   2 rows (EC, IC)  ×  4 columns (MOI 0, 1, 10, 100)
#   x = carprofen (categorical), bars = Media only vs +Carprofen
# ─────────────────────────────────────────────────────────────────────────────
COND_ORDER  = ['Media only', '+Carprofen']
FRAC_ORDER  = ['EC', 'IC']
MOI_LEVELS  = [0, 1, 10, 100]
MOI_LABELS  = ['MOI 0', 'MOI 1', 'MOI 10', 'MOI 100']
CARP_CATS   = [0, 1, 10, 100]
CARP_XLABS  = ['0', '1', '10', '100']
COND_COLORS = {'Media only': '#5B9BD5', '+Carprofen': '#E06060'}
BAR_W       = 0.35
PANEL_LET4  = [['A','B','C','D'],['E','F','G','H']]

fig4, axes4 = plt.subplots(2, 4, figsize=(10, 5.5))
fig4.subplots_adjust(wspace=0.35, hspace=0.55)

rng = np.random.default_rng(42)

for ri, frac in enumerate(FRAC_ORDER):
    for ci, moi in enumerate(MOI_LEVELS):
        ax = axes4[ri, ci]

        x_pos = np.arange(len(CARP_CATS))

        for bi, cond in enumerate(COND_ORDER):
            offsets = x_pos + (bi - 0.5) * BAR_W
            col     = COND_COLORS[cond]

            means, sems, donor_vals = [], [], []
            for carp in CARP_CATS:
                vals = fc_df[
                    (fc_df['fraction']  == frac) &
                    (fc_df['condition'] == cond) &
                    (fc_df['moi']       == moi)  &
                    (fc_df['carprofen'] == carp)
                ]['fold_change'].dropna().values
                m = np.mean(vals) if len(vals) else np.nan
                s = (np.std(vals, ddof=1)/np.sqrt(len(vals))
                     if len(vals) > 1 else 0)
                means.append(m); sems.append(s); donor_vals.append(vals)

            # Bars
            ax.bar(offsets, means, width=BAR_W*0.88, color=col,
                   zorder=2, label=cond)
            # Error bars (SEM)
            ax.errorbar(offsets, means, yerr=sems,
                        fmt='none', ecolor='#222222',
                        elinewidth=0.8, capsize=2.5, capthick=0.8, zorder=3)
            # Individual donor dots
            for xi, vals in zip(offsets, donor_vals):
                if len(vals):
                    jit = rng.uniform(-BAR_W*0.18, BAR_W*0.18, len(vals))
                    ax.scatter(xi + jit, vals, s=10, color='black',
                               zorder=4, linewidths=0, alpha=0.75)

        # Reference y = 1
        ax.axhline(1, color='#888888', linewidth=0.6,
                   linestyle='--', zorder=1)

        ax.set_xticks(x_pos)
        ax.set_xticklabels(CARP_XLABS, fontsize=6.5)
        ax.set_ylim(bottom=0)
        ax.yaxis.grid(True, which='major', color='#E8E8E8',
                      linewidth=0.3, zorder=0)
        ax.set_axisbelow(True)

        # Column title = MOI (top row only)
        if ri == 0:
            ax.set_title(MOI_LABELS[ci], pad=3, fontsize=7.5)
        # Row label = fraction (left column only)
        if ci == 0:
            ax.set_ylabel(f'{frac}\nFold change (norm. 4 hpi)', fontsize=7)
        # x label (bottom row only)
        if ri == 1:
            ax.set_xlabel('Carprofen (µg/mL)', fontsize=7)

        # Panel letter
        ax.text(-0.20 if ci==0 else -0.12, 1.07,
                PANEL_LET4[ri][ci], transform=ax.transAxes,
                fontsize=9, fontweight='bold', va='top')

# Shared legend
import matplotlib.patches as mpatches
handles4 = [mpatches.Patch(color=COND_COLORS[c], label=c)
            for c in COND_ORDER]
fig4.legend(handles=handles4, title_fontsize=7,
            loc='lower center', ncol=2, bbox_to_anchor=(0.5, -0.04),
            fontsize=7, handlelength=1.2)

fig4.suptitle('Figure 4 | BCG growth assay: fold change in extracellular and intracellular BCG',
              fontsize=9, y=1.01)
fig4.savefig('/home/user/elisa/Figure4_BCG_growth_assay.pdf')
fig4.savefig('/home/user/elisa/Figure4_BCG_growth_assay.png')
plt.close(fig4)
print('Saved Figure4_BCG_growth_assay')

print('\nAll figures done.')
