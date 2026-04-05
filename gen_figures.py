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

def draw_elisa_panel(ax, df_panel, cytokine_col, cyt_label, palette,
                     dose_levels, anova_txt='',
                     show_xlabel=True, show_ylabel=True, title=''):
    """Draw dose-response lines for one panel (one donor or average)."""
    all_vals = []  # track actual values to set y-limits from data, not fit

    for dose in dose_levels:
        d   = df_panel[df_panel['stim_num'] == dose].copy()
        col = palette[dose]

        good = d.dropna(subset=['value'])
        good = good[good['value'] > 0]
        if len(good):
            ax.scatter(good['carp_x'], good['value'],
                       color=col, s=8, alpha=0.55, zorder=3, linewidths=0)
            all_vals.extend(good['value'].tolist())

        # Only draw smooth if ≥4 data points (avoids wild extrapolation)
        fit = log_log_fit(good['carp_x'].values, good['value'].values, X_SMOOTH) \
              if len(good) >= 4 else None
        if fit:
            yf, yl, yu = fit
            ax.plot(X_SMOOTH, yf, color=col, linewidth=0.9, zorder=4,
                    label=str(dose))
            ax.fill_between(X_SMOOTH, yl, yu, color=col, alpha=0.12, zorder=2)
        elif len(good) >= 2:
            # Too few points for CI — just connect them with a thin line
            pts = good.sort_values('carp_x')
            ax.plot(pts['carp_x'], pts['value'], color=col, linewidth=0.7,
                    linestyle='--', zorder=4, label=str(dose), alpha=0.7)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks(X_BREAKS)
    ax.set_xticklabels(X_LABELS)
    ax.set_xlim(0.4, 135)
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())

    # Y-axis limits: based only on actual data (not the fitted curve)
    if all_vals:
        ylo = max(0.05, np.nanmin(all_vals) * 0.3)
        yhi = np.nanmax(all_vals) * 5
        ax.set_ylim(ylo, yhi)

    if show_xlabel:
        ax.set_xlabel('Carprofen (µg/mL)')
    if show_ylabel:
        ax.set_ylabel(cyt_label)
    ax.set_title(title, pad=3)

    if anova_txt:
        ax.text(0.03, 0.02, anova_txt, transform=ax.transAxes,
                fontsize=5.5, color='#444444', va='bottom', ha='left',
                linespacing=1.4)

    ax.yaxis.grid(True, which='major', color='#E8E8E8', linewidth=0.3, zorder=0)
    ax.set_axisbelow(True)

# ─────────────────────────────────────────────────────────────────────────────
# Build full ELISA figure: 3 rows (cytokines) × 5 cols (D1-D4 + Average)
# ─────────────────────────────────────────────────────────────────────────────
def make_elisa_figure(stimulus, palette, dose_title, dose_levels,
                      fig_label, filename):
    df = elisa[elisa['Stimulus'] == stimulus].copy()
    donors     = ['D1','D2','D3','D4']
    col_titles = ['Donor 1','Donor 2','Donor 3','Donor 4','Average (n=4)']
    panel_letters = 'ABCDEFGHIJKLMNO'

    fig, axes = plt.subplots(3, 5, figsize=(13, 7.5))
    fig.subplots_adjust(wspace=0.38, hspace=0.52)

    for ri, (col_id, cyt_label) in enumerate(CYTOKINES):
        # Pre-compute ANOVA on full dataset for the average panel
        df_cyt = df[df['cytokine'] == col_id].copy()
        anova_txt = run_anova(df_cyt)

        for ci, (donor, col_title) in enumerate(zip(donors + [None], col_titles)):
            ax    = axes[ri, ci]
            is_avg = donor is None
            df_p  = df if is_avg else df[df['Donor'] == donor]

            df_p2 = (df_cyt if is_avg else df_cyt[df_cyt['Donor'] == donor])
            draw_elisa_panel(
                ax          = ax,
                df_panel    = df_p2,
                cytokine_col= col_id,
                cyt_label   = cyt_label,
                palette     = palette,
                dose_levels = dose_levels,
                anova_txt   = anova_txt if is_avg else '',
                show_xlabel = (ri == 2),
                show_ylabel = (ci == 0),
                title       = col_title,
            )

            # Panel letter
            letter = panel_letters[ri*5 + ci]
            ax.text(-0.18 if ci==0 else -0.12, 1.06, letter,
                    transform=ax.transAxes,
                    fontsize=9, fontweight='bold', va='top')

    # Row labels (cytokine names) on the right
    for ri, (_, cyt_label) in enumerate(CYTOKINES):
        axes[ri, -1].annotate(
            cyt_label, xy=(1.18, 0.5), xycoords='axes fraction',
            rotation=270, va='center', ha='left', fontsize=7, color='#333333')

    # Legend (top-right panel)
    leg_ax = axes[0, 4]
    handles = [plt.Line2D([0],[0], color=palette[d], linewidth=1.4,
                          label=str(d)) for d in dose_levels]
    leg_ax.legend(handles=handles, title=dose_title, title_fontsize=6.5,
                  loc='lower right', fontsize=6.5, handlelength=1.5)

    fig.suptitle(fig_label, fontsize=9, y=1.005, fontweight='normal')
    fig.savefig(filename + '.pdf')
    fig.savefig(filename + '.png')
    plt.close(fig)
    print(f'Saved {filename}')

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — LPS
# ─────────────────────────────────────────────────────────────────────────────
make_elisa_figure(
    stimulus    = 'LPS',
    palette     = LPS_PAL,
    dose_title  = 'LPS (ng/mL)',
    dose_levels = [0, 1, 10, 100],
    fig_label   = 'Figure 2 | ELISA: Effect of carprofen on LPS-induced cytokine production in human MDMs',
    filename    = '/home/user/elisa/Figure2_ELISA_LPS',
)

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — BCG
# ─────────────────────────────────────────────────────────────────────────────
make_elisa_figure(
    stimulus    = 'BCG',
    palette     = BCG_PAL,
    dose_title  = 'BCG (MOI)',
    dose_levels = [0, 1, 10, 100],
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
# Figure 4 — BCG growth assay: 2 rows (EC, IC) × 2 cols (conditions)
# ─────────────────────────────────────────────────────────────────────────────
COND_ORDER = ['Media only', '+Carprofen']
FRAC_ORDER = ['EC', 'IC']
MOI_LEVELS = [0, 1, 10, 100]
PANEL_LET  = ['A','B','C','D']

fig4, axes4 = plt.subplots(2, 2, figsize=(7.5, 6.0))
fig4.subplots_adjust(wspace=0.38, hspace=0.48)

for ri, frac in enumerate(FRAC_ORDER):
    for ci, cond in enumerate(COND_ORDER):
        ax  = axes4[ri, ci]
        d   = fc_df[(fc_df['fraction']==frac) & (fc_df['condition']==cond)]

        for moi in MOI_LEVELS:
            dm  = d[d['moi']==moi]
            col = BCG_PAL[moi]

            # Scatter
            if len(dm):
                ax.scatter(dm['carp_x'], dm['fold_change'],
                           color=col, s=10, alpha=0.55, zorder=3, linewidths=0)

            # Smooth + CI (fold change: fit on log-x, linear y)
            xs, ys = dm['carp_x'].values, dm['fold_change'].values
            mask   = (xs > 0) & np.isfinite(ys)
            if mask.sum() >= 3:
                lx    = np.log10(xs[mask])
                ly    = ys[mask]
                n     = len(lx)
                sl,ic,*_ = stats.linregress(lx, ly)
                ly_hat = ic + sl*lx
                se_res = np.sqrt(np.sum((ly-ly_hat)**2)/(n-2))
                lx_pred= np.log10(X_SMOOTH)
                y_pred = ic + sl*lx_pred
                lx_m   = lx.mean()
                ss_x   = np.sum((lx-lx_m)**2)
                t95    = stats.t.ppf(0.975, n-2)
                if ss_x > 0:
                    se_m = se_res*np.sqrt(1/n + (lx_pred-lx_m)**2/ss_x)
                else:
                    se_m = np.zeros_like(lx_pred)
                ax.plot(X_SMOOTH, y_pred, color=col, linewidth=0.9, zorder=4,
                        label=str(moi))
                ax.fill_between(X_SMOOTH, y_pred-t95*se_m, y_pred+t95*se_m,
                                color=col, alpha=0.12, zorder=2)

        # Reference y=1
        ax.axhline(1, color='#888888', linewidth=0.6, linestyle='--', zorder=1)

        ax.set_xscale('log')
        ax.set_xticks(X_BREAKS)
        ax.set_xticklabels(X_LABELS)
        ax.set_xlim(0.4, 135)
        ax.xaxis.set_minor_formatter(NullFormatter())

        # Y limits: from data, never negative
        panel_d = fc_df[(fc_df['fraction']==frac) & (fc_df['condition']==cond)]
        fc_vals = panel_d['fold_change'].dropna().values
        if len(fc_vals):
            ylo = max(0, np.percentile(fc_vals, 2) * 0.7)
            yhi = np.percentile(fc_vals, 98) * 1.5
            ax.set_ylim(ylo, max(yhi, 1.5))

        ax.set_title(cond, pad=3)
        if ri == 1:
            ax.set_xlabel('Carprofen (µg/mL)')
        if ci == 0:
            ax.set_ylabel(f'{frac} fold change\n(normalised to 4 hpi)')

        ax.yaxis.grid(True, which='major', color='#E8E8E8', linewidth=0.3, zorder=0)
        ax.set_axisbelow(True)

        # Panel letter
        ax.text(-0.15, 1.06, PANEL_LET[ri*2+ci], transform=ax.transAxes,
                fontsize=10, fontweight='bold', va='top')

# Shared legend
handles4 = [plt.Line2D([0],[0], color=BCG_PAL[m], linewidth=1.4,
                        label=str(m)) for m in MOI_LEVELS]
fig4.legend(handles=handles4, title='BCG (MOI)', title_fontsize=7,
            loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.04),
            fontsize=7, handlelength=1.5)

fig4.suptitle('Figure 4 | BCG growth assay: fold change in extracellular and intracellular BCG',
              fontsize=9, y=1.01)
fig4.savefig('/home/user/elisa/Figure4_BCG_growth_assay.pdf')
fig4.savefig('/home/user/elisa/Figure4_BCG_growth_assay.png')
plt.close(fig4)
print('Saved Figure4_BCG_growth_assay')

print('\nAll figures done.')
