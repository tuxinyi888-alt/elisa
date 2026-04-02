"""Generate ELISA cytokine plots using matplotlib (Python equivalent of the R script)."""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os

matplotlib.use('Agg')

os.makedirs('plots', exist_ok=True)

# Sheet definitions
sheets = [
    {'name': 'IL1b', 'col': 'IL1b_pgml', 'label': 'IL-1\u03b2 (pg/ml)'},
    {'name': 'IL6',  'col': 'IL6_pgml',  'label': 'IL-6 (pg/ml)'},
    {'name': 'TNFa', 'col': 'TNFa_pgml', 'label': 'TNF-\u03b1 (pg/ml)'},
]

carp_order = ['C0', 'C1', 'C10', 'C100']
carp_labels = ['0', '1', '10', '100']

lps_levels = ['L0', 'L1', 'L10', 'L100']
lps_labels = ['LPS 0 ng/ml', 'LPS 1 ng/ml', 'LPS 10 ng/ml', 'LPS 100 ng/ml']
bcg_levels = ['B0', 'B1', 'B10', 'B100']
bcg_labels = ['BCG 0 MOI', 'BCG 1 MOI', 'BCG 10 MOI', 'BCG 100 MOI']

colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a']


def make_plot(data, value_col, y_label, stim_type, title, filename, show_errorbars=False, sd_col=None):
    fig, ax = plt.subplots(figsize=(6, 4))

    if stim_type == 'LPS':
        levels, labels = lps_levels, lps_labels
    else:
        levels, labels = bcg_levels, bcg_labels

    for i, (lev, lab) in enumerate(zip(levels, labels)):
        subset = data[data['Stim_Conc'] == lev].copy()
        subset['Carprofen'] = pd.Categorical(subset['Carprofen'], categories=carp_order, ordered=True)
        subset = subset.sort_values('Carprofen')

        x = [carp_order.index(c) for c in subset['Carprofen']]
        y = subset[value_col].fillna(0).values

        ax.plot(x, y, marker='o', color=colors[i], label=lab, linewidth=1.2, markersize=5)

        if show_errorbars and sd_col is not None:
            sd = subset[sd_col].fillna(0).values
            ax.errorbar(x, y, yerr=sd, fmt='none', ecolor=colors[i], capsize=3, linewidth=0.8)

    ax.set_xticks(range(len(carp_labels)))
    ax.set_xticklabels(carp_labels)
    ax.set_xlabel('Carprofen (\u03bcM)')
    ax.set_ylabel(y_label)
    ax.set_title(title, fontweight='bold', fontsize=11)
    ax.legend(fontsize=8, frameon=True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f'Saved: {filename}')


donors = ['D1', 'D2', 'D3', 'D4']
stimuli = ['LPS', 'BCG']

for s in sheets:
    df = pd.read_excel('ELISA.xlsx', sheet_name=s['name'], usecols=[0, 1, 2, 3, 4])
    df.columns = ['Donor', 'Carprofen', 'Stimulus', 'Stim_Conc', s['col']]

    for stim in stimuli:
        stim_data = df[df['Stimulus'] == stim]

        # Per-donor plots
        for d in donors:
            donor_data = stim_data[stim_data['Donor'] == d]
            title = f"{s['label']} - {d} ({stim} + Carprofen)"
            fname = f"plots/{s['name']}_{stim}_{d}.png"
            make_plot(donor_data, s['col'], s['label'], stim, title, fname)

        # Average plot
        avg = stim_data.groupby(['Carprofen', 'Stimulus', 'Stim_Conc'])[s['col']].agg(['mean', 'std']).reset_index()
        avg.columns = ['Carprofen', 'Stimulus', 'Stim_Conc', s['col'], 'sd']
        avg['sd'] = avg['sd'].fillna(0)

        title = f"{s['label']} - Average of 4 Donors ({stim} + Carprofen)"
        fname = f"plots/{s['name']}_{stim}_Average.png"
        make_plot(avg, s['col'], s['label'], stim, title, fname, show_errorbars=True, sd_col='sd')

print('\nAll plots saved in the plots/ directory!')
