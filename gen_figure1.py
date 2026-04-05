"""
Figure 1 — Mechanism of carprofen as host-directed therapy in human MDMs.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse
import matplotlib.patheffects as pe

# ── palette ───────────────────────────────────────────────────────────────────
CARP_C   = '#2E86AB'
LPS_C    = '#C9392B'
BCG_C    = '#2D6A4F'
PHAGO_C  = '#B7E4C7'
IL1_C    = '#E76F51'
NEUTRAL_C= '#7A8FA6'
SIGNAL_C = '#457B9D'
COX_C    = '#9B72CF'
PGE_C    = '#C8A2C8'
CELL_BG  = '#FDFBF5'
CELL_EDGE= '#A08C6E'
DARK     = '#333333'

# ── figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(15, 9.5))
ax.set_xlim(0, 15)
ax.set_ylim(0, 9.5)
ax.axis('off')
fig.patch.set_facecolor('white')

# ── helpers ───────────────────────────────────────────────────────────────────
def box(cx, cy, w, h, fc, text='', fs=9, tc='white',
        ec=None, lw=1.8, z=6, ls='-', alpha=1.0):
    ec2 = ec or fc
    ax.add_patch(FancyBboxPatch(
        (cx - w/2, cy - h/2), w, h,
        boxstyle='round,pad=0.1', facecolor=fc, edgecolor=ec2,
        linewidth=lw, linestyle=ls, zorder=z, alpha=alpha))
    if text:
        ax.text(cx, cy, text, ha='center', va='center', fontsize=fs,
                color=tc, fontweight='bold', zorder=z+1, linespacing=1.35)

def arr(x1, y1, x2, y2, col=DARK, lw=1.6, rad=0.0, z=5, ms=13,
        style='->'):
    cs = f'arc3,rad={rad}'
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, color=col, lw=lw,
                                mutation_scale=ms, connectionstyle=cs),
                zorder=z)

def inhibit(x1, y1, x2, y2, col='#C9392B', lw=2.0, z=5):
    """Flat-head inhibition symbol ⊣."""
    arr(x1, y1, x2, y2, col=col, lw=lw, z=z, style='-')
    dx, dy = x2-x1, y2-y1
    L = np.hypot(dx, dy) or 1
    px, py = -dy/L * 0.18, dx/L * 0.18
    ax.plot([x2-px, x2+px], [y2-py, y2+py], color=col, lw=2.8,
            solid_capstyle='round', zorder=z)

def txt(x, y, s, fs=8, col=DARK, ha='center', va='center',
        bold=False, italic=False, z=8):
    ax.text(x, y, s, ha=ha, va=va, fontsize=fs, color=col,
            fontweight='bold' if bold else 'normal',
            style='italic' if italic else 'normal',
            zorder=z, linespacing=1.35)

def rod(cx, cy, ang=20, sc=1.0, z=9):
    ax.add_patch(Ellipse((cx, cy), 0.44*sc, 0.17*sc, angle=ang,
                          facecolor=BCG_C, edgecolor='#1B4332',
                          linewidth=0.8, zorder=z))

# ═════════════════════════════════════════════════════════════════════════════
# MACROPHAGE
# ═════════════════════════════════════════════════════════════════════════════
ax.add_patch(Ellipse((7.5, 4.3), width=13.8, height=7.6,
                      facecolor=CELL_BG, edgecolor=CELL_EDGE,
                      linewidth=3, zorder=1))
txt(7.5, 0.55, 'Human monocyte-derived macrophage (MDM)',
    fs=9, col=CELL_EDGE, bold=True)
# nucleus hint
ax.add_patch(Ellipse((7.5, 4.0), 2.4, 1.6,
                      facecolor='#EDE8D8', edgecolor=CELL_EDGE,
                      linewidth=1.2, linestyle='--', zorder=2, alpha=0.5))

# ═════════════════════════════════════════════════════════════════════════════
# CARPROFEN  (top centre)
# ═════════════════════════════════════════════════════════════════════════════
box(7.5, 8.75, 2.8, 0.65, CARP_C, 'Carprofen', fs=12, lw=2.5, z=10)
txt(7.5, 8.20, '(microbicidal NSAID)', fs=7.5, col=CARP_C, italic=True)

# ═════════════════════════════════════════════════════════════════════════════
# LEFT BRANCH — LPS / immunomodulation
# Column A: x=1.9  Column B (COX): x=4.6
# ═════════════════════════════════════════════════════════════════════════════

# LPS (extracellular)
box(1.9, 8.55, 1.5, 0.55, LPS_C, 'LPS', fs=10.5, lw=2, z=10)
txt(1.9, 8.08, '(extracellular)', fs=7, col=LPS_C, italic=True)

# TLR4
box(1.9, 7.05, 1.65, 0.52, '#7F5539', 'TLR4', fs=9.5, z=8)
arr(1.9, 8.27, 1.9, 7.33, col=LPS_C, lw=2.0)

# NF-κB
box(1.9, 5.80, 1.72, 0.54, SIGNAL_C, 'NF-κB', fs=10, z=8)
arr(1.9, 6.79, 1.9, 6.08, col=SIGNAL_C)

# IL-6, TNF-α  (right-branch of NF-κB)
box(4.0, 5.80, 2.0, 0.58, NEUTRAL_C, 'IL-6, TNF-α', fs=9, z=8)
arr(2.77, 5.80, 3.0, 5.80, col=NEUTRAL_C, lw=1.6)
txt(4.0, 5.24, 'not significantly\naffected by carprofen',
    fs=7, col=NEUTRAL_C)

# NF-κB → NLRP3 priming
box(1.9, 4.50, 2.1, 0.66, '#E9C46A',
    'NLRP3\ninflammasome', fs=8.5, tc='#333', z=8)
arr(1.9, 5.53, 1.9, 4.84, col=SIGNAL_C)

# NLRP3 → IL-1β
box(1.9, 3.22, 1.75, 0.55, IL1_C, 'IL-1β ↑', fs=11, z=8)
arr(1.9, 4.17, 1.9, 3.51, col=IL1_C, lw=2.2)
txt(2.68, 3.83, 'caspase-1', fs=7, col=IL1_C, italic=True)

# synergy callout
txt(0.72, 3.22, '~4.6×\nsynergy\nwith LPS', fs=7, col=IL1_C)
ax.annotate('', xy=(1.03, 3.22), xytext=(1.48, 3.22),
            arrowprops=dict(arrowstyle='<-', color=IL1_C, lw=1.3,
                            mutation_scale=11), zorder=7)

# —— COX arm ——
# Carprofen → COX (inhibition)
box(4.6, 7.25, 1.85, 0.54, COX_C, 'COX-1/2', fs=9, z=8)
arr(6.13, 8.58, 5.05, 7.52, col=CARP_C, lw=2.0)
inhibit(4.6, 6.98, 4.6, 6.60, col=CARP_C, lw=2.0)

# ↓ PGE₂
box(4.6, 6.24, 1.85, 0.54, PGE_C, '↓ PGE₂', fs=9.5, tc='#333', z=8)

# ↓PGE₂ relieves NLRP3 suppression (curved arrow)
arr(3.67, 6.18, 2.95, 4.72, col=PGE_C, lw=1.7, rad=-0.30)
txt(3.72, 5.38, 'relieves\nEP2/EP4\nsuppression', fs=7, col=PGE_C)

# ═════════════════════════════════════════════════════════════════════════════
# RIGHT BRANCH — BCG / direct antimycobacterial
# ═════════════════════════════════════════════════════════════════════════════

# BCG extracellular
for (bx, by, ang) in [(11.55, 8.62, 12), (12.25, 8.70, -8), (11.9, 8.22, 28)]:
    rod(bx, by, ang=ang)
txt(11.9, 7.90, 'BCG (extracellular)', fs=8.5, col=BCG_C, bold=True)

# Carprofen → ETC disruption → EC BCG (direct arrow)
arr(8.87, 8.75, 11.1, 8.65, col=CARP_C, lw=2.0)
txt(9.98, 9.0, 'ETC disruption', fs=8, col=CARP_C)

# ↓ EC BCG growth box
box(11.9, 7.18, 2.45, 0.58, BCG_C, '↓ EC BCG growth', fs=8.5, z=8)
txt(11.9, 6.73, '(transient or sustained\nexposure)', fs=7, col=BCG_C,
    italic=True)
arr(11.9, 7.90, 11.9, 7.49, col=BCG_C, lw=1.6)

# —— Phagosome ——
ax.add_patch(Ellipse((11.1, 4.45), 3.4, 2.7,
                      facecolor=PHAGO_C, edgecolor=BCG_C,
                      linewidth=2.0, zorder=3, alpha=0.7))
txt(11.1, 5.42, 'Phagosome', fs=8.5, col='#1B4332', italic=True)
for (bx, by, ang) in [(10.75, 4.50, 10), (11.35, 4.42, -18),
                       (11.05, 3.98, 22), (10.6, 3.95, 5)]:
    rod(bx, by, ang=ang, z=9)
txt(11.1, 3.55, 'BCG (intracellular)', fs=8.5, col=BCG_C, bold=True, z=9)

# Carprofen → phagosome (sustained accumulation)
arr(8.25, 8.35, 9.92, 5.82, col=CARP_C, lw=1.8, rad=0.22)
txt(8.60, 7.0, 'drug accumulates\nin phagosome\n(sustained)', fs=7.5,
    col=CARP_C, italic=True)

# Phagosome → outcome boxes
arr(10.30, 3.1, 9.55, 2.82, col=DARK, lw=1.5)
arr(11.70, 3.1, 12.55, 2.82, col=DARK, lw=1.5)

# Outcome: transient
box(9.05, 2.44, 2.7, 0.90, '#FFF3CD',
    'Transient (4 h)\n↑ IC BCG', fs=8.5,
    tc='#7D4E00', ec='#E9C46A', lw=2.0, z=8)

# Outcome: sustained
box(12.95, 2.44, 2.7, 0.90, '#D8F3DC',
    'Sustained (7 d)\n↓ IC BCG', fs=8.5,
    tc='#1B4332', ec=BCG_C, lw=2.0, z=8)

# ═════════════════════════════════════════════════════════════════════════════
# LEGEND
# ═════════════════════════════════════════════════════════════════════════════
lx, ly = 0.3, 2.2
txt(lx, ly + 0.65, 'Key:', fs=8.5, col=DARK, bold=True, ha='left')
ax.annotate('', xy=(lx + 0.6, ly + 0.18), xytext=(lx + 0.05, ly + 0.18),
            arrowprops=dict(arrowstyle='->', color=DARK, lw=1.5,
                            mutation_scale=12), zorder=8)
txt(lx + 1.35, ly + 0.18, 'Activation / induction',
    fs=8, col=DARK, ha='left')

inhibit(lx + 0.05, ly - 0.32, lx + 0.60, ly - 0.32, col='#C9392B', lw=2.0)
txt(lx + 1.35, ly - 0.32, 'Inhibition', fs=8, col='#C9392B', ha='left')

# ═════════════════════════════════════════════════════════════════════════════
# TITLE
# ═════════════════════════════════════════════════════════════════════════════
ax.text(7.5, 9.35,
        'Figure 1 | Proposed mechanism of carprofen as host-directed therapy in human MDMs',
        ha='center', va='center', fontsize=10.5, color='#333333')

for ext in ('pdf', 'png'):
    fig.savefig(f'/home/user/elisa/Figure1_mechanism.{ext}',
                bbox_inches='tight', dpi=300)
plt.close(fig)
print('Saved Figure1_mechanism')
