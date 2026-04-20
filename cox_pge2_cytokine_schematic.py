"""Schematic: LPS -> TLR4 -> MyD88 -> NF-kB -> COX-2/PGE2 & cytokines.

Renders a publication-style diagram showing the hypothesised modulatory
link between PGE2 and pro-inflammatory cytokine responses
(TNF-alpha, IL-1beta, IL-6), and the point of action of COX inhibitors.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D


def box(ax, xy, w, h, text, facecolor, edgecolor="#333333",
        textcolor="#111111", fontsize=10, fontweight="bold"):
    x, y = xy
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        linewidth=1.4, facecolor=facecolor, edgecolor=edgecolor,
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center",
            fontsize=fontsize, fontweight=fontweight, color=textcolor)
    return (x, y, w, h)


def arrow(ax, start, end, color="#222222", lw=1.8, style="-|>",
          connectionstyle="arc3,rad=0.0", linestyle="-"):
    a = FancyArrowPatch(
        start, end,
        arrowstyle=style, mutation_scale=16,
        linewidth=lw, color=color,
        connectionstyle=connectionstyle, linestyle=linestyle,
        shrinkA=6, shrinkB=6,
    )
    ax.add_patch(a)


def inhibit_arrow(ax, start, end, color="#b30000", lw=2.0):
    a = FancyArrowPatch(
        start, end,
        arrowstyle="-[", mutation_scale=14,
        linewidth=lw, color=color, shrinkA=4, shrinkB=4,
    )
    ax.add_patch(a)


fig, ax = plt.subplots(figsize=(12, 9.5))
ax.set_xlim(0, 12)
ax.set_ylim(0, 9)
ax.set_aspect("equal")
ax.axis("off")

# ---- Cellular compartments ---------------------------------------------
# Extracellular band
ax.add_patch(Rectangle((0.2, 7.7), 11.6, 1.1,
                       facecolor="#eaf4ff", edgecolor="none", zorder=0))
ax.text(0.35, 8.55, "Extracellular", fontsize=9, style="italic",
        color="#446")

# Plasma membrane (double line)
for y in (7.55, 7.45):
    ax.add_patch(Rectangle((0.2, y), 11.6, 0.03,
                           facecolor="#f6c26b", edgecolor="none", zorder=0))
ax.text(0.35, 7.30, "Plasma membrane", fontsize=9, style="italic",
        color="#7a5a20")

# Cytoplasm
ax.add_patch(Rectangle((0.2, 2.2), 11.6, 5.1,
                       facecolor="#fbfbf3", edgecolor="none", zorder=0))
ax.text(0.35, 7.10, "Cytoplasm", fontsize=9, style="italic", color="#555")

# Nucleus
ax.add_patch(FancyBboxPatch((3.2, 4.4), 5.6, 1.3,
                            boxstyle="round,pad=0.02,rounding_size=0.25",
                            facecolor="#f3ecfb", edgecolor="#8a6fb3",
                            linewidth=1.2, zorder=0))
ax.text(3.35, 5.55, "Nucleus", fontsize=9, style="italic", color="#553b83")

# ---- Nodes -------------------------------------------------------------
# LPS (extracellular)
box(ax, (6.0, 8.25), 1.4, 0.55, "LPS",
    facecolor="#d7ebff", fontsize=11)

# TLR4 / TLRs (membrane)
box(ax, (6.0, 7.05), 2.2, 0.55, "TLR4 / TLRs",
    facecolor="#ffe6bf", fontsize=11)

# MyD88 (adapter)
box(ax, (6.0, 6.15), 1.8, 0.55, "MyD88", facecolor="#ffd9d9", fontsize=11)

# NF-kB (in nucleus)
box(ax, (6.0, 5.05), 2.2, 0.65, r"NF-$\kappa$B",
    facecolor="#e7dbf5", fontsize=12)

# COX-2 (left branch)
box(ax, (3.2, 3.2), 1.8, 0.6, "COX-2",
    facecolor="#d9f0d9", fontsize=11)

# PGE2
box(ax, (3.2, 1.7), 1.8, 0.6, r"PGE$_2$",
    facecolor="#bfe3bf", fontsize=12)

# Cytokines (right branch)
box(ax, (8.8, 3.2), 2.4, 0.6,
    "Pro-inflammatory\ngene transcription",
    facecolor="#ffe3c2", fontsize=9.5, fontweight="normal")

box(ax, (8.8, 1.7), 3.0, 0.8,
    r"TNF-$\alpha$   IL-1$\beta$   IL-6",
    facecolor="#ffc9a8", fontsize=12)

# EP receptor (for PGE2 autocrine/paracrine loop)
box(ax, (6.0, 1.1), 2.2, 0.5, "EP1-4 receptors",
    facecolor="#eafbe4", fontsize=9.5, fontweight="normal")

# ---- Canonical signalling arrows --------------------------------------
arrow(ax, (6.0, 7.97), (6.0, 7.33))          # LPS -> TLR4
arrow(ax, (6.0, 6.77), (6.0, 6.43))          # TLR4 -> MyD88
arrow(ax, (6.0, 5.87), (6.0, 5.40))          # MyD88 -> NF-kB
# NF-kB -> COX-2 (left)
arrow(ax, (5.1, 4.85), (3.6, 3.52))
# COX-2 -> PGE2
arrow(ax, (3.2, 2.88), (3.2, 2.02))
# NF-kB -> cytokine gene transcription (right)
arrow(ax, (6.9, 4.85), (8.35, 3.52))
# transcription -> cytokines
arrow(ax, (8.8, 2.88), (8.8, 2.12))

# ---- Hypothesised modulatory link (PGE2 -> cytokine response) ---------
# PGE2 -> EP receptors (autocrine/paracrine)
arrow(ax, (3.95, 1.55), (5.05, 1.20),
      color="#2a7a2a", lw=1.6, linestyle="--")
# EP receptors -> cytokine output (modulation, dashed)
arrow(ax, (7.05, 1.25), (7.50, 1.55),
      color="#2a7a2a", lw=1.6, linestyle="--")

# Feedback from PGE2/EP to NF-kB activity (dashed, curved)
arrow(ax, (6.0, 1.38), (6.0, 4.70),
      color="#2a7a2a", lw=1.4, linestyle=":",
      connectionstyle="arc3,rad=0.4")
ax.text(1.6, 0.55,
        r"Hypothesised modulatory link: PGE$_2$ / EP signalling can"
        "\nenhance or dampen pro-inflammatory cytokine output",
        fontsize=9, color="#2a7a2a", style="italic")

# ---- COX inhibitor (e.g. carprofen) ------------------------------------
box(ax, (1.3, 3.2), 1.4, 0.55, "COX\ninhibitor",
    facecolor="#ffe1e1", edgecolor="#b30000",
    textcolor="#7a0000", fontsize=9.5)
inhibit_arrow(ax, (2.00, 3.2), (2.68, 3.2))
ax.text(1.3, 2.55, "(e.g. carprofen)",
        ha="center", fontsize=8, color="#7a0000", style="italic")

# ---- Title / caption ---------------------------------------------------
ax.text(6.0, 8.85,
        r"Hypothesised mechanism: COX inhibition / PGE$_2$ reduction"
        r" and cytokine responses",
        ha="center", fontsize=13, fontweight="bold", color="#222")

# Legend
legend_handles = [
    Line2D([0], [0], color="#222222", lw=1.8, label="Canonical signalling"),
    Line2D([0], [0], color="#2a7a2a", lw=1.6, linestyle="--",
           label=r"PGE$_2$ / EP-mediated modulation (hypothesised)"),
    Line2D([0], [0], color="#b30000", lw=2.0,
           label="Pharmacological inhibition"),
]
ax.legend(handles=legend_handles, loc="lower right",
          frameon=True, fontsize=8.5, bbox_to_anchor=(0.995, 0.01))

plt.tight_layout()
out = "cox_pge2_cytokine_schematic.png"
plt.savefig(out, dpi=300, bbox_inches="tight")
plt.savefig("cox_pge2_cytokine_schematic.pdf", bbox_inches="tight")
print(f"Saved: {out}")
