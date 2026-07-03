#!/usr/bin/env python3
"""Generate a subway-graph-style overview of the PPI splitting pipeline."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

OUT = "../pipeline_overview.png"

C_IN    = '#455A64'
C_FETCH = '#2E7D32'
C_MAIN  = '#1565C0'  # KaHIP path + common spine
C_ILP   = '#AD1457'  # ILP path
C_ANA   = '#E65100'
C_HEAT  = '#00695C'
C_OUT   = '#6A1B9A'

R   = 0.32
LW  = 5.5
LWb = 3.8

MY    = 6.5   # main spine y
KH_Y  = 9.5   # KaHIP branch y
ILP_Y = 3.5   # ILP branch y
BA_Y  = 12.0  # bias analysis / heatmap y
TC_Y  = 1.2   # train classifier y

# x-coordinates
X_PPIS    = 1.0
X_FETCH   = 3.0
X_BLAST   = 5.2

X_LENGTHS = 7.2   # common prefix: shared by both split methods
X_METIS   = 9.4

X_KAHIP_A = 11.6  # KaHIP branch: RUN_KAHIP (k=3)
X_SORT    = 13.8

X_KAHIP_B = 11.6  # ILP branch: RUN_KAHIP (k=100)
X_SOLVE   = 13.8

X_CDHIT2D = 15.8  # common suffix
X_REMOVE  = 17.8
X_SAMPLE  = 19.8
X_EMBED   = 21.8
X_MULTIQC = 29.5

X_BIAS    = 21.8
X_COLLECT = 24.0
X_CLF     = 24.0

fig, ax = plt.subplots(figsize=(36, 14))
ax.set_xlim(0, 36)
ax.set_ylim(0, 14)
ax.axis('off')
fig.patch.set_facecolor('white')


def station(x, y, col):
    ax.add_patch(Circle((x, y), R,       color=col,     zorder=5))
    ax.add_patch(Circle((x, y), R * 0.5, color='white', zorder=6))


def lbl(x, y, text, col, side='above', sub=None):
    yo = (R + 0.18) if side == 'above' else -(R + 0.18)
    va = 'bottom'   if side == 'above' else 'top'
    ax.text(x, y + yo, text, ha='center', va=va,
            fontsize=18, fontweight='bold', color=col, linespacing=1.3)
    if sub:
        yso = yo + (0.55 if side == 'above' else -0.55)
        ax.text(x, y + yso, sub, ha='center', va=va,
                fontsize=16.5, color='#90A4AE', style='italic')


def track(x1, y1, x2, y2, col, lw=LW):
    ax.plot([x1, x2], [y1, y2], color=col, linewidth=lw,
            solid_capstyle='round', zorder=2)


def hcorner(x1, y1, x2, y2, col, lw=LW):
    """Horizontal first, then vertical."""
    ax.plot([x1, x2, x2], [y1, y1, y2], color=col, linewidth=lw,
            solid_capstyle='round', solid_joinstyle='round', zorder=2)


def vcorner(x1, y1, x2, y2, col, lw=LW):
    """Vertical first, then horizontal."""
    ax.plot([x1, x1, x2], [y1, y2, y2], color=col, linewidth=lw,
            solid_capstyle='round', solid_joinstyle='round', zorder=2)


# ── common prefix: ppis.csv → FETCH_DATA → RUN_BLAST → GET_LENGTHS → MAKE_METIS
track(X_PPIS,        MY, X_FETCH - R, MY, C_IN)
track(X_FETCH - R,   MY, X_BLAST - R, MY, C_FETCH)
track(X_BLAST - R,   MY, X_METIS + R, MY, C_MAIN)

station(X_PPIS,  MY, C_IN)
lbl(X_PPIS,  MY, 'ppis.csv',    C_IN,    side='below', sub='input')
station(X_FETCH, MY, C_FETCH)
lbl(X_FETCH, MY, 'FETCH\nDATA', C_FETCH, side='above', sub='seqs · GO · species')
station(X_BLAST, MY, C_MAIN)
lbl(X_BLAST, MY, 'RUN\nBLAST',  C_MAIN,  side='below', sub='all-vs-all BLASTp')
station(X_LENGTHS, MY, C_MAIN)
lbl(X_LENGTHS, MY, 'GET\nLENGTHS', C_MAIN, side='above')
station(X_METIS, MY, C_MAIN)
lbl(X_METIS, MY, 'MAKE\nMETIS', C_MAIN, side='below', sub='similarity graph')

# ── KaHIP branch (up): RUN_KAHIP (k=3) → SORT_PPIS ───────────────────────────
vcorner(X_METIS + R, MY, X_KAHIP_A - R, KH_Y, C_MAIN)
track(X_KAHIP_A - R, KH_Y, X_SORT + R, KH_Y, C_MAIN)
hcorner(X_SORT + R, KH_Y, X_CDHIT2D, MY + R, C_MAIN)

for x, txt, sub, side in [
    (X_KAHIP_A, 'RUN\nKAHIP', 'graph partition · k=3', 'above'),
    (X_SORT,    'SORT\nPPIS', None,                     'below'),
]:
    station(x, KH_Y, C_MAIN)
    lbl(x, KH_Y, txt, C_MAIN, side=side, sub=sub)

ax.text((X_KAHIP_A + X_SORT) / 2, KH_Y + 1.55,
        'KaHIP path', ha='center', va='center',
        fontsize=17, color=C_MAIN, fontstyle='italic', fontweight='bold')

# ── ILP branch (down): RUN_KAHIP (k=100) → SOLVE_ILP ─────────────────────────
vcorner(X_METIS + R, MY, X_KAHIP_B - R, ILP_Y, C_ILP)
track(X_KAHIP_B - R, ILP_Y, X_SOLVE + R, ILP_Y, C_ILP)
hcorner(X_SOLVE + R, ILP_Y, X_CDHIT2D, MY - R, C_ILP)

for x, txt, sub, side in [
    (X_KAHIP_B, 'RUN\nKAHIP',   'graph partition · k=100',    'below'),
    (X_SOLVE,   'SOLVE\nILP',   'gurobi license optional',    'above'),
]:
    station(x, ILP_Y, C_ILP)
    lbl(x, ILP_Y, txt, C_ILP, side=side, sub=sub)

ax.text((X_KAHIP_B + X_SOLVE) / 2, ILP_Y - 1.45,
        'ILP path', ha='center', va='center',
        fontsize=17, color=C_ILP, fontstyle='italic', fontweight='bold')

# ── common suffix: CDHIT2D → … → EMBED_SEQUENCES ─────────────────────────────
track(X_CDHIT2D + R, MY, X_EMBED, MY, C_MAIN)

for x, txt, sub, side in [
    (X_CDHIT2D, 'CDHIT2D\n×2',      'train↔val · train↔test', 'above'),
    (X_REMOVE,  'REMOVE\nREDUNDANT', None,                      'below'),
    (X_SAMPLE,  'SAMPLE\nNEGATIVES', 'balanced · realistic',    'above'),
    (X_EMBED,   'EMBED\nSEQUENCES',  'esm2 · prot_t5 · …',     'below'),
]:
    station(x, MY, C_MAIN)
    lbl(x, MY, txt, C_MAIN, side=side, sub=sub)

# ── SIMILARITY_HEATMAP — branches UP from CDHIT2D ────────────────────────────
track(X_CDHIT2D, MY + R, X_CDHIT2D, BA_Y, C_HEAT, lw=LWb)
station(X_CDHIT2D, BA_Y, C_HEAT)
lbl(X_CDHIT2D, BA_Y, 'SIMILARITY\nHEATMAP', C_HEAT, side='above')

# ── BIAS_ANALYSIS + COLLECT_BIAS — branch UP from EMBED_SEQUENCES ─────────────
track(X_BIAS, MY + R, X_BIAS, BA_Y, C_ANA)
station(X_BIAS, BA_Y, C_ANA)
lbl(X_BIAS, BA_Y, 'BIAS\nANALYSIS ×6–7', C_ANA, side='above', sub='parallel')

track(X_BIAS + R, BA_Y, X_COLLECT, BA_Y, C_ANA)
station(X_COLLECT, BA_Y, C_ANA)
lbl(X_COLLECT, BA_Y, 'COLLECT\nBIAS', C_ANA, side='above', sub='scatter plot')

# ── TRAIN_CLASSIFIER — branches DOWN from EMBED_SEQUENCES ─────────────────────
vcorner(X_EMBED, MY - R, X_CLF, TC_Y, C_ANA)
station(X_CLF, TC_Y, C_ANA)
lbl(X_CLF, TC_Y, 'TRAIN\nCLASSIFIER', C_ANA, side='below')

# ── all branches converge to MULTIQC ─────────────────────────────────────────
hcorner(X_CDHIT2D + R, BA_Y, X_MULTIQC, MY + R, C_OUT, lw=LWb)   # heatmap
hcorner(X_COLLECT + R, BA_Y, X_MULTIQC, MY + R, C_OUT, lw=LW)     # collect bias
hcorner(X_CLF + R,     TC_Y, X_MULTIQC, MY - R, C_OUT, lw=LW)     # classifier
track(X_EMBED + R,     MY,   X_MULTIQC, MY,      C_OUT, lw=LW)    # main spine

station(X_MULTIQC, MY, C_OUT)
lbl(X_MULTIQC, MY, 'MULTIQC', C_OUT, side='above', sub='multiqc_report.html')

# ── legend ────────────────────────────────────────────────────────────────────
legend_items = [
    (C_IN,    'Input'),
    (C_FETCH, 'UniProt fetch'),
    (C_MAIN,  'KaHIP path'),
    (C_ILP,   'ILP path'),
    (C_ANA,   'Analysis'),
    (C_HEAT,  'Similarity heatmap'),
    (C_OUT,   'Report'),
]
lx0, ly0 = 0.4, 6.0
ax.text(lx0, ly0 + 0.6, 'Legend', fontsize=18, fontweight='bold', color='#455A64')
for i, (col, lab) in enumerate(legend_items):
    y = ly0 - i * 0.62
    ax.add_patch(Circle((lx0 + 0.2, y), 0.18, color=col, zorder=5))
    ax.text(lx0 + 0.52, y, lab, va='center', fontsize=18,
            color=col, fontweight='bold')

# ── title ─────────────────────────────────────────────────────────────────────
ax.text(18.0, 13.3, 'PPI Splitting Pipeline', ha='center', va='center',
        fontsize=18, fontweight='bold', color='#1A237E')

plt.tight_layout(pad=0.3)
plt.savefig(OUT, dpi=180, bbox_inches='tight', facecolor='white')
print(f"Saved {OUT}")