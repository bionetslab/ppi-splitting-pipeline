#!/usr/bin/env python3
"""Generate a subway-graph-style overview of the PPI splitting pipeline."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

OUT = "../pipeline_overview.png"

C_IN    = '#455A64'
C_FETCH = '#2E7D32'
C_MAIN  = '#1565C0'
C_ANA   = '#E65100'
C_HEAT  = '#00695C'
C_OUT   = '#6A1B9A'

R   = 0.32
LW  = 5.5
LWb = 3.8
MY  = 5.5   # main spine y

fig, ax = plt.subplots(figsize=(26, 11))
ax.set_xlim(0, 26)
ax.set_ylim(0, 11)
ax.axis('off')
fig.patch.set_facecolor('white')


def station(x, y, col):
    ax.add_patch(Circle((x, y), R,       color=col,     zorder=5))
    ax.add_patch(Circle((x, y), R * 0.5, color='white', zorder=6))


def lbl(x, y, text, col, side='above', sub=None):
    yo  = (R + 0.18)  if side == 'above' else -(R + 0.18)
    va  = 'bottom'    if side == 'above' else 'top'
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


# ── main spine nodes ──────────────────────────────────────────────────────────
# (x, label, colour, side, sublabel)
nodes = [
    (1.2,  'ppis.csv',            C_IN,    'below', 'input'),
    (3.0,  'FETCH\nDATA',         C_FETCH, 'above', 'seqs · GO · species'),
    (4.8,  'GET\nLENGTHS',        C_MAIN,  'below', None),
    (6.6,  'RUN\nBLAST',          C_MAIN,  'above', 'all-vs-all BLASTp'),
    (8.4,  'MAKE\nMETIS',         C_MAIN,  'below', None),
    (10.2, 'RUN\nKAHIP',          C_MAIN,  'above', 'graph partition'),
    (12.0, 'SORT\nPPIS',          C_MAIN,  'below', None),
    (13.8, 'CDHIT\n×2',           C_MAIN,  'above', 'train↔val · train↔test'),
    (15.6, 'REMOVE\nREDUNDANT',   C_MAIN,  'below', None),
    (17.4, 'SAMPLE\nNEGATIVES',   C_MAIN,  'above', 'balanced · realistic'),
    (19.2, 'EMBED\nSEQUENCES',    C_MAIN,  'below', 'esm2 · prot_t5 · …'),
]

# coloured spine segments
track(1.2,       MY, 3.0 + R, MY, C_IN,    lw=LW)
track(3.0 - R,   MY, 4.8 + R, MY, C_FETCH, lw=LW)
track(4.8 - R,   MY, 19.2,    MY, C_MAIN,  lw=LW)

for x, txt, col, side, sub in nodes:
    station(x, MY, col)
    lbl(x, MY, txt, col, side, sub)

# ── SIMILARITY_HEATMAP — branches UP from CDHIT (13.8) ───────────────────────
SH_X, SH_Y = 13.8, 8.8
track(SH_X, MY + R, SH_X, SH_Y, C_HEAT, lw=LWb)
station(SH_X, SH_Y, C_HEAT)
lbl(SH_X, SH_Y, 'SIMILARITY\nHEATMAP', C_HEAT, side='above')

# ── BIAS_ANALYSIS + COLLECT_BIAS — branch UP from EMBED_SEQUENCES (19.2) ─────
BA_X, BA_Y = 19.2, 8.8
CA_X, CA_Y = 21.2, 8.8

track(BA_X, MY + R, BA_X, BA_Y, C_ANA)
station(BA_X, BA_Y, C_ANA)
lbl(BA_X, BA_Y, 'BIAS\nANALYSIS ×6–7', C_ANA, side='above', sub='parallel')

track(BA_X + R, BA_Y, CA_X, BA_Y, C_ANA)
station(CA_X, CA_Y, C_ANA)
lbl(CA_X, CA_Y, 'COLLECT\nBIAS', C_ANA, side='above', sub='scatter plot')

# ── TRAIN_CLASSIFIER — branch DOWN from EMBED_SEQUENCES (19.2) ───────────────
TC_X, TC_Y = 21.2, 2.0
vcorner(19.2, MY - R, TC_X, TC_Y, C_ANA)
station(TC_X, TC_Y, C_ANA)
lbl(TC_X, TC_Y, 'TRAIN\nCLASSIFIER', C_ANA, side='below')

# ── converge all branches to MULTIQC ─────────────────────────────────────────
MQC_X = 24.2

# SIM_HEATMAP → right along y=SH_Y, then down to MY
hcorner(SH_X + R, SH_Y, MQC_X, MY + R, C_OUT, lw=LWb)

# COLLECT_BIAS → right along y=CA_Y, then down to MY
hcorner(CA_X + R, CA_Y, MQC_X, MY + R, C_OUT, lw=LW)

# TRAIN_CLASSIFIER → right along y=TC_Y, then up to MY
hcorner(TC_X + R, TC_Y, MQC_X, MY - R, C_OUT, lw=LW)

# main spine continues right to MULTIQC
track(19.2 + R, MY, MQC_X, MY, C_OUT, lw=LW)

station(MQC_X, MY, C_OUT)
lbl(MQC_X, MY, 'MULTIQC', C_OUT, side='above', sub='multiqc_report.html')

# ── legend (bottom-left) ──────────────────────────────────────────────────────
legend_items = [
    (C_IN,    'Input'),
    (C_FETCH, 'UniProt fetch'),
    (C_MAIN,  'Core pipeline'),
    (C_ANA,   'Analysis'),
    (C_HEAT,  'Similarity heatmap'),
    (C_OUT,   'Report'),
]
lx0, ly0 = 0.4, 4.0
ax.text(lx0, ly0 + 0.6, 'Legend', fontsize=18, fontweight='bold', color='#455A64')
for i, (col, lab) in enumerate(legend_items):
    y = ly0 - i * 0.62
    ax.add_patch(Circle((lx0 + 0.2, y), 0.18, color=col, zorder=5))
    ax.text(lx0 + 0.52, y, lab, va='center', fontsize=18,
            color=col, fontweight='bold')

# ── title ─────────────────────────────────────────────────────────────────────
ax.text(13.0, 10.4, 'PPI Splitting Pipeline', ha='center', va='center',
        fontsize=18, fontweight='bold', color='#1A237E')

plt.tight_layout(pad=0.3)
plt.savefig(OUT, dpi=180, bbox_inches='tight', facecolor='white')
print(f"Saved {OUT}")
