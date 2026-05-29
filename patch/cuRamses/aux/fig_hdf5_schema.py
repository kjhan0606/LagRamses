import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

fig, ax = plt.subplots(1, 1, figsize=(3.4, 5.0))

# Colors
C_ROOT = '#1E293B'
C_GRP  = '#1D4ED8'
C_DSET = '#059669'
C_LV   = '#7C3AED'
C_LINE = '#94A3B8'

LW = 1.5  # line weight for tree connectors

# Approximate monospace char width in data coords
# figsize 3.4", xlim [0, 2.8] → 2.8 data units = 244.8 pt
def char_w(fs):
    return 0.6 * fs * 2.8 / 244.8

# X positions
X_TRUNK = 0.10   # main vertical trunk
X_GRP   = 0.22   # group text
X_LV    = 0.78   # level text (and non-level datasets)
X_DS    = 1.15   # datasets under level groups

# Font sizes
FS_ROOT = 9
FS_GRP  = 8
FS_LV   = 7
FS_DS   = 7

y = 4.6  # start y
GAP = 0.10  # vertical gap below text before connector starts

def hline(x1, x2, yy):
    ax.plot([x1, x2], [yy, yy], color=C_LINE, lw=LW, solid_capstyle='round')

def vline(x, y1, y2):
    ax.plot([x, x], [y1, y2], color=C_LINE, lw=LW, solid_capstyle='round')

def txt(x, yy, text, color, fs):
    ax.text(x, yy, text, ha='left', va='center', fontsize=fs,
            color=color, fontweight='bold', fontfamily='monospace')

# ============ Root box ============
rw = len('output_XXXXX.h5') * char_w(FS_ROOT) + 0.15
r = mpatches.FancyBboxPatch((X_TRUNK - 0.03, y - 0.14), rw, 0.28,
    boxstyle='round,pad=0.03', fc='#E2E8F0', ec=C_ROOT, lw=1.2)
ax.add_patch(r)
txt(X_TRUNK + 0.02, y, 'output_XXXXX.h5', C_ROOT, FS_ROOT)
y_root = y - 0.14

# ============ Build tree data ============
groups = [
    ('/amr/',       True,  ['xg (positions)', 'son (child flags)', 'cpu_map']),
    ('/hydro/',     True,  [u'\u03c1, \u03c1v, E, Z, ...']),
    ('/gravity/',   True,  [u'\u03c6 (potential), f (force)']),
    ('/particles/', False, ['xp, vp (pos, vel)', 'mp, idp (mass, ID)',
                             'tp, zp (birth, Z)']),
    ('/sinks/',     False, ['msink, xsink, vsink, ...']),
]

# Pre-compute y positions
dy_grp  = 0.40
dy_item = 0.24
positions = []
cy = y - 0.45

for name, has_level, datasets in groups:
    grp_y = cy
    lv_y = None
    ds_list = []
    if has_level:
        cy -= dy_item
        lv_y = cy
        for ds in datasets:
            cy -= dy_item
            ds_list.append((cy, ds))
    else:
        for ds in datasets:
            cy -= dy_item
            ds_list.append((cy, ds))
    positions.append((grp_y, lv_y, ds_list))
    cy -= dy_grp

y_last_grp = positions[-1][0]

# ============ Draw main trunk ============
vline(X_TRUNK, y_root, y_last_grp)

# ============ Draw each group ============
for i, (grp_y, lv_y, ds_list) in enumerate(positions):
    name, has_level, _ = groups[i]

    # Horizontal stub: main trunk → group name
    hline(X_TRUNK, X_GRP - 0.02, grp_y)
    txt(X_GRP, grp_y, name, C_GRP, FS_GRP)

    # Center x of group name
    xm = X_GRP + len(name) * char_w(FS_GRP) / 2

    if has_level:
        # ㄴ from group center down to level
        vline(xm, grp_y - GAP, lv_y)
        hline(xm, X_LV - 0.02, lv_y)
        txt(X_LV, lv_y, 'level_{l}/', C_LV, FS_LV)

        # ㄴ from level center down to datasets
        lv_xm = X_LV + len('level_{l}/') * char_w(FS_LV) / 2
        if ds_list:
            vline(lv_xm, lv_y - GAP, ds_list[-1][0])
            for ds_y, ds_text in ds_list:
                hline(lv_xm, X_DS - 0.02, ds_y)
                txt(X_DS, ds_y, ds_text, C_DSET, FS_DS)
    else:
        # ㄴ from group center down to datasets
        if ds_list:
            vline(xm, grp_y - GAP, ds_list[-1][0])
            for ds_y, ds_text in ds_list:
                hline(xm, X_LV - 0.02, ds_y)
                txt(X_LV, ds_y, ds_text, C_DSET, FS_DS)

# ============ Legend (color swatches) ============
ly = cy + 0.20
sw = 0.12  # swatch width
sh = 0.12  # swatch height
for sx, sc, sl in [(0.10, C_GRP, 'Group'),
                    (0.85, C_DSET, 'Dataset'),
                    (1.75, C_LV, 'Per-level')]:
    r = mpatches.Rectangle((sx, ly - sh/2), sw, sh, fc=sc, ec='none', alpha=0.25)
    ax.add_patch(r)
    ax.text(sx + sw/2, ly, sl, ha='center', va='center',
            fontsize=6, color=sc, fontweight='bold')
ax.text(0.10, ly - 0.22, 'l = levelmin, ..., levelmax',
        fontsize=6, color=C_LINE, style='italic')

ax.set_xlim(-0.05, 2.8)
ax.set_ylim(ly - 0.40, y + 0.25)
ax.axis('off')

plt.tight_layout()
plt.savefig('/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/misc/fig_hdf5_schema.pdf',
            dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("fig_hdf5_schema.pdf generated successfully")
