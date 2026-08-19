"""Figure 06: concatenated per-level storage of the hierarchy
(ialv / iintf / ncolf / ncolc offset arithmetic of MD_MG_coord)."""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from style_common import setup, save, C_FINE, C_COARSE, C_GC, C_GHOST

fig, ax = setup(figsize=(12, 6.2))
ax.set_xlim(0, 13); ax.set_ylim(-0.8, 6.6); ax.axis("off")

# level sizes (example): nnode=6 (fine), lev2=3.4, lev3=1.8, lev4=1.0 units wide
lv_owned = [4.2, 2.3, 1.2]
lv_ghost = [1.3, 0.8, 0.5]
names = ["level 1 (fine)", "level 2", "level 3 = nlevel"]

def bar(y, segs, label, h=0.62):
    x = 1.8
    ax.text(1.65, y + h / 2, label, ha="right", va="center", fontsize=10)
    xs = []
    for (w, fc, txt) in segs:
        ax.add_patch(Rectangle((x, y), w, h, fc=fc, ec="k", lw=0.7))
        if txt:
            ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center",
                    fontsize=7.6, color="white")
        xs.append(x)
        x += w
    xs.append(x)
    return xs

# --- ncolf-space arrays: rt, et, iai(P rows)  (levels 1..nlevel)
segs = []
cols = [C_FINE, C_COARSE, C_GC]
for i in range(3):
    segs.append((lv_owned[i], cols[i], f"owned\niintf({i+1})"))
    segs.append((lv_ghost[i], C_GHOST, "ghost"))
xs_f = bar(4.9, segs, "rt, et, iai/jai (P)\n[ncolf space]")
ax.text(xs_f[-1] + 0.15, 5.2, "ncolf =\nialv(nlevel+1)-1", fontsize=8.5, va="center")

# ialv pointers
acc = 1.8
ptr = [1.8]
for i in range(3):
    acc += lv_owned[i] + lv_ghost[i]
    ptr.append(acc)
for i, x in enumerate(ptr):
    lab_i = f"ialv({i+1})" if i < 3 else "ialv(nlevel+1)"
    ax.annotate(lab_i, (x, 5.55), xytext=(x - 0.35, 6.15),
                fontsize=9, color="#B03A48",
                arrowprops=dict(arrowstyle="->", color="#B03A48", lw=1.0))

# --- ncolc-space arrays: e, rc, iac/auc, diagrc, iar(R rows) (levels 2..nlevel)
segs = []
for i in range(1, 3):
    segs.append((lv_owned[i], cols[i], "owned"))
    segs.append((lv_ghost[i], C_GHOST, "ghost"))
x0 = 1.8 + lv_owned[0] + lv_ghost[0]
# shift: draw aligned under the coarse part of the top bar
ax.text(1.65, 3.0 + 0.31, "e, rc, auc/iac (Ac),\ndiagrc, iar/jar (R)\n[ncolc space]",
        ha="right", va="center", fontsize=10)
x = x0
for (w, fc, txt) in segs:
    ax.add_patch(Rectangle((x, 3.0), w, 0.62, fc=fc, ec="k", lw=0.7))
    if txt:
        ax.text(x + w / 2, 3.31, txt, ha="center", va="center", fontsize=7.6, color="white")
    x += w
ax.text(x + 0.15, 3.31, "ncolc = ncolf - nnode", fontsize=8.5, va="center")

# alignment arrows between the two spaces
ax.add_patch(FancyArrowPatch((x0, 4.9), (x0, 3.62), arrowstyle="-|>",
                             mutation_scale=12, color="#555", lw=1.2, linestyle="--"))
ax.text(x0 + 0.08, 4.15, "index shift by nnode:\nncolc idx = ncolf idx - nnode\n"
        "(e.g. ista = ialv(ilv) - nnode)", fontsize=9, color="#333")

# --- fine-level arrays
ax.text(1.65, 1.55 + 0.31, "u, b, r, diagr, au/ia (A)\n[fine space]",
        ha="right", va="center", fontsize=10)
ax.add_patch(Rectangle((1.8, 1.55), lv_owned[0], 0.62, fc=C_FINE, ec="k", lw=0.7))
ax.text(1.8 + lv_owned[0]/2, 1.86, "owned  1..nintf", ha="center", va="center",
        fontsize=8, color="white")
ax.add_patch(Rectangle((1.8 + lv_owned[0], 1.55), lv_ghost[0], 0.62, fc=C_GHOST, ec="k", lw=0.7))
ax.text(1.8 + lv_owned[0] + lv_ghost[0]/2, 1.86, "ghost\nnintf+1..nnode",
        ha="center", va="center", fontsize=7.2, color="white")
ax.text(1.8 + lv_owned[0] + lv_ghost[0] + 0.15, 1.86,
        "(+ 2nd halo layer up to nnodegl,\nused only for Galerkin RAP columns)",
        fontsize=8.5, va="center", color="#555")

# --- ownership note
ax.text(1.8, 0.45,
        "Within each level segment:  owned nodes first (1..iintf), ghost nodes last, ghosts grouped\n"
        "in contiguous blocks per neighbor rank (in nbdom order) - this makes rpt offsets directly\n"
        "addressable by MPI_IRECV without extra copies.",
        fontsize=9.5, color="#333")

ax.set_title("Concatenated level storage: all coarse levels share single CSR / vector arrays\n"
             "(MD_MG_coord: ialv offsets;  2_read_mesh_MPI.f90 allocation)",
             fontsize=12, pad=10)
save(fig, "fig06_memory_layout.png")
