"""Figure 12: worked 1D example of the parallel Galerkin RAP (stiff_coarse_P).

6 fine nodes, 2 ranks, A = tridiag(-1,2,-1), C points at fine 1,3,5.
All numbers in the figure are COMPUTED here with numpy (not hand-typed),
so the picture is guaranteed consistent with Ac = P^T A P.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from style_common import setup, save, RANK_COLORS, RANK_LIGHT, C_GC, C_COARSE

# ---------------- exact operators ----------------
n = 6
A = np.zeros((n, n))
for i in range(n):
    A[i, i] = 2.0
    if i > 0: A[i, i-1] = -1.0
    if i < n-1: A[i, i+1] = -1.0

P = np.zeros((n, 3))
P[0, 0] = 1.0
P[1, 0] = P[1, 1] = 0.5
P[2, 1] = 1.0
P[3, 1] = P[3, 2] = 0.5
P[4, 2] = 1.0
P[5, 2] = 1.0          # only one C neighbor -> injection weight 1
R = P.T
Ac = R @ A @ P

# vi accumulation for coarse row I=2 (0-based 1): k in {2,3,4} (1-based) w/ weights R[1,k]
ks = [1, 2, 3]                     # 0-based fine rows 2,3,4
wk = [R[1, k] for k in ks]
contrib = [w * A[k, :] for w, k in zip(wk, ks)]
vi = np.sum(contrib, axis=0)
# contraction checks
assert abs(vi @ P[:, 0] - Ac[1, 0]) < 1e-14
assert abs(vi @ P[:, 1] - Ac[1, 1]) < 1e-14
assert abs(vi @ P[:, 2] - Ac[1, 2]) < 1e-14

owner_f = [0, 0, 0, 1, 1, 1]       # fine ownership
owner_c = [0, 0, 1]                # coarse ownership (inherits C-point owner)

# ---------------- figure ----------------
fig = plt.figure(figsize=(13.5, 10))
plt.rcParams.update({"font.size": 11, "font.family": "DejaVu Sans"})

# ===== panel A: grid + operators + dependencies =====
axA = fig.add_axes([0.04, 0.56, 0.92, 0.40]); axA.axis("off")
axA.set_xlim(0.2, 6.9); axA.set_ylim(-2.7, 2.6)

Cpts = [0, 2, 4]
for i in range(n):
    fc = RANK_LIGHT[owner_f[i]]
    if i in Cpts:
        axA.add_patch(Rectangle((i+1-0.17, -0.17), 0.34, 0.34, fc=fc,
                                ec=RANK_COLORS[owner_f[i]], lw=2.4, zorder=5))
    else:
        axA.add_patch(plt.Circle((i+1, 0), 0.15, fc=fc,
                                 ec="#888", lw=1.0, zorder=5))
    axA.text(i+1, -0.52, f"fine {i+1}", ha="center", fontsize=9)
axA.plot([1, 6], [0, 0], color="#CCCCCC", lw=1.2, zorder=1)
axA.axvline(3.5, color="k", lw=1.2, ls="--")
axA.text(2.0, 2.35, "rank 0  (owns fine 1-3, coarse C1,C2)", color=RANK_COLORS[0],
         fontsize=10.5, ha="center")
axA.text(5.1, 2.35, "rank 1  (owns fine 4-6, coarse C3)", color=RANK_COLORS[1],
         fontsize=10.5, ha="center")
# C labels
for ci, i in enumerate(Cpts):
    axA.text(i+1, 0.34, f"C{ci+1}", ha="center", fontsize=10, weight="bold",
             color=RANK_COLORS[owner_f[i]])
# P arrows (from C to F) with weights
for f in range(n):
    for c in range(3):
        if P[f, c] > 0 and f != Cpts[c]:
            src, dst = Cpts[c]+1, f+1
            axA.add_patch(FancyArrowPatch((src, 0.16), (dst, 0.16),
                          connectionstyle="arc3,rad=-0.45", arrowstyle="-|>",
                          mutation_scale=9, color="#9AA5B1", lw=1.1, zorder=3))
            axA.text((src+dst)/2, 0.78, f"{P[f,c]:.1f}", ha="center",
                     fontsize=8, color="#6B7682")
axA.text(0.35, 0.95, "P (weights)", fontsize=9, color="#6B7682")

# R-row of C2 brace: weights over fine 2,3,4
yb = -1.05
axA.plot([2, 4], [yb, yb], color=C_GC, lw=2.0)
for k, w in zip([2, 3, 4], wk):
    axA.plot([k, k], [yb, -0.30], color=C_GC, lw=1.2, ls=":")
    axA.text(k, yb-0.32, f"{w:.1f}", ha="center", fontsize=9.5, color=C_GC)
axA.text(3.0, yb-0.75, "R-row of C2 = P column 2 transposed:  fine {2, 3, 4}",
         color=C_GC, fontsize=10, ha="center")

# ghost row import annotation
axA.annotate("A(4,:) = ghost ROW for rank 0\nvalues imported from rank 1\n(send_receive_mtf, before RAP)",
             (4, 0.18), xytext=(4.75, 1.35), fontsize=9, color=C_GC,
             arrowprops=dict(arrowstyle="->", color=C_GC))
# 2nd halo annotation
axA.annotate("fine 5 = 2nd-halo index (nnodegl):\nappears only as a COLUMN of ghost row A(4,:)\n-> needs a local index for vi, but no value",
             (5, -0.2), xytext=(4.9, -1.9), fontsize=9, color="#7A4A12",
             arrowprops=dict(arrowstyle="->", color="#7A4A12"))
# rank0 local numbering
axA.text(0.35, 1.55, "rank 0 local ids:", fontsize=8.5, color="#555")
loc0 = ["1", "2", "3", "4 (ghost)", "5 (2nd halo)", "-"]
for i, t in enumerate(loc0):
    axA.text(i+1, 1.25, t, ha="center", fontsize=8, color="#555")
axA.set_title("Parallel RAP, worked 1D example:  rank 0 computes ITS OWNED coarse row Ac(C2, :)\n"
              r"$A_c(I,J)=\sum_k R(I,k)\,\sum_l A(k,l)\,P(l,J)$   with  $A$ = tridiag(-1, 2, -1),  C points at fine 1, 3, 5",
              fontsize=12)

# ===== panel B: vi accumulation table =====
axB = fig.add_axes([0.115, 0.06, 0.44, 0.44]); axB.axis("off")
rows = []
labels = []
for k, w, c in zip([2, 3, 4], wk, contrib):
    labels.append(f"{w:.1f} x A({k},:)")
    rows.append([f"{v:+.1f}" if abs(v) > 1e-14 else "·" for v in c[:5]])
labels.append("vi (sum)")
rows.append([f"{v:+.1f}" if abs(v) > 1e-14 else "·" for v in vi[:5]])
tbl = axB.table(cellText=rows, rowLabels=labels,
                colLabels=[f"fine {j+1}" for j in range(5)],
                loc="center", cellLoc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(10); tbl.scale(1, 1.7)
# highlight ghost-row contribution (row index 2 -> table row 3) and 2nd-halo column (col 4)
for j in range(5):
    tbl[3, j].set_facecolor("#F6E3D3")          # ghost row A(4,:) contribution
for i in range(1, 5):
    tbl[i, 4].set_facecolor("#F3D9DB")          # column fine5 (2nd halo)
tbl[3, 4].set_facecolor("#EFC3C7")
axB.set_title("STEP 1 - dense accumulator (code: vi):\n"
              r"$v_i(\cdot) = \sum_{k \in R\mathrm{-row}(C2)} R(C2,k)\; A(k,\cdot)$"
              "\norange = uses imported ghost row  |  red = lands on 2nd-halo index",
              fontsize=10.5)

# ===== panel C: contraction + Ac ownership =====
axC = fig.add_axes([0.58, 0.06, 0.40, 0.44]); axC.axis("off")
axC.set_xlim(0, 10); axC.set_ylim(0, 10)
lines = [
    r"STEP 2 - contract vi with R-rows of neighbor coarse J:",
    r"$A_c(C2,J) = \sum_l v_i(l)\, R(J,l)$",
    "",
    f"J=C1 (owned):    1.0·vi(1) + 0.5·vi(2)              = {Ac[1,0]:+.1f}",
    f"J=C2 (owned):    0.5·vi(2) + 1.0·vi(3) + 0.5·vi(4)  = {Ac[1,1]:+.1f}",
    f"J=C3 (GHOST coarse): 0.5·vi(4) + 1.0·vi(5) + 1.0·vi(6)",
    f"                 vi(6): no local index -> nj<=0, SKIPPED",
    f"                 (safe: vi(6)=0 by construction)     = {Ac[1,2]:+.1f}",
]
for i, t in enumerate(lines):
    axC.text(0.2, 9.4 - i*0.62, t, fontsize=9.6, family="monospace" if i >= 3 else None,
             color=C_GC if "GHOST" in t or "SKIP" in t or "vi(6)" in t else "#222")
# Ac matrix with row ownership
y0 = 3.9
axC.text(0.2, y0 + 0.75, "Result  Ac = P$^T$AP  (each row complete on its OWNER, no reduction):",
         fontsize=10)
for i in range(3):
    for j in range(3):
        v = Ac[i, j]
        axC.add_patch(Rectangle((1.4 + j*1.5, y0 - i*0.85 - 0.55), 1.5, 0.85,
                                fc=RANK_LIGHT[owner_c[i]], ec="#888", lw=0.7))
        axC.text(1.4 + j*1.5 + 0.75, y0 - i*0.85 - 0.12,
                 f"{v:+.1f}" if abs(v) > 1e-14 else "0",
                 ha="center", va="center", fontsize=10.5)
    axC.text(6.2, y0 - i*0.85 - 0.12, f"row C{i+1}: rank {owner_c[i]} computes",
             fontsize=9.5, va="center", color=RANK_COLORS[owner_c[i]])
axC.text(0.2, 0.35, "Ghost coarse rows of Ac are NOT computed here - they are imported\n"
                    "by MD_S_R_MT right before the NEXT level's RAP (recursion).",
         fontsize=9.5, color="#444")

save(fig, "fig12_rap_example.png")
print("Ac =\n", Ac)
