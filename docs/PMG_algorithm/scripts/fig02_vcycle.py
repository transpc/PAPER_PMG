"""Figure 02: V-cycle structure exactly as implemented in SOLVER_NEW (7_SOLVE_GMG.f90)."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
from style_common import setup, save, C_FINE, C_COARSE, C_GC, C_SEND

fig, ax = setup(figsize=(12, 7))
ax.set_xlim(-0.7, 12.2); ax.set_ylim(-0.6, 5.4); ax.axis("off")

# levels: y positions (level 1 top ... level L bottom -> GC)
LV = {1: 4.5, 2: 3.3, 3: 2.1, 4: 0.9}   # 4 shown local levels; level 4 = coarsest (GC)
GCY = 0.0

# node positions along the cycle
down_x = {1: 1.7, 2: 2.9, 3: 4.1}
gc_x = 5.8
up_x = {3: 7.4, 2: 8.8, 1: 10.2}

def node(x, y, color, r=0.22):
    ax.add_patch(Circle((x, y), r, fc=color, ec="#222222", zorder=5))

def arr(p, q, color="#444444", lw=1.6):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=13,
                                 color=color, lw=lw, zorder=4,
                                 shrinkA=10, shrinkB=10))

# level guide lines + labels
for lv, y in LV.items():
    ax.axhline(y, color="#DDDDDD", lw=0.8, zorder=0)
lab = {1: "level 1 (finest)\nlocal A, n = nintf",
       2: "level 2\nGalerkin A$_2$=R$_1$A$_1$P$_1$",
       3: "level ilv = 2..nlevel-1\n(loop over coarse levels)",
       4: "level nlevel (local coarsest)"}
for lv, y in LV.items():
    ax.text(-0.65, y, lab[lv], fontsize=8.5, va="center", ha="left", color="#333333")
ax.text(-0.65, GCY, "GLOBAL COARSE (GC)\ngathered on all ranks", fontsize=8.5,
        va="center", ha="left", color=C_GC)

# --- downward leg
node(down_x[1], LV[1], C_FINE)
ax.annotate("PRE-SMOOTH  itergs(1)x\npoly_cheb_smooth (Chebyshev)\n+ halo S&R(A) each sweep\nthen r = b - Au ; S&R(R) on r",
            (down_x[1], LV[1]), xytext=(down_x[1]-0.55, LV[1]+0.42), fontsize=8, color=C_FINE)

for lv in (2, 3):
    node(down_x[lv], LV[lv], C_COARSE)
    arr((down_x[lv-1], LV[lv-1]), (down_x[lv], LV[lv]))
ax.annotate("restrict rc = R rt (matrix_vec_N)\nsmooth e: itergs(ilv)x Smooth_GS2\n  [each sweep: S&R(A, level ilv)]\nrt = rc - A e ; S&R(R, level ilv)",
            (down_x[2], LV[2]), xytext=(down_x[2]-0.9, LV[2]-0.95), fontsize=8, color="#7A4A12")

# --- coarsest
node(gc_x, GCY, C_GC, r=0.28)
arr((down_x[3], LV[3]), (gc_x, GCY), color=C_GC)
ax.annotate("restrict to local coarsest rs\nMPI_ALLGATHERV -> global rG (all ranks)\nSOLVE_EXACT: eG = A$_{GC}^{-1}$ rG  (dense inverse, redundant)\nor SOLVE_COARSE: serial sub-V-cycle (nlv_glo>0)\ncopy back local es = eG(imapG)",
            (gc_x, GCY), xytext=(gc_x-1.3, -0.62+0.9), fontsize=8, color=C_GC,
            ha="left")

# --- upward leg
arr((gc_x, GCY), (up_x[3], LV[3]), color=C_GC)
for lv in (3, 2):
    node(up_x[lv], LV[lv], C_COARSE)
arr((up_x[3], LV[3]), (up_x[2], LV[2]))
arr((up_x[2], LV[2]), (up_x[1], LV[1]))
node(up_x[1], LV[1], C_FINE)
ax.annotate("prolong et = P e (matrix_vec_N)\ne += et ; then itergs(ilv)x:\n  S&R(A) -> Smooth_GS2\nfinally S&R(P) on e",
            (up_x[2], LV[2]), xytext=(up_x[2]-0.4, LV[2]-1.0), fontsize=8, color="#7A4A12")
ax.annotate("u += P e ; S&R(A)\nPOST-SMOOTH itergs(1)x Chebyshev\nresidual norm + MPI_ALLREDUCE\n(converged? res < crit -> exit)",
            (up_x[1], LV[1]), xytext=(up_x[1]-0.8, LV[1]+0.42), fontsize=8, color=C_FINE)

ax.set_title("One V-cycle of SOLVER_NEW (ncycle=1 per preconditioner application)   "
             "S&R = MPI halo exchange", fontsize=12.5, pad=14)

# legend
from matplotlib.lines import Line2D
leg = [Line2D([], [], marker="o", ls="", mfc=C_FINE, mec="k", label="finest level (Chebyshev smoothing)"),
       Line2D([], [], marker="o", ls="", mfc=C_COARSE, mec="k", label="coarse level (Gauss-Seidel smoothing)"),
       Line2D([], [], marker="o", ls="", mfc=C_GC, mec="k", label="global coarse (gathered serial solve)")]
ax.legend(handles=leg, loc="lower right", fontsize=8.5, frameon=False)
save(fig, "fig02_vcycle.png")
