"""Figure 01: overall solver component stack (call hierarchy) of the PMG-preconditioned BiCGSTAB."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from style_common import setup, save, C_FINE, C_COARSE, C_GC, C_SEND

fig, ax = setup(figsize=(11, 7.5))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, text, fc, fontsize=10, tc="white", lw=0.8):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06",
                       fc=fc, ec="#333333", lw=lw)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center",
            fontsize=fontsize, color=tc)

def arrow(x1, y1, x2, y2, color="#333333", style="-|>", ls="-"):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                        mutation_scale=14, color=color, lw=1.4, linestyle=ls)
    ax.add_patch(a)

# Level 0: entry
box(2.6, 9.0, 4.8, 0.8,
    "SOLVE_GMG(icase)\n[7_SOLVE_GMG.f90]  pressure eq. entry", "#37474F", fontsize=10)

# Level 1: setup-per-call + outer solver
box(0.3, 7.4, 3.4, 1.0,
    "icase=1 only:\nstiffness_MG (Galerkin RAP)\neig_value (Gershgorin bound)", C_COARSE, fontsize=9)
box(4.3, 7.4, 5.2, 1.0,
    "solve_pbcg_mg  [6_solver_pbcg_mg.f90]\nBiCGSTAB outer iteration (unpreconditioned ops on A)",
    C_FINE, fontsize=10)
arrow(4.0, 9.0, 2.2, 8.5)
arrow(5.8, 9.0, 6.6, 8.5)

# Level 2: pieces used by BiCGSTAB
box(0.3, 5.4, 2.9, 1.2,
    "amux0P\nSpMV  y = A x\n(fine level, CSR)", "#546E7A", fontsize=9)
box(3.5, 5.4, 2.9, 1.2,
    "communicate_s\nhalo exchange of x\n(MPI_ISEND/IRECV)", C_SEND, fontsize=9)
box(6.7, 5.4, 3.0, 1.2,
    "SOLVER_NEW  (x2 per iter.)\nPMG V-cycle preconditioner\ny = M$^{-1}$p ,  z = M$^{-1}$s", C_COARSE, fontsize=9)
arrow(5.5, 7.4, 1.8, 6.7)
arrow(6.6, 7.4, 5.0, 6.7)
arrow(7.6, 7.4, 8.2, 6.7)

# Level 3: V-cycle internals
box(0.3, 3.2, 2.6, 1.3,
    "smoothing_fine\npoly_cheb_smooth\n(4th-kind Chebyshev,\nfinest level only)", C_FINE, fontsize=8.5)
box(3.1, 3.2, 2.2, 1.3,
    "matrix_vec_N\nrestriction rc=R rt\nprolongation et=P e", "#546E7A", fontsize=8.5)
box(5.5, 3.2, 2.2, 1.3,
    "Smooth_GS2\nGauss-Seidel(Jacobi-\nhybrid) smoothing,\ncoarse levels 2..L-1", C_COARSE, fontsize=8.5)
box(7.9, 3.2, 1.9, 1.3,
    "SOLVE_GC_all\ncoarsest level\n(gathered, serial)", C_GC, fontsize=8.5)
for xx in (1.6, 4.2, 6.6, 8.85):
    arrow(8.2, 5.4, xx, 4.5)

# Level 4: communication used inside V-cycle
box(0.9, 1.0, 3.6, 1.2,
    "send_receive / send_receive_C\nfine: A-pattern halo (nnbdA...)\ncoarse: MD_S_R_NEW(id=1,2,3)\nA / R / P per-level lists", C_SEND, fontsize=8.5)
box(5.1, 1.0, 4.0, 1.2,
    "MPI_ALLGATHERV (residual gather)\nSOLVE_EXACT: e = A$_{GC}^{-1}$ r (dense inverse,\nredundant on every rank)", C_GC, fontsize=8.5)
arrow(1.6, 3.2, 2.2, 2.2)
arrow(6.6, 3.2, 5.6, 2.2)
arrow(8.85, 3.2, 7.6, 2.2)

ax.set_title("PMG-preconditioned BiCGSTAB: component stack (file / subroutine map)",
             fontsize=13, pad=12)
save(fig, "fig01_solver_stack.png")
