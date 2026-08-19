"""Figure 04: one outer BiCGSTAB iteration (solve_pbcg_mg) with the two PMG
preconditioner applications and every MPI synchronization point marked."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from style_common import setup, save, C_FINE, C_COARSE, C_SEND, C_GC

fig, ax = setup(figsize=(9.5, 11))
ax.set_xlim(0, 10); ax.set_ylim(0, 15.2); ax.axis("off")

steps = [
    # (text, color, kind)  kind: local / mpi / mg
    (r"$\rho = (\hat r_0, r)$   + MPI_ALLREDUCE", "mpi"),
    (r"breakdown guard ($|\rho_{old}|,|\omega| < 10^{-300}$ ?)"
     "\n" r"$\beta = (\rho/\rho_{old})(\alpha/\omega)$;   $p = r + \beta(p - \omega v)$", "local"),
    ("PMG:  u=y (warm start), b = scale_mg * p\n"
     r"SOLVER_NEW (1 V-cycle)  ->  $y = M^{-1}p$ / scale_mg", "mg"),
    ("communicate_s(y)   halo exchange (fine A-pattern)", "halo"),
    (r"$v = A\,y$   (amux0P, local rows 1..nintf)", "local"),
    (r"$\alpha_d=(\hat r_0, v)$  + MPI_ALLREDUCE;   $\alpha = \rho/\alpha_d$", "mpi"),
    (r"$s = r - \alpha v$", "local"),
    ("PMG:  u=z (warm start), b = scale_mg * s\n"
     r"SOLVER_NEW (1 V-cycle)  ->  $z = M^{-1}s$ / scale_mg", "mg"),
    ("communicate_s(z)   halo exchange", "halo"),
    (r"$t = A\,z$   (amux0P)", "local"),
    (r"$\omega = (t,s)/(t,t)$   + 2x MPI_ALLREDUCE", "mpi"),
    (r"$x = x + \alpha y + \omega z$;   $r = s - \omega t$;   $\|r\|$ + MPI_ALLREDUCE", "mpi"),
    (r"NaN guard; converged if $\|r\|/\|r_0\| \leq$ crit_bcg_mg", "local"),
]

colors = {"local": "#546E7A", "mpi": C_GC, "mg": C_COARSE, "halo": C_SEND}
label  = {"local": "local compute", "mpi": "global reduction (blocking)",
          "mg": "PMG V-cycle preconditioner", "halo": "neighbor halo exchange"}

y = 14.4
H = {"local": 0.75, "mpi": 0.75, "mg": 0.95, "halo": 0.6}
prev_bottom = None
for text, kind in steps:
    h = H[kind]
    y -= h + 0.28
    p = FancyBboxPatch((1.6, y), 6.8, h, boxstyle="round,pad=0.05",
                       fc=colors[kind], ec="#333", lw=0.7)
    ax.add_patch(p)
    ax.text(5.0, y + h/2, text, ha="center", va="center", fontsize=9.2, color="white")
    if prev_bottom is not None:
        ax.add_patch(FancyArrowPatch((5.0, prev_bottom), (5.0, y + h),
                                     arrowstyle="-|>", mutation_scale=12,
                                     color="#333", lw=1.1, shrinkA=1, shrinkB=1))
    prev_bottom = y

# loop-back arrow
ax.add_patch(FancyArrowPatch((1.6, y + 0.35), (0.85, y + 0.35),
                             arrowstyle="-", color="#333", lw=1.1))
ax.add_patch(FancyArrowPatch((0.85, y + 0.35), (0.85, 13.9),
                             arrowstyle="-", color="#333", lw=1.1))
ax.add_patch(FancyArrowPatch((0.85, 13.9), (1.6, 13.9),
                             arrowstyle="-|>", mutation_scale=12, color="#333", lw=1.1))
ax.text(0.62, 7.5, "next iteration (its < maxit)", rotation=90,
        fontsize=9, va="center", color="#333")

# legend
yl = 14.9
for kind in ["local", "halo", "mpi", "mg"]:
    ax.add_patch(FancyBboxPatch((1.6 + ["local","halo","mpi","mg"].index(kind)*2.2, yl-0.12),
                                0.32, 0.32, boxstyle="round,pad=0.02",
                                fc=colors[kind], ec="#333", lw=0.5))
    ax.text(2.0 + ["local","halo","mpi","mg"].index(kind)*2.2, yl+0.04, label[kind],
            fontsize=7.8, va="center")

ax.set_title("One BiCGSTAB iteration in solve_pbcg_mg\n"
             "(2 V-cycles + 5 ALLREDUCEs + 2 halo exchanges per iteration)",
             fontsize=12, pad=16)
save(fig, "fig04_bicgstab_flow.png")
