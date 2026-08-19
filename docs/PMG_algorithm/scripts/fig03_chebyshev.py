"""Figure 03: error-damping polynomial of the 4th-kind Chebyshev smoother
as literally implemented in poly_cheb_smooth (poly_smooth.f90, method-2, a=0.3).

We reproduce the recurrence numerically for a scalar 'matrix' A = lambda:
starting from an error e0 = 1 (b = 0), the smoother update gives the
amplification factor |e_k(lambda)| after max_iter iterations.
"""
import numpy as np
from style_common import setup, save, C_FINE, C_COARSE, C_GC

def cheb_damp(lam_ratio, max_iter, a=0.3):
    """lam_ratio = lambda / eig_max in (0, 1]; returns |error factor|."""
    x = 1.0            # current error (acts on e since b=0 => r = -A x)
    lam = lam_ratio    # A/eig_max combined: r = (b - A x)/eig_max = -lam*x
    r = -lam * x
    z = 2.0 / (1.0 + a) * r
    ro = (1.0 - a) / (1.0 + a)
    for k in range(1, max_iter + 1):
        x = x + z
        if k < max_iter:
            ro_new = 1.0 / (2.0 * (1.0 + a) / (1.0 - a) - ro)
            y = lam * z                       # (A z)/eig_max
            r = r - y
            alpha = ro * ro_new
            beta = 4.0 * ro_new / (1.0 - a)
            z = alpha * z + beta * r
            ro = ro_new
    return x

lam = np.linspace(1e-4, 1.0, 800)
fig, ax = setup(figsize=(9.5, 5.5))

colors = [C_FINE, C_COARSE, C_GC, "#2E8B57"]
for i, m in enumerate([1, 2, 3, 4]):
    damp = np.array([cheb_damp(l, m) for l in lam])
    ax.plot(lam, damp, color=colors[i], lw=2,
            label=f"icheb(1) = {m} iterations")

ax.axhline(0, color="#999999", lw=0.7)
ax.axvspan(0.5, 1.0, color="#F0F0F0", zorder=0)
ax.text(0.75, 0.88, "high-frequency band\n(target of the smoother)",
        ha="center", fontsize=9, color="#555555")
ax.set_xlabel(r"$\lambda\,/\,\lambda_{\max}$   ($\lambda_{\max}$: Gershgorin row-sum bound, eig_value)")
ax.set_ylabel(r"error amplification factor  $p_m(\lambda)$")
ax.set_title("4th-kind Chebyshev polynomial smoother (poly_cheb_smooth, a = 0.3)\n"
             "error factor after m smoothing iterations (code recurrence reproduced numerically)",
             fontsize=11.5)
ax.legend(frameon=False, fontsize=9.5)
ax.set_xlim(0, 1); ax.set_ylim(-0.35, 1.0)
ax.grid(alpha=0.25)
save(fig, "fig03_chebyshev.png")
