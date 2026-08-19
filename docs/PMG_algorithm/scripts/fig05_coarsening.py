"""Figure 05: the distance-based MIS coarsening of coarsening_semi.f90 and the
distance-weighted interpolation of interpolation.f90 (P_distance), reproduced
on a small perturbed 2D point cloud.

The algorithm below is a line-by-line miniature of coarsening_semi.f90:
  - sweep nodes in order; an unclassified node becomes C
  - dmin = min distance to graph neighbors;  dmin <- dmin / teta
  - unclassified neighbors within dmin become F
  - afterwards, F nodes with no C neighbor are promoted to C
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.spatial import Delaunay
from style_common import setup, save, C_FINE, C_COARSE, C_GC

rng = np.random.default_rng(3)
nx, ny = 9, 7
pts = np.array([[i + 0.30 * rng.uniform(-1, 1), j + 0.30 * rng.uniform(-1, 1)]
                for j in range(ny) for i in range(nx)])
n = len(pts)
tri = Delaunay(pts)
nbr = [set() for _ in range(n)]
for s in tri.simplices:
    for a in s:
        for b in s:
            if a != b:
                nbr[a].add(b)

def mis_coarsen(teta):
    state = np.full(n, 2)          # 2 = unclassified, 1 = C, 0 = F
    for j in range(n):
        if state[j] != 2:
            continue
        state[j] = 1               # C point
        d = {k: np.linalg.norm(pts[k] - pts[j]) for k in nbr[j]}
        dmin = min(d.values()) / teta          # radius scaling by teta
        for k, dk in d.items():
            if state[k] == 2 and dk <= dmin:
                state[k] = 0       # F point
    # promotion of isolated F points
    for j in range(n):
        if state[j] == 0 and not any(state[k] == 1 for k in nbr[j]):
            state[j] = 1
    return state

fig, axes = plt.subplots(1, 3, figsize=(15, 5.2))
plt.rcParams.update({"font.size": 11})

for ax, teta in zip(axes[:2], (0.6, 0.9)):
    st = mis_coarsen(teta)
    # edges
    for a in range(n):
        for b in nbr[a]:
            if a < b:
                ax.plot(*zip(pts[a], pts[b]), color="#DDDDDD", lw=0.7, zorder=1)
    F = st == 0; C = st == 1
    ax.scatter(*pts[F].T, s=28, c="#B0B7C0", zorder=3, label="F point (fine only)")
    ax.scatter(*pts[C].T, s=95, c=C_COARSE, ec="k", zorder=4, marker="s",
               label="C point (kept on coarse level)")
    ax.set_title(f"teta = {teta}  ->  {C.sum()}/{n} C points\n"
                 f"(radius = nearest-neighbor dist / teta)", fontsize=10.5)
    ax.set_aspect("equal"); ax.axis("off")
    ax.legend(loc="upper right", fontsize=8, frameon=True)

# --- interpolation stencil panel
ax = axes[2]
st = mis_coarsen(0.6)
C_idx = np.where(st == 1)[0]
# pick an F node with >= 4 C candidates in 1-ring for illustration
center = pts.mean(axis=0)
cands = [(j, [k for k in nbr[j] if st[k] == 1]) for j in np.where(st == 0)[0]]
cands = [(j, cn) for j, cn in cands if 3 <= len(cn) <= 4]
fpick, cnb = min(cands, key=lambda t: np.linalg.norm(pts[t[0]] - center))
for a in range(n):
    for b in nbr[a]:
        if a < b:
            ax.plot(*zip(pts[a], pts[b]), color="#EEEEEE", lw=0.7, zorder=1)
ax.scatter(*pts[st == 0].T, s=20, c="#D5D9DE", zorder=2)
ax.scatter(*pts[st == 1].T, s=80, c=C_COARSE, ec="k", marker="s", zorder=3)
ax.scatter(*pts[fpick], s=130, c=C_FINE, ec="k", zorder=5)
ax.annotate("F node i", pts[fpick], xytext=(pts[fpick][0] - 2.4, pts[fpick][1] + 0.8),
            fontsize=10, color=C_FINE, arrowprops=dict(arrowstyle="->", color=C_FINE))
# inverse-squared-distance weights, as P_distance (dx = squared distance)
d2 = np.array([np.sum((pts[k] - pts[fpick]) ** 2) for k in cnb])
w = (1.0 / d2) / np.sum(1.0 / d2)
for k, wk in zip(cnb, w):
    ax.plot(*zip(pts[fpick], pts[k]), color=C_GC, lw=2.2, zorder=4)
    mid = 0.55 * pts[k] + 0.45 * pts[fpick]
    ax.text(*mid, f"{wk:.2f}", fontsize=9.5, color=C_GC, zorder=6,
            bbox=dict(fc="white", ec=C_GC, lw=0.6, pad=1.5))
ax.set_title("P_distance weights (ip_inter=1):\n"
             r"$w_c = (1/d_c^2)\,/\,\sum_j 1/d_j^2$   (row sum = 1)"
             "\ncandidates: C points in 1-ring (ip_lev=1), max ip_nmax=4",
             fontsize=10)
ax.set_aspect("equal"); ax.axis("off")

fig.subplots_adjust(top=0.72)
fig.suptitle("coarsening_semi.f90: distance-based MIS C/F splitting (no matrix values used)  "
             "+  interpolation.f90 weights", fontsize=12.5, y=0.97)
save(fig, "fig05_coarsening.png")
