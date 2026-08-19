"""Generate all demo-based figures (fig07..fig11) from the 2D miniature PMG."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from style_common import save, RANK_COLORS, RANK_LIGHT, C_GC, C_SEND, C_RECV
from pmg2d import PMG2D, halo_exchange, Traffic

NX = 8
mg = PMG2D(nx=NX, ny=NX, nlevel=3)          # small: clear node-level pictures
mgL = PMG2D(nx=16, ny=16, nlevel=3)         # larger: convergence history

# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def draw_grid(ax, level, owners, size=90, marker="s"):
    xy = mg.coords[level]
    for rk in range(4):
        m = owners == rk
        ax.scatter(xy[m, 0], xy[m, 1], s=size, marker=marker,
                   c=RANK_LIGHT[rk], ec="#888", lw=0.5, zorder=2)
    ax.set_aspect("equal"); ax.axis("off")

# ======================================================================
# fig07: partition + interior/interface/ghost classification + numbering
# ======================================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 6.4))

ax = axes[0]
owners0 = mg.owners[0]
xy = mg.coords[0]
for rk in range(4):
    m = owners0 == rk
    ax.scatter(xy[m, 0], xy[m, 1], s=210, marker="s", c=RANK_LIGHT[rk],
               ec="#777", lw=0.6, zorder=2)
ax.plot([NX/2, NX/2], [0, NX], color="k", lw=1.4, zorder=3)
ax.plot([0, NX], [NX/2, NX/2], color="k", lw=1.4, zorder=3)
for rk, (cx, cy) in enumerate([(NX/4, NX/4), (3*NX/4, NX/4),
                               (NX/4, 3*NX/4), (3*NX/4, 3*NX/4)]):
    ax.text(cx, cy, f"rank {rk}", fontsize=15, ha="center", va="center",
            color=RANK_COLORS[rk], weight="bold", zorder=4)
ax.set_title(f"domain decomposition: {NX}x{NX} cells, 4 ranks\n"
             "(real code: cell ownership celem() from host partitioner,\n"
             "coarse levels inherit fine ownership)", fontsize=10.5)
ax.set_aspect("equal"); ax.axis("off")

ax = axes[1]
rk = 0
rl = mg.levels[0][rk]
own = set(rl.owned)
# interface = owned nodes sent to at least one neighbor (union over A pattern)
send_union = set()
for nb, gids in rl.commA.sintf.items():
    if gids is not None:
        send_union |= set(gids)
interior = [g for g in rl.owned if g not in send_union]
interface = [g for g in rl.owned if g in send_union]
ghosts = rl.ghost
ax.scatter(xy[interior, 0], xy[interior, 1], s=210, marker="s",
           c=RANK_COLORS[0], ec="k", lw=0.5, label=f"interior (1..nintr) : {len(interior)}")
ax.scatter(xy[interface, 0], xy[interface, 1], s=210, marker="s",
           c="#8FB7E8", ec="k", lw=0.8, label=f"interface (nintr+1..nintf) : {len(interface)}")
ax.scatter(xy[ghosts, 0], xy[ghosts, 1], s=190, marker="s",
           c="#CCCCCC", ec="#777", lw=0.8, label=f"ghost (nintf+1..nnode) : {len(ghosts)}")
# local numbering annotation: number a few nodes in local order
order = list(interior) + list(interface) + list(ghosts)
for loc, g in enumerate(order):
    ax.text(xy[g, 0], xy[g, 1], str(loc + 1), fontsize=6.5,
            ha="center", va="center", color="white" if g in own else "#444")
others = [g for g in range(len(owners0)) if owners0[g] != rk and g not in set(ghosts)]
ax.scatter(xy[others, 0], xy[others, 1], s=40, marker="s", c="#F2F2F2",
           ec="#DDD", lw=0.4)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.02), fontsize=9, ncol=2, frameon=False)
ax.set_title("rank 0 local numbering (Domain_infor_FVM_fine ordering):\n"
             "interior first, then interface(owned/sent), then ghosts\n"
             "ghosts contiguous per neighbor rank (nbdom order)", fontsize=10.5)
ax.set_aspect("equal"); ax.axis("off")
save(fig, "fig07_partition_ordering.png")

# ======================================================================
# fig08: halo exchange steps (pack -> ISEND/IRECV -> unpack)
# ======================================================================
tr = Traffic()
vec = [np.zeros(mg.A[0].shape[0]) for _ in range(4)]
for r_ in range(4):
    vec[r_][mg.levels[0][r_].owned] = r_ + 1        # rank-identifying values
halo_exchange(vec, mg.levels[0], "commA", 0, tr)

fig, axes = plt.subplots(1, 3, figsize=(16, 5.6))

# (a) pack: rank0's send lists per neighbor
ax = axes[0]
draw_grid(ax, 0, owners0, size=150)
rl = mg.levels[0][0]
for i, nb in enumerate(sorted(rl.commA.sintf)):
    gids = rl.commA.sintf[nb]
    if gids is None:
        continue
    ax.scatter(xy[gids, 0], xy[gids, 1], s=170, marker="s",
               fc="none", ec=RANK_COLORS[nb], lw=2.4, zorder=5)
    for k, g in enumerate(gids):
        ax.text(xy[g, 0], xy[g, 1], str(k + 1), fontsize=7, ha="center",
                va="center", zorder=6)
ax.set_title("(a) PACK on rank 0:  svar(i) = u(sintf(i))\n"
             "one contiguous buffer segment per neighbor\n"
             "(outline color = destination rank, number = buffer position)",
             fontsize=10)

# (b) messages
ax = axes[1]
draw_grid(ax, 0, owners0, size=150)
centers = [(NX/4, NX/4), (3*NX/4, NX/4), (NX/4, 3*NX/4), (3*NX/4, 3*NX/4)]
for m in tr.log:
    if m["pattern"] != "commA":
        continue
    s, d = m["src"], m["dst"]
    p1 = np.array(centers[s]); p2 = np.array(centers[d])
    p1s = p1 + 0.12 * (p2 - p1); p2s = p2 - 0.12 * (p2 - p1)
    ax.add_patch(FancyArrowPatch(p1s, p2s, arrowstyle="-|>", mutation_scale=15,
                                 color=RANK_COLORS[s], lw=1.8, zorder=6,
                                 connectionstyle="arc3,rad=0.18"))
    mid = 0.5 * (p1s + p2s) + 0.32 * np.array([p2s[1]-p1s[1], p1s[0]-p2s[0]]) / np.linalg.norm(p2s-p1s)
    ax.text(*mid, f"{len(m['gids'])}", fontsize=9, color=RANK_COLORS[s],
            ha="center", va="center", zorder=7,
            bbox=dict(fc="white", ec=RANK_COLORS[s], lw=0.6, pad=1.2))
ax.set_title("(b) MPI_ISEND / MPI_IRECV per neighbor pair\n"
             "(numbers = message length; all posted non-blocking,\n"
             "then MPI_WAIT on all - no compute/comm overlap)", fontsize=10)

# (c) unpack: rank0 ghosts filled with senders' values
ax = axes[2]
draw_grid(ax, 0, owners0, size=150)
gh = mg.levels[0][0].commA
for nb in gh.nbdom:
    gids = gh.rintf.get(nb)
    if gids is None:
        continue
    ax.scatter(xy[gids, 0], xy[gids, 1], s=200, marker="s",
               c=RANK_COLORS[nb], ec="k", lw=1.0, zorder=5)
    for g in gids:
        ax.text(xy[g, 0], xy[g, 1], f"{vec[0][g]:.0f}", fontsize=8,
                ha="center", va="center", color="white", zorder=6)
ax.set_title("(c) UNPACK on rank 0:  u(rintf(i)) = rvar(i)\n"
             "ghost values overwritten (pure copy - no summation);\n"
             "ghost blocks are contiguous, so rvar maps 1:1", fontsize=10)

fig.suptitle("send_receive.f90: one halo exchange of the fine-level A-pattern (2D miniature, values = owner rank+1)",
             fontsize=12.5, y=1.03)
save(fig, "fig08_halo_steps.png")

# ======================================================================
# fig09: the three communication patterns A / R / P
# ======================================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5.8))
rk = 0
lvl = 1     # first coarse level
own_c = mg.levels[lvl][rk].owned
own_f = mg.levels[0][rk].owned
xyc = mg.coords[lvl]; xyf = mg.coords[0]
oc = mg.owners[lvl]; of = mg.owners[0]

# --- A pattern on coarse level
ax = axes[0]
draw_grid(ax, lvl, oc, size=300)
ax.scatter(xyc[own_c, 0], xyc[own_c, 1], s=300, marker="s", c=RANK_COLORS[0],
           ec="k", zorder=4)
gA = [g for nb in mg.levels[lvl][rk].commA.nbdom
      for g in (mg.levels[lvl][rk].commA.rintf.get(nb) if mg.levels[lvl][rk].commA.rintf.get(nb) is not None else [])]
ax.scatter(xyc[gA, 0], xyc[gA, 1], s=330, marker="s", fc="none",
           ec=C_GC, lw=2.6, zorder=5)
ax.set_title("(A)  smoothing / residual on level l\n"
             "ghosts = off-rank columns of my rows of $A_l$\n"
             "used by Smooth_GS2, resi;  MD_S_R_NEW(id=1)", fontsize=10)

# --- R pattern: fine ghosts needed by my owned coarse rows of R
ax = axes[1]
draw_grid(ax, 0, of, size=110)
ax.scatter(xyf[own_f, 0], xyf[own_f, 1], s=110, marker="s", c=RANK_LIGHT[0],
           ec="#777", zorder=3)
ax.scatter(xyc[own_c, 0], xyc[own_c, 1], s=340, marker="s", fc="none",
           ec=RANK_COLORS[0], lw=2.2, zorder=5)
pat = mg.levels[0][rk].commR
gR = [g for nb in pat.nbdom for g in (pat.rintf.get(nb) if pat.rintf.get(nb) is not None else [])]
ax.scatter(xyf[gR, 0], xyf[gR, 1], s=140, marker="s", c=C_GC, ec="k", zorder=4)
ax.set_title("(R)  restriction  $r_c = R\\,r_f$\n"
             "red = off-rank FINE nodes referenced by R-rows of my\n"
             "owned coarse points (blue squares);  MD_S_R_NEW(id=2)", fontsize=10)

# --- P pattern: coarse ghosts needed by P rows of my owned fine nodes
ax = axes[2]
draw_grid(ax, 0, of, size=110)
ax.scatter(xyf[own_f, 0], xyf[own_f, 1], s=110, marker="s", c=RANK_LIGHT[0],
           ec="#777", zorder=3)
pat = mg.levels[lvl][rk].commP
gP = [g for nb in pat.nbdom for g in (pat.rintf.get(nb) if pat.rintf.get(nb) is not None else [])]
ax.scatter(xyc[oc == 0, 0][:], xyc[oc == 0, 1][:], s=340, marker="s",
           fc="none", ec=RANK_COLORS[0], lw=2.2, zorder=5)
ax.scatter(xyc[gP, 0], xyc[gP, 1], s=360, marker="D", c=C_RECV, ec="k", zorder=6)
ax.set_title("(P)  prolongation  $e_f = P\\,e_c$\n"
             "purple = off-rank COARSE points referenced by P-rows of\n"
             "my owned fine nodes;  MD_S_R_NEW(id=3)", fontsize=10)

fig.suptitle("Three per-level communication patterns (Neighbor_node_ARP): base halo list = A ∪ R ∪ P,\n"
             "each operator exchanges only its own minimal subset", fontsize=12.5, y=1.06)
save(fig, "fig09_arp_patterns.png")

# ======================================================================
# fig10: coarsest-level ALLGATHERV + redundant solve
# ======================================================================
fig, ax = plt.subplots(figsize=(12, 6.4))
ax.set_xlim(0, 12); ax.set_ylim(0, 10); ax.axis("off")
lc = 2
nG = mgL.A[lc].shape[0]
counts = [len(mgL.levels[lc][r_].owned) for r_ in range(4)]
W = 7.2; x0 = 0.6
# per-rank local coarsest vectors
y0 = 8.4
ax.text(x0 + W/2, 9.55, "local coarsest residual r$_s$ on each rank (owned part, length nintfs)",
        ha="center", fontsize=10)
xx = x0
for r_ in range(4):
    w = W * counts[r_] / nG
    ax.add_patch(Rectangle((xx, y0), w, 0.6, fc=RANK_COLORS[r_], ec="k"))
    ax.text(xx + w/2, y0 + 0.3, f"rank {r_}\n({counts[r_]})", ha="center",
            va="center", fontsize=8, color="white")
    xx += w + 0.35
# gathered rG0
y1 = 6.2
xx2 = x0 + 0.5
ax.text(x0 + W/2 + 0.7, y1 + 0.95, "MPI_ALLGATHERV(r, nsengatR, ..., rG0, irevgatR, idispR, ...)",
        ha="center", fontsize=10, color=C_GC)
xx = xx2
for r_ in range(4):
    w = W * counts[r_] / nG
    ax.add_patch(Rectangle((xx, y1), w, 0.6, fc=RANK_COLORS[r_], ec="k"))
    ax.add_patch(FancyArrowPatch((x0 + sum(W*counts[q]/nG + 0.35 for q in range(r_)) + W*counts[r_]/nG/2, y0),
                                 (xx + w/2, y1 + 0.62), arrowstyle="-|>",
                                 mutation_scale=11, color=RANK_COLORS[r_], lw=1.2))
    xx += w
ax.text(xx + 0.15, y1 + 0.3, "rG0 (rank-ordered\nconcatenation)", fontsize=8.5, va="center")
# permutation to rG
y2 = 4.4
ax.text(x0 + W/2 + 0.7, y2 + 0.95, "rG(imapgatR(j)) = rG0(j)   (permutation to global coarse numbering)",
        ha="center", fontsize=10)
ax.add_patch(Rectangle((xx2, y2), W, 0.6, fc="#7A8894", ec="k"))
ax.text(xx2 + W/2, y2 + 0.3, f"rG  (global numbering, length nnodeG = {nG})",
        ha="center", va="center", fontsize=9, color="white")
for f_ in (0.1, 0.35, 0.6, 0.85):
    ax.add_patch(FancyArrowPatch((xx2 + W*f_, y1), (xx2 + W*((f_*7+2)%10)/10*1.0, y2 + 0.62),
                                 arrowstyle="-|>", mutation_scale=9, color="#666",
                                 lw=0.9, linestyle="--"))
# redundant solve on every rank
y3 = 1.6
for r_ in range(4):
    bx = 0.9 + r_ * 2.75
    ax.add_patch(Rectangle((bx, y3), 2.4, 1.5, fc=RANK_LIGHT[r_], ec=RANK_COLORS[r_], lw=1.6))
    ax.text(bx + 1.2, y3 + 0.75,
            f"rank {r_}\neG = A$_{{GC}}^{{-1}}$ rG\n(dense inverse Ainv,\nidentical on all ranks)",
            ha="center", va="center", fontsize=8.4)
    ax.add_patch(FancyArrowPatch((xx2 + W*(0.15 + 0.23*r_), y2), (bx + 1.2, y3 + 1.52),
                                 arrowstyle="-|>", mutation_scale=11, color=RANK_COLORS[r_], lw=1.2))
ax.text(6.0, 0.75, "no broadcast needed afterwards: every rank already holds the full eG;\n"
                   "local copy-back  e(i) = eG(imapG(i))  also fills ghosts consistently -> no halo exchange",
        ha="center", fontsize=9.5, color="#333")
ax.set_title("Coarsest 'Global Coarse' level (SOLVE_GC_all, igather=1):\n"
             "gather-to-all + redundant serial solve replaces all coarse-level halo traffic",
             fontsize=12)
save(fig, "fig10_coarsest_gather.png")

# ======================================================================
# fig11: convergence history PMG vs plain BiCGSTAB (+ V-cycle traffic table)
# ======================================================================
rng = np.random.default_rng(0)
b = rng.standard_normal(mgL.A[0].shape[0])
x1, h1 = mgL.bicgstab(b, precond=True)
x2, h2 = mgL.bicgstab(b, precond=False)

# count messages of one recorded V-cycle on the small case
mg.traffic = Traffic()
bs = [b[:mg.A[0].shape[0]].copy() * 0 + 1 for _ in range(4)]
xs = [np.zeros(mg.A[0].shape[0]) for _ in range(4)]
mg.vcycle(bs, xs, record=True)
from collections import Counter
cnt = Counter((m["level"], m["pattern"]) for m in mg.traffic.log)
vol = Counter()
for m in mg.traffic.log:
    vol[(m["level"], m["pattern"])] += len(m["gids"])

fig, axes = plt.subplots(1, 2, figsize=(13, 5.4))
ax = axes[0]
ax.semilogy(np.arange(len(h2)), h2 / h2[0], "-o", ms=3.5, color="#8A94A2",
            label=f"plain BiCGSTAB ({len(h2)-1} iters)")
ax.semilogy(np.arange(len(h1)), h1 / h1[0], "-o", ms=3.5, color=RANK_COLORS[0],
            label=f"PMG-preconditioned ({len(h1)-1} iters)")
ax.axhline(1e-8, color="#BBB", lw=0.8, ls="--")
ax.set_xlabel("BiCGSTAB iteration"); ax.set_ylabel(r"$\|r\|/\|r_0\|$")
ax.set_title("2D miniature, 16x16 Poisson, 4 ranks, 3 levels", fontsize=10.5)
ax.legend(frameon=False); ax.grid(alpha=0.25)

ax = axes[1]
ax.axis("off")
rows = []
name = {"commA": "A (smooth)", "commR": "R (restrict)",
        "commP": "P (prolong)", "GATHER": "ALLGATHERV"}
for (lvl_, pat), c in sorted(cnt.items()):
    rows.append([f"level {lvl_+1}", name[pat], str(c), str(vol[(lvl_, pat)])])
tbl = ax.table(cellText=rows,
               colLabels=["level", "pattern", "messages", "values moved"],
               loc="center", cellLoc="center")
tbl.auto_set_font_size(False); tbl.set_fontsize(9.5); tbl.scale(1, 1.5)
ax.set_title("message count per ONE V-cycle (8x8 miniature, 4 ranks)\n"
             "recorded by the demo's traffic logger", fontsize=10.5)
save(fig, "fig11_convergence_traffic.png")
print("demo figures done")
