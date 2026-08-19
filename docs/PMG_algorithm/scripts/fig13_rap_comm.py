"""Figure 13: communication timeline of the parallel RAP (owner-computes rows)
vs. the alternative additive (partial-sum + reduce) scheme."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from style_common import setup, save, RANK_COLORS, RANK_LIGHT, C_GC, C_SEND

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(14, 7.2), width_ratios=[1.35, 1])
for a in (ax, ax2):
    a.axis("off")

# ================= left: actual scheme timeline =================
ax.set_xlim(0, 10); ax.set_ylim(0, 12)
X0, X1, W = 0.35, 5.35, 4.3

def box(a, x, y, w, h, text, fc, tc="white", fs=9):
    a.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06",
                               fc=fc, ec="#333", lw=0.7))
    a.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, color=tc)

ax.text(X0 + W/2, 11.55, "rank 0", fontsize=12, ha="center", color=RANK_COLORS[0], weight="bold")
ax.text(X1 + W/2, 11.55, "rank 1", fontsize=12, ha="center", color=RANK_COLORS[1], weight="bold")
for x in (X0 + W/2, X1 + W/2):
    ax.plot([x, x], [0.6, 11.3], color="#DDDDDD", lw=1.0, zorder=0)

# setup once
box(ax, X0, 10.2, W, 0.85, "SETUP (once)\nsend_receive_csr: exchange ghost-row\ncolumn coordinates -> reorder ja", "#7A8894", fs=8.2)
box(ax, X1, 10.2, W, 0.85, "SETUP (once)\nsame: column ORDER of every ghost row\nnow matches the owner's storage", "#7A8894", fs=8.2)

# step 1: row-value exchange
box(ax, X0, 8.3, W, 1.0, "1) import ghost-row VALUES of A\nsend_receive_mtf / MD_S_R_MT(ilv)\nrecv A(4,:) | send A(3,:)", C_SEND, fs=8.6)
box(ax, X1, 8.3, W, 1.0, "1) import ghost-row VALUES of A\nrecv A(3,:) | send A(4,:)", C_SEND, fs=8.6)
ax.add_patch(FancyArrowPatch((X1, 8.95), (X0 + W + 0.05, 9.05), arrowstyle="-|>",
             mutation_scale=13, color=RANK_COLORS[1], lw=1.6))
ax.text(5.0, 9.42, "A(4,:) values", fontsize=8, color=RANK_COLORS[1], ha="center")
ax.add_patch(FancyArrowPatch((X0 + W + 0.05, 8.55), (X1, 8.45), arrowstyle="-|>",
             mutation_scale=13, color=RANK_COLORS[0], lw=1.6))
ax.text(5.0, 8.02, "A(3,:) values", fontsize=8, color=RANK_COLORS[0], ha="center")

# step 2: local SpGEMM
box(ax, X0, 6.1, W, 1.5,
    "2) local SpGEMM - OWNED coarse rows only\nfor I in {C1, C2}:\n  vi = SUM_k R(I,k) A(k,:)\n  Ac(I,J) = SUM_l vi(l) R(J,l)\n(OpenMP over I, vi firstprivate)",
    RANK_COLORS[0], fs=8.2)
box(ax, X1, 6.1, W, 1.5,
    "2) local SpGEMM - OWNED coarse rows only\nfor I in {C3}:\n  vi = SUM_k R(I,k) A(k,:)\n  Ac(I,J) = SUM_l vi(l) R(J,l)",
    RANK_COLORS[1], fs=8.2)

# step 3: no reduce
box(ax, X0, 4.6, W + (X1 - X0), 0.85,
    "3) NO communication after RAP: every owned Ac row is already complete\n"
    "(all fine rows k with R(I,k) != 0 are local after step 1)", "#3d6b47", fs=9)

# step 4: recursion
box(ax, X0, 2.9, W + (X1 - X0), 1.1,
    "4) recursion to level ilv+1: ghost COARSE rows of Ac are imported by\n"
    "MD_S_R_MT(ilv+1) right before the next RAP - same pattern, one level down\n"
    "(this replaces any Ac assembly communication)", C_GC, fs=9)

ax.text(5.0, 1.9, "per V-cycle cost: 0 (RAP runs only when the matrix changes, icase=1)\n"
                  "per RAP cost: ONE row-value halo exchange per level",
        fontsize=9.5, ha="center", color="#333")
ax.set_title("Actual scheme (stiffness_MG): owner-computes rows,\n'move A-rows to the computation'", fontsize=11.5)

# ================= right: rejected additive alternative =================
ax2.set_xlim(0, 10); ax2.set_ylim(0, 12)
box(ax2, 1.0, 9.2, 8.0, 1.3,
    "each rank k-splits the sum by FINE-row ownership:\n"
    "Ac_partial(I,J) = SUM_{k owned} R(I,k) [A(k,:) P](J)\n"
    "-> contributions to NON-owned coarse rows appear", "#7A8894", fs=8.8)
box(ax2, 1.0, 7.0, 8.0, 1.2,
    "additive exchange required:\nAc(I,J) = SUM over ranks of partial sums\n"
    "(FEM-assembly-style reduce on shared coarse rows)", "#a05252", fs=8.8)
box(ax2, 1.0, 5.2, 8.0, 1.2,
    "needs a second, ADDITIVE halo primitive\n(sum into ghosts, then return) -\n"
    "does not exist anywhere in this codebase", "#a05252", fs=8.8)
for y1, y2 in [(9.2, 8.2), (7.0, 6.4)]:
    ax2.add_patch(FancyArrowPatch((5.0, y1), (5.0, y2 + 0.0), arrowstyle="-|>",
                  mutation_scale=13, color="#333", lw=1.2))
ax2.text(5.0, 2.35, "trade-off taken by the code:\n"
                   "import ghost A-ROWS (values, one exchange)\n"
                   "instead of reducing partial Ac entries.\n\n"
                   "-> keeps the single 'overwrite-only' halo model\n"
                   "used everywhere else (vectors and matrices),\n"
                   "at the price of the 2nd halo layer (nnodegl)\n"
                   "for ghost-row column indexing.",
         fontsize=9.5, ha="center", color="#333",
         bbox=dict(fc="#F7F8FA", ec="#BBB", lw=0.8, pad=6))
ax2.set_title("Rejected alternative: partial sums + reduce\n('move partial results to the owner')", fontsize=11.5)

save(fig, "fig13_rap_comm.png")
