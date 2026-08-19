"""2D miniature of the CUPID PMG-preconditioned BiCGSTAB solver.

This is a faithful *structural* miniature of code/Source/GMG:
  - domain decomposition with owned-first / ghost-last local numbering
    (owned 0..nintf-1  <->  Fortran 1..nintf;  ghosts nintf..nnode-1),
  - halo-exchange lists per neighbor rank (nbdom, spt/rpt, sintf/rintf)
    with explicit pack -> "ISEND/IRECV" -> unpack staging, mirroring
    send_receive.f90,
  - three distinct communication patterns per level (A / R / P) derived
    from the sparsity of the level matrix, the restriction and the
    prolongation, mirroring Neighbor_node_ARP.f90 / MD_S_R_NEW,
  - a V-cycle with 4th-kind Chebyshev smoothing on the finest level and
    Gauss-Seidel on coarse levels (SOLVER_NEW in 7_SOLVE_GMG.f90),
  - a gathered, redundantly-solved global-coarse level (SOLVE_GC.f90:
    MPI_ALLGATHERV + dense solve on every rank),
  - a BiCGSTAB outer loop applying the V-cycle twice per iteration
    (6_solver_pbcg_mg.f90).

Ranks are simulated in a single process; all "MPI" traffic goes through
explicit send/recv buffers so each communication step can be visualized.
Grids are structured 2D only to keep the miniature readable - the real
code is unstructured and matrix-driven (strength-of-connection coarsening),
which is documented in the accompanying notes.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field

# ----------------------------------------------------------------------
# global problem: 5-point Poisson on an nx x ny cell grid (Dirichlet)
# ----------------------------------------------------------------------

def poisson2d(nx, ny):
    n = nx * ny
    A = np.zeros((n, n))
    def gid(i, j): return j * nx + i
    for j in range(ny):
        for i in range(nx):
            g = gid(i, j)
            A[g, g] = 4.0
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < nx and 0 <= jj < ny:
                    A[g, gid(ii, jj)] = -1.0
    return A

# ----------------------------------------------------------------------
# geometric hierarchy (miniature of coarsening + interpolation setup)
# ----------------------------------------------------------------------

def build_hierarchy(nx, ny, nlevel):
    """Returns per-level: grid dims, coordinates, A (Galerkin), P (lev->lev-1).

    P uses inverse-distance weights to the ip_nmax=4 nearest coarse points
    (ip_inter=1 analogue); R = row-normalized P^T (the real code builds
    Xrest separately, but with the same neighbor/weight logic).
    """
    dims = [(nx, ny)]
    coords = [np.array([[i + 0.5, j + 0.5] for j in range(ny) for i in range(nx)])]
    A = [poisson2d(nx, ny)]
    P, R = [], []
    for l in range(1, nlevel):
        fx, fy = dims[-1]
        cx, cy = (fx + 1) // 2, (fy + 1) // 2
        cc = np.array([[2.0 * i + 1.0, 2.0 * j + 1.0] for j in range(cy) for i in range(cx)])
        cc *= 2 ** (l - 1)
        fc = coords[-1]
        nf, nc = len(fc), len(cc)
        Pl = np.zeros((nf, nc))
        for f in range(nf):
            d = np.linalg.norm(cc - fc[f], axis=1)
            near = np.argsort(d)[:4]                      # ip_nmax = 4
            if d[near[0]] < 1e-12:
                Pl[f, near[0]] = 1.0
                continue
            w = 1.0 / d[near]                             # distance weights
            Pl[f, near] = w / w.sum()
        Rl = Pl.T.copy()
        Rl /= Rl.sum(axis=1, keepdims=True)               # row-normalized
        Ac = Rl @ A[-1] @ Pl                              # Galerkin RAP
        dims.append((cx, cy)); coords.append(cc)
        A.append(Ac); P.append(Pl); R.append(Rl)
    return dims, coords, A, P, R

# ----------------------------------------------------------------------
# domain decomposition: 2x2 ranks, owned-first ghost-last numbering
# ----------------------------------------------------------------------

@dataclass
class CommPattern:
    """One halo-exchange pattern (the analogue of one (spt,rpt,sintf,rintf,nbdom) set)."""
    nbdom: list = field(default_factory=list)   # neighbor rank ids
    sintf: dict = field(default_factory=dict)   # nb -> local indices to SEND (owned)
    rintf: dict = field(default_factory=dict)   # nb -> local indices to RECV (ghost)

@dataclass
class RankLevel:
    owned: np.ndarray          # global ids of owned nodes (local 0..nintf-1)
    ghost: np.ndarray          # global ids of ghost nodes (local nintf..nnode-1)
    g2l: dict                  # global -> local
    commA: CommPattern
    commR: CommPattern = None  # ghosts of the *finer* level needed by R rows
    commP: CommPattern = None  # ghosts of *this* level needed by P rows of finer

def owner_of(coords, dims, npx=2, npy=2, extent=None):
    """Block ownership by coordinate (levels inherit the same spatial cut)."""
    ex, ey = extent
    px = np.minimum((coords[:, 0] / ex * npx).astype(int), npx - 1)
    py = np.minimum((coords[:, 1] / ey * npy).astype(int), npy - 1)
    return py * npx + px

def ghost_set(M_rows_owned, owned_mask):
    """Column dependencies outside the owned set = required ghosts."""
    cols = np.where(np.abs(M_rows_owned).sum(axis=0) > 1e-14)[0]
    return np.array([c for c in cols if not owned_mask[c]], dtype=int)

def build_pattern(rank, owners, need_ghost_gids, owned_sorted):
    """Make send/recv lists; receiver-side ordering fixes the buffer order,
    the matching owner sends in that same order (as the Fortran setup does)."""
    pat = CommPattern()
    for nb in range(owners.max() + 1):
        if nb == rank:
            continue
        gids = [g for g in need_ghost_gids if owners[g] == nb]
        if gids:
            pat.nbdom.append(nb)
            pat.rintf[nb] = np.array(gids, dtype=int)      # global ids (recv side)
    return pat

def decompose(dims, coords, A, P, R, npr=(2, 2)):
    """Per level, per rank: owned/ghost sets and A/R/P communication patterns."""
    nlevel = len(A)
    nranks = npr[0] * npr[1]
    extent = (dims[0][0], dims[0][1])
    levels = []           # levels[l][rank] = RankLevel
    owners_all = []
    for l in range(nlevel):
        owners = owner_of(coords[l], dims[l], npr[0], npr[1], extent)
        owners_all.append(owners)
    for l in range(nlevel):
        owners = owners_all[l]
        per_rank = []
        for rk in range(nranks):
            owned = np.where(owners == rk)[0]
            omask = np.zeros(len(owners), bool); omask[owned] = True
            # --- A pattern: stencil of my owned rows of A_l
            gA = ghost_set(A[l][owned, :], omask)
            # --- R pattern (l < nlevel-1): my owned *coarse* rows of R need fine ghosts
            # --- P pattern (l >= 1): rows of P (fine, owned by me on level l-1) need
            #     coarse ghosts of level l
            ghost = list(dict.fromkeys(list(gA)))
            gR = gP = np.array([], int)
            if l < nlevel - 1:
                co = np.where(owners_all[l + 1] == rk)[0]
                fmask = omask
                gR_all = ghost_set(R[l][co, :], fmask)     # fine-level ghosts for restriction
                gR = gR_all
                ghost += [g for g in gR if g not in ghost]
            if l >= 1:
                fo = np.where(owners_all[l - 1] == rk)[0]
                cmask = omask
                gP_all = ghost_set(P[l - 1][fo, :], cmask)  # coarse-level ghosts for prolongation
                gP = gP_all
                ghost += [g for g in gP if g not in ghost]
            ghost = np.array(ghost, int)
            g2l = {g: i for i, g in enumerate(owned)}
            g2l.update({g: len(owned) + i for i, g in enumerate(ghost)})
            rl = RankLevel(owned=owned, ghost=ghost, g2l=g2l,
                           commA=build_pattern(rk, owners, gA, owned),
                           commR=build_pattern(rk, owners, gR, owned),
                           commP=build_pattern(rk, owners, gP, owned))
            per_rank.append(rl)
        # fill matching send lists (owner side) in receiver order
        # (mirrors the setup rule: sender packs in the order fixed by the
        #  receiver's rintf, so no index list travels at runtime)
        for rk in range(nranks):
            for pat_name in ("commA", "commR", "commP"):
                pat = getattr(per_rank[rk], pat_name)
                for nb in pat.nbdom:
                    gids = pat.rintf.get(nb)
                    if gids is None:
                        continue
                    opat = getattr(per_rank[nb], pat_name)
                    opat.sintf[rk] = gids             # owner nb -> rk, receiver order
                    if rk not in opat.nbdom:
                        opat.nbdom.append(rk)
        levels.append(per_rank)
    return levels, owners_all

# ----------------------------------------------------------------------
# halo exchange through explicit buffers (send_receive.f90 miniature)
# ----------------------------------------------------------------------

class Traffic:
    """Records every message for visualization: (level, pattern, src, dst, nvals)."""
    def __init__(self):
        self.log = []
    def add(self, level, pattern, src, dst, gids, vals):
        self.log.append(dict(level=level, pattern=pattern, src=src, dst=dst,
                             gids=np.array(gids), vals=np.array(vals)))

def halo_exchange(vecs, per_rank, pattern, level, traffic=None):
    """vecs[rank] is the GLOBAL-length array view of that rank's data.
    Only owned entries of the sender are trusted; ghost entries of the
    receiver are overwritten - exactly the contract of send_receive."""
    nranks = len(per_rank)
    # pack + send  (svar = u(sintf))
    inbox = {rk: [] for rk in range(nranks)}
    for rk in range(nranks):
        pat = getattr(per_rank[rk], pattern)
        for nb, gids in pat.sintf.items():
            if gids is None:
                continue
            svar = vecs[rk][gids]                  # pack
            inbox[nb].append((rk, gids, svar))     # "MPI_ISEND"
            if traffic is not None:
                traffic.add(level, pattern, rk, nb, gids, svar)
    # recv + unpack  (u(rintf) = rvar)
    for rk in range(nranks):
        for (src, gids, rvar) in inbox[rk]:
            vecs[rk][gids] = rvar                  # unpack into ghosts

# ----------------------------------------------------------------------
# smoothers (poly_smooth.f90 / Relax_GSP.f90 miniatures)
# ----------------------------------------------------------------------

def cheb4_smooth(A, b, x, owned_masks, per_rank, level, eig_max, m=2, a=0.3,
                 traffic=None):
    """4th-kind Chebyshev smoothing, method-2 recurrence of poly_cheb_smooth.
    x is a list of global-length arrays (one per rank)."""
    nranks = len(per_rank)
    halo_exchange(x, per_rank, "commA", level, traffic)
    r, z = [], []
    for rk in range(nranks):
        ow = owned_masks[rk]
        rr = np.zeros_like(x[rk]); rr[ow] = (b[rk][ow] - (A @ x[rk])[ow]) / eig_max
        zz = np.zeros_like(rr);    zz[ow] = 2.0 / (1.0 + a) * rr[ow]
        r.append(rr); z.append(zz)
    ro = (1.0 - a) / (1.0 + a)
    for k in range(1, m + 1):
        for rk in range(nranks):
            ow = owned_masks[rk]
            x[rk][ow] += z[rk][ow]
        if k < m:
            ro_new = 1.0 / (2.0 * (1.0 + a) / (1.0 - a) - ro)
            halo_exchange(z, per_rank, "commA", level, traffic)
            for rk in range(nranks):
                ow = owned_masks[rk]
                y = (A @ z[rk])[ow] / eig_max
                r[rk][ow] -= y
                z[rk][ow] = ro * ro_new * z[rk][ow] + 4.0 * ro_new / (1.0 - a) * r[rk][ow]
            ro = ro_new
    return x

def gs_smooth(A, b, x, owned, diag, relax_ghost=True):
    """Jacobi-across-ranks / Gauss-Seidel-within-rank sweep on owned rows
    (Smooth_GS2: u_i += (b_i - sum a_ij u_j) / a_ii, ghosts frozen)."""
    for i in owned:
        x[i] += (b[i] - A[i] @ x) / diag[i]

# ----------------------------------------------------------------------
# V-cycle (SOLVER_NEW miniature) and BiCGSTAB (solve_pbcg_mg miniature)
# ----------------------------------------------------------------------

class PMG2D:
    def __init__(self, nx=16, ny=16, nlevel=3, npr=(2, 2), itergs=(1, 1, 2)):
        self.dims, self.coords, self.A, self.P, self.R = build_hierarchy(nx, ny, nlevel)
        self.levels, self.owners = decompose(self.dims, self.coords,
                                             self.A, self.P, self.R, npr)
        self.nranks = npr[0] * npr[1]
        self.nlevel = nlevel
        self.itergs = itergs
        self.eig_max = max(np.abs(self.A[0]).sum(axis=1))   # Gershgorin (eig_value)
        self.diag = [np.diag(a).copy() for a in self.A]
        self.traffic = Traffic()

    def masks(self, l):
        return [np.isin(np.arange(self.A[l].shape[0]), self.levels[l][rk].owned)
                for rk in range(self.nranks)]

    def vcycle(self, b_ranks, x_ranks, record=False):
        tr = self.traffic if record else None
        L, A, P, R = self.nlevel, self.A, self.P, self.R
        masks = [self.masks(l) for l in range(L)]
        r = [[np.zeros(A[l].shape[0]) for _ in range(self.nranks)] for l in range(L)]
        e = [[np.zeros(A[l].shape[0]) for _ in range(self.nranks)] for l in range(L)]

        # ---- finest: pre-smooth (Chebyshev) + residual + S&R(R-pattern)
        cheb4_smooth(A[0], b_ranks, x_ranks, masks[0], self.levels[0], 0,
                     self.eig_max, m=self.itergs[0] * 2, traffic=tr)
        halo_exchange(x_ranks, self.levels[0], "commA", 0, tr)
        for rk in range(self.nranks):
            ow = masks[0][rk]
            r[0][rk][ow] = (b_ranks[rk] - A[0] @ x_ranks[rk])[ow]
        halo_exchange(r[0], self.levels[0], "commR", 0, tr)

        # ---- down-leg over middle levels
        for l in range(1, L - 1):
            for rk in range(self.nranks):
                ow = masks[l][rk]
                r[l][rk][ow] = (R[l - 1] @ r[l - 1][rk])[ow]     # rc = R rt
            for _ in range(self.itergs[min(l, len(self.itergs) - 1)]):
                for rk in range(self.nranks):
                    gs_smooth(A[l], r[l][rk], e[l][rk],
                              self.levels[l][rk].owned, self.diag[l])
                halo_exchange(e[l], self.levels[l], "commA", l, tr)   # S&R(A) each sweep
            for rk in range(self.nranks):
                ow = masks[l][rk]
                r[l][rk][ow] = (r[l][rk] - A[l] @ e[l][rk])[ow]      # rt = rc - A e
            halo_exchange(r[l], self.levels[l], "commR", l, tr)

        # ---- coarsest: restrict + ALLGATHERV + redundant direct solve
        lc = L - 1
        for rk in range(self.nranks):
            ow = masks[lc][rk]
            r[lc][rk][ow] = (R[lc - 1] @ r[lc - 1][rk])[ow]
        nG = A[lc].shape[0]
        rG = np.zeros(nG)
        for rk in range(self.nranks):                      # ALLGATHERV analogue
            ow = self.levels[lc][rk].owned
            rG[ow] = r[lc][rk][ow]
            if tr is not None:
                tr.add(lc, "GATHER", rk, -1, ow, r[lc][rk][ow])
        eG = np.linalg.solve(A[lc], rG)                    # SOLVE_EXACT on every rank
        for rk in range(self.nranks):
            e[lc][rk][:] = eG                              # every rank holds full eG

        # ---- up-leg
        for l in range(L - 2, 0, -1):
            for rk in range(self.nranks):
                ow = masks[l][rk]
                e[l][rk][ow] += (P[l] @ e[l + 1][rk])[ow]  # e += P e_coarse
            for _ in range(self.itergs[min(l, len(self.itergs) - 1)]):
                halo_exchange(e[l], self.levels[l], "commA", l, tr)
                for rk in range(self.nranks):
                    gs_smooth(A[l], r[l][rk], e[l][rk],
                              self.levels[l][rk].owned, self.diag[l])
            halo_exchange(e[l], self.levels[l], "commP", l, tr)   # S&R(P) before prolong

        # ---- finest: correct + post-smooth
        for rk in range(self.nranks):
            ow = masks[0][rk]
            x_ranks[rk][ow] += (P[0] @ e[1][rk])[ow]
        halo_exchange(x_ranks, self.levels[0], "commA", 0, tr)
        cheb4_smooth(A[0], b_ranks, x_ranks, masks[0], self.levels[0], 0,
                     self.eig_max, m=self.itergs[0] * 2, traffic=tr)
        return x_ranks

    # ---- outer BiCGSTAB with the V-cycle as M^{-1} (solve_pbcg_mg miniature)
    def bicgstab(self, b, tol=1e-8, maxit=60, precond=True):
        A = self.A[0]; n = A.shape[0]
        masks = self.masks(0)
        def Minv(v):
            if not precond:
                return v.copy()
            bs = [v.copy() for _ in range(self.nranks)]
            xs = [np.zeros(n) for _ in range(self.nranks)]
            self.vcycle(bs, xs)
            out = np.zeros(n)
            for rk in range(self.nranks):
                out[masks[rk]] = xs[rk][masks[rk]]
            return out
        x = np.zeros(n); r = b - A @ x; rb = r.copy()
        rho = alpha = omega = 1.0; p = np.zeros(n); v = np.zeros(n)
        hist = [np.linalg.norm(r)]
        for it in range(maxit):
            rho, rhold = rb @ r, rho
            beta = (rho / rhold) * (alpha / omega) if it else 0.0
            p = r + beta * (p - omega * v) if it else r.copy()
            y = Minv(p); v = A @ y
            alpha = rho / (rb @ v)
            s = r - alpha * v
            z = Minv(s); t = A @ z
            omega = (t @ s) / (t @ t)
            x += alpha * y + omega * z
            r = s - omega * t
            hist.append(np.linalg.norm(r))
            if hist[-1] / hist[0] < tol:
                break
        return x, np.array(hist)


if __name__ == "__main__":
    mg = PMG2D(nx=16, ny=16, nlevel=3)
    rng = np.random.default_rng(0)
    b = rng.standard_normal(mg.A[0].shape[0])
    x, hist = mg.bicgstab(b)
    print("PMG-BiCGSTAB iters:", len(hist) - 1, " final rel.res:", hist[-1] / hist[0])
    x2, hist2 = mg.bicgstab(b, precond=False)
    print("plain BiCGSTAB iters:", len(hist2) - 1, " final rel.res:", hist2[-1] / hist2[0])
