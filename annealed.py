"""Annealed SYK partition function Z_n(z) = E_g tr exp(z H_SYK(g)) from sampled ED spectra, and its zeros.

Estimator. Write g = r * ghat with r = |g| ~ chi_D independent of the direction ghat. For each sampled g we
diagonalize A = sum_S ghat_S Gamma_S (eigenvalues a_j) and integrate r out exactly with a Gauss rule for chi_D
(Rao-Blackwellization):

    Z(z) = E_ghat (1/N) sum_j E_r cosh(z sigma r a_j),   E_r cosh(t r) = 1F1(D/2; 1/2; t^2/2).

(Z is even since g -> -g flips the spectrum.) The atoms lam = sigma r_q |a_j| with weights w_q / N are stored as
per-bin Taylor moments S[batch, b, p] = sum w (lam - c_b)^p, c_b = (b + 1/2) h, so that

    Z(z) = sum_b sum_p S[b, p] z^p cosh^{(p)}(z c_b) / p!

can be evaluated anywhere afterwards (truncation error ~ (|z| h / 2)^P / P!). Batches give jackknife errors.
"""
import math

import numpy as np


def chi_rule(D, Q):
    """Gauss rule (nodes r_q, weights w_q summing to 1) for r ~ chi_D, exact for polynomials in r^2 of degree 2Q-1."""
    a = D / 2 - 1  # u = r^2/2 ~ Gamma(D/2): generalized Laguerre weight u^a e^{-u}, Golub-Welsch
    i = np.arange(1, Q + 1, dtype=float)
    J = np.diag(2 * i - 1 + a) + np.diag(np.sqrt(i[:-1] * (i[:-1] + a)), 1) + np.diag(np.sqrt(i[:-1] * (i[:-1] + a)), -1)
    u, vec = np.linalg.eigh(J)
    w = vec[0] ** 2
    return np.sqrt(2 * u), w / w.sum()


class BinnedSpectrum:
    """Accumulates S[batch, b, p] for atoms lam >= 0 with weights w; bins grow on demand."""

    def __init__(self, h, P, batches, nbins=64):
        self.h, self.P, self.K = h, P, batches
        self.S = np.zeros((batches, nbins, P))
        self.count = np.zeros(batches, np.int64)

    def add(self, batch, lam, w, nsamples):
        b = np.floor(lam / self.h).astype(np.int64)
        if b.max() >= self.S.shape[1]:
            grow = max(b.max() + 1, 2 * self.S.shape[1]) - self.S.shape[1]
            self.S = np.concatenate([self.S, np.zeros((self.K, grow, self.P))], axis=1)
        d = lam - (b + 0.5) * self.h
        wp = w.copy()
        for p in range(self.P):
            self.S[batch, :, p] += np.bincount(b, weights=wp, minlength=self.S.shape[1])
            wp *= d
        self.count[batch] += nsamples


class AnnealedZ:
    """Evaluate Z(z), Z'(z) from binned moments S (NB, P) normalized by the number of samples."""

    def __init__(self, S_batches, counts, h, meta=None):
        self.Sb = np.asarray(S_batches, float)  # (K, NB, P), unnormalized sums
        self.counts = np.asarray(counts, float)
        self.h = float(h)
        self.meta = dict(meta or {})
        self.K, self.NB, self.P = self.Sb.shape
        self.c = (np.arange(self.NB) + 0.5) * self.h
        self.S = self.Sb.sum(0) / self.counts.sum()
        self.fact = np.array([math.factorial(p) for p in range(self.P + 1)], float)

    @classmethod
    def load(cls, paths, groups=16):
        """Merge sampler outputs (same n, k, h, P); regroup all batches into `groups` jackknife groups."""
        Sbs, counts, meta, h = [], [], None, None
        for path in paths:
            d = np.load(path, allow_pickle=False)
            m = {key: d[key].item() for key in d.files if d[key].ndim == 0}
            if meta is None:
                meta, h = m, float(d["h"])
            assert (m["n"], m["k"], m["P"]) == (meta["n"], meta["k"], meta["P"]) and float(d["h"]) == h, path
            Sbs.append(d["S"])
            counts.append(d["count"])
        NB = max(S.shape[1] for S in Sbs)
        Sb = np.concatenate([np.pad(S, ((0, 0), (0, NB - S.shape[1]), (0, 0))) for S in Sbs])
        counts = np.concatenate(counts)
        keep = counts > 0
        Sb, counts = Sb[keep], counts[keep]
        g = np.arange(len(counts)) % min(groups, len(counts))
        Sg = np.stack([Sb[g == i].sum(0) for i in range(g.max() + 1)])
        cg = np.array([counts[g == i].sum() for i in range(g.max() + 1)])
        meta["samples"] = int(counts.sum())
        meta["files"] = len(paths)
        return cls(Sg, cg, h, meta)

    def jackknife(self, i):
        """Estimator with jackknife group i left out."""
        keep = np.arange(self.K) != i
        return AnnealedZ(self.Sb[keep], self.counts[keep], self.h, self.meta)

    def rebin(self, factor):
        """Exact re-expansion of the moments about the centres of bins `factor` times wider (faster evaluation;
        Taylor error (|z| h / 2)^P / P! grows accordingly)."""
        if factor <= 1:
            return self
        NB2 = -(-self.NB // factor)
        Sb = np.pad(self.Sb, ((0, 0), (0, NB2 * factor - self.NB), (0, 0))).reshape(self.K, NB2, factor, self.P)
        h2 = self.h * factor
        shift = (np.arange(factor) + 0.5) * self.h - h2 / 2  # old centre - new centre
        out = np.zeros((self.K, NB2, self.P))
        for p in range(self.P):
            for q in range(p + 1):
                out[:, :, p] += math.comb(p, q) * np.einsum("kbf,f->kb", Sb[:, :, :, q], shift ** (p - q))
        return AnnealedZ(out, self.counts, h2, self.meta)

    def for_window(self, zmax, tol=0.5):
        """Coarsest rebinning with |z| h / 2 <= tol on |z| <= zmax (tol 0.5: Taylor error < 1e-18 for P = 16)."""
        return self.rebin(max(1, int(2 * tol / (zmax * self.h))))

    def groups(self, z, chunk=4096):
        """Unnormalized values of each jackknife group's sum at z: array (K, len(z)). Z = sum_k / M."""
        z = np.asarray(z, complex).ravel()
        Sk = np.moveaxis(self.Sb, 0, 2)  # (NB, P, K)
        Se, So = Sk[:, 0::2].reshape(self.NB, -1), Sk[:, 1::2].reshape(self.NB, -1)
        pe = np.arange(0, self.P, 2)
        po = np.arange(1, self.P, 2)
        out = np.empty((self.K, len(z)), complex)
        for s in range(0, len(z), chunk):
            zz = z[s:s + chunk]
            E = np.exp(np.outer(zz, self.c))
            Ei = 1 / E
            A = ((E + Ei) / 2 @ Se).reshape(len(zz), len(pe), self.K)
            B = ((E - Ei) / 2 @ So).reshape(len(zz), len(po), self.K)
            out[:, s:s + chunk] = (np.einsum("gpk,gp->kg", A, zz[:, None] ** pe / self.fact[pe])
                                   + np.einsum("gpk,gp->kg", B, zz[:, None] ** po / self.fact[po]))
        return out

    def with_errors(self, z):
        """Z(z) and its jackknife standard error."""
        G = self.groups(z)
        M, Mk = self.counts.sum(), self.counts
        full = G.sum(0) / M
        jk = (G.sum(0)[None] - G) / (M - Mk)[:, None]
        err = np.sqrt((self.K - 1) / self.K * np.sum(np.abs(jk - jk.mean(0)) ** 2, axis=0))
        return full, err

    def __call__(self, z, deriv=False, chunk=8192):
        z = np.asarray(z, complex)
        shape = z.shape
        z = z.ravel()
        out = np.empty(len(z), complex)
        dout = np.empty(len(z), complex) if deriv else None
        Se, So = self.S[:, 0::2], self.S[:, 1::2]  # even / odd p
        if deriv:
            cSe, cSo = self.c[:, None] * Se, self.c[:, None] * So
        pe = np.arange(0, self.P, 2)
        po = np.arange(1, self.P, 2)
        for s in range(0, len(z), chunk):
            zz = z[s:s + chunk]
            E = np.exp(np.outer(zz, self.c))
            Ei = 1 / E
            ch, sh = (E + Ei) / 2, (E - Ei) / 2
            A, B = ch @ Se, sh @ So  # sum_b S[b, p] cosh^{(p)}(z c_b)
            out[s:s + chunk] = (A * zz[:, None] ** pe / self.fact[pe]).sum(1) + (B * zz[:, None] ** po / self.fact[po]).sum(1)
            if deriv:
                # d/dz [z^p f_p(zc)/p!] = z^{p-1} f_p(zc)/(p-1)! + c z^p f_{p+1}(zc)/p!
                d1 = (B * zz[:, None] ** (po - 1) / self.fact[po - 1]).sum(1)
                if len(pe) > 1:
                    d1 += (A[:, 1:] * zz[:, None] ** (pe[1:] - 1) / self.fact[pe[1:] - 1]).sum(1)
                d2 = ((sh @ cSe) * zz[:, None] ** pe / self.fact[pe]).sum(1) + ((ch @ cSo) * zz[:, None] ** po / self.fact[po]).sum(1)
                dout[s:s + chunk] = d1 + d2
        if deriv:
            return out.reshape(shape), dout.reshape(shape)
        return out.reshape(shape)

    def moment(self, two_l):
        """E tr H^{2l} from the binned atoms (exact sum of (c + d)^{2l})."""
        tot = 0.0
        for p in range(min(two_l, self.P - 1) + 1):
            tot += math.comb(two_l, p) * np.sum(self.c ** (two_l - p) * self.S[:, p])
        return tot


# ---------------------------------------------------------------------------------------------- zero finding

def _contour(rect, n_x, n_y):
    x0, x1, y0, y1 = rect
    tx = np.linspace(0, 1, n_x, endpoint=False)
    ty = np.linspace(0, 1, n_y, endpoint=False)
    return np.concatenate([x0 + (x1 - x0) * tx + 1j * y0, x1 + 1j * (y0 + (y1 - y0) * ty),
                           x1 - (x1 - x0) * tx + 1j * y1, x0 + 1j * (y1 - (y1 - y0) * ty)])


def winding(f, rect, n_edge=48, max_rounds=30, dphi=0.5):
    """Number of zeros of f inside rect = (x0, x1, y0, y1), by the argument principle with adaptive refinement of
    the boundary until consecutive phase steps are below dphi. For an AnnealedZ the initial sampling resolves the
    fastest phase rotation away from zeros, |d arg Z / dz| <= max eigenvalue atom (else steps of 2 pi alias)."""
    rate = (f.c[-1] + f.h) if hasattr(f, "c") else 0.0
    zs = _contour(rect, max(n_edge, int(4 * rate * (rect[1] - rect[0]))), max(n_edge, int(4 * rate * (rect[3] - rect[2]))))
    fs = f(zs)
    for _ in range(max_rounds):
        steps = np.angle(np.roll(fs, -1) / fs)
        bad = np.nonzero(np.abs(steps) > dphi)[0]
        if len(bad) == 0:
            break
        mid = (zs[bad] + np.roll(zs, -1)[bad]) / 2
        fm = f(mid)
        zs = np.insert(zs, bad + 1, mid)
        fs = np.insert(fs, bad + 1, fm)
    total = np.sum(np.angle(np.roll(fs, -1) / fs)) / (2 * np.pi)
    return int(round(total)), total, np.min(np.abs(fs))


def newton(F, z0, atol, maxit=60):
    """Newton iteration; converged once a step is below atol (choose atol above the rounding floor of F, which at
    large |Re z| is ~1e-9 |z| because sum_b S_b cosh(z c_b) cancels)."""
    z = complex(z0)
    for _ in range(maxit):
        v, dv = F(np.array([z]), deriv=True)
        step = v[0] / dv[0]
        z -= step
        if abs(step) < atol:
            return z, True
    return z, False


def find_zeros(F, rect, min_size, split=0.4813, _depth=0):
    """All zeros of F in rect by recursive argument-principle subdivision + Newton (to min_size / 10).
    Returns list of (z, multiplicity); boxes below min_size with winding > 1 are reported as multiple zeros."""
    x0, x1, y0, y1 = rect
    nz, _, _ = winding(F, rect)
    if nz <= 0:
        return []
    if nz == 1:
        z, ok = newton(F, complex((x0 + x1) / 2, (y0 + y1) / 2), atol=min_size / 10)
        sx, sy = 1e-9 * (x1 - x0), 1e-9 * (y1 - y0)
        if ok and x0 - sx <= z.real <= x1 + sx and y0 - sy <= z.imag <= y1 + sy:
            return [(z, 1)]
    if max(x1 - x0, y1 - y0) < min_size:
        return [(complex((x0 + x1) / 2, (y0 + y1) / 2), nz)]
    xm, ym = x0 + split * (x1 - x0), y0 + split * (y1 - y0)
    if (x1 - x0) > 2 * (y1 - y0):
        subs = [(x0, xm, y0, y1), (xm, x1, y0, y1)]
    elif (y1 - y0) > 2 * (x1 - x0):
        subs = [(x0, x1, y0, ym), (x0, x1, ym, y1)]
    else:
        subs = [(x0, xm, y0, ym), (xm, x1, y0, ym), (x0, xm, ym, y1), (xm, x1, ym, y1)]
    out = []
    for r in subs:
        out += find_zeros(F, r, min_size, split, _depth + 1)
    return out


def axis_zeros(F, Y, n_grid=20000):
    """Zeros of the real function y -> Z(iy) on (0, Y] (sign changes + bisection)."""
    from scipy.optimize import brentq
    ys = np.linspace(Y / n_grid, Y, n_grid)
    v = F(1j * ys).real
    roots = []
    for i in np.nonzero(np.sign(v[:-1]) * np.sign(v[1:]) < 0)[0]:
        roots.append(brentq(lambda y: F(np.array([1j * y]))[0].real, ys[i], ys[i + 1], xtol=1e-15, rtol=1e-15))
    return np.array(roots)


def cell_winding(Zg):
    """Winding number of the phase of Zg (rows = increasing Im, cols = increasing Re) around each grid cell."""
    ph = np.angle(Zg)
    d = lambda a, b: np.angle(np.exp(1j * (b - a)))
    w = (d(ph[:-1, :-1], ph[:-1, 1:]) + d(ph[:-1, 1:], ph[1:, 1:]) + d(ph[1:, 1:], ph[1:, :-1])
         + d(ph[1:, :-1], ph[:-1, :-1])) / (2 * np.pi)
    return np.rint(w).astype(int)


def quadrant_zeros(F, X, Y, dx, dy, Zg=None):
    """Zeros with Re >= 0 and 0 < Im <= Y, |Re| <= X. The grid x_j = (j - 1/2) dx, j = 0.., straddles the imaginary
    axis so that axis zeros are cell interiors; cells with winding 1 are polished by Newton, others re-searched by
    adaptive argument-principle subdivision. Mirror images (Z is even and real on R) are folded away.
    Returns (zeros, on_axis flags, xs, ys, Zg)."""
    xs = (np.arange(int(np.ceil(X / dx)) + 2) - 0.5) * dx
    ys = np.linspace(1e-3 * dy, Y, int(np.ceil(Y / dy)) + 1)
    if Zg is None:
        Zg = F(xs[None, :] + 1j * ys[:, None])
    W = cell_winding(Zg)
    found = []
    for i, j in zip(*np.nonzero(W)):
        rect = (xs[j], xs[j + 1], ys[i], ys[i + 1])
        z = None
        if W[i, j] == 1:
            z, ok = newton(F, complex((xs[j] + xs[j + 1]) / 2, (ys[i] + ys[i + 1]) / 2), atol=1e-7 * min(dx, dy))
            if not (ok and rect[0] - 1e-9 * dx <= z.real <= rect[1] + 1e-9 * dx and rect[2] - 1e-9 * dy <= z.imag <= rect[3] + 1e-9 * dy):
                z = None
        if z is not None:
            found.append(z)
        else:
            found += [r for r, mult in find_zeros(F, rect, min_size=1e-6 * min(dx, dy)) for _ in range(mult)]
    # verify the count in horizontal strips by the adaptive argument principle; re-search mismatched strips
    strip = 24
    for i0 in range(0, len(ys) - 1, strip):
        i1 = min(i0 + strip, len(ys) - 1)
        rect = (xs[0], xs[-1], ys[i0], ys[i1])
        inside = [z for z in found if rect[0] <= z.real <= rect[1] and rect[2] <= z.imag < rect[3]]
        nw, _, _ = winding(F, rect)
        if nw != len(inside):
            found = [z for z in found if z not in inside]
            found += [r for r, mult in find_zeros(F, rect, min_size=1e-6 * min(dx, dy)) for _ in range(mult)]
    scale = max(X, Y)
    zs, axis = [], []
    for z in sorted(found, key=lambda z: (z.imag, z.real)):
        if abs(z.real) < 1e-5 * dx:
            z = complex(0.0, z.imag)
        elif z.real < 0:
            continue
        if any(abs(z - u) < 1e-5 * dx for u in zs):
            continue
        zs.append(z)
        axis.append(z.real == 0.0)
    return np.array(zs, complex), np.array(axis, bool), xs, ys, Zg


def jackknife_zeros(Z, zeros):
    """Jackknife standard errors (Re, Im) of each zero: Newton on every leave-one-group-out estimator."""
    jk = np.array([[newton(Z.jackknife(i), z, atol=1e-8 * abs(z))[0] for z in zeros] for i in range(Z.K)]).reshape(Z.K, len(zeros))
    f = (Z.K - 1) / Z.K
    jk_re = np.abs(jk.real)  # a mirror image is the same zero
    return np.sqrt(f * np.sum((jk_re - jk_re.mean(0)) ** 2, 0)), np.sqrt(f * np.sum((jk.imag - jk.imag.mean(0)) ** 2, 0))
