"""Checks of the Majorana construction, the fast Hamiltonian builder and the annealed estimator.

    python -m pytest tests -q      (or: python tests/test_syk.py)
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from annealed import AnnealedZ, BinnedSpectrum, chi_rule  # noqa: E402
from sykmaj import SYK, majorana_paulis, pauli_matrix  # noqa: E402


def test_majoranas_anticommute():
    for n in (4, 6, 8):
        N = 1 << (n // 2)
        G = [pauli_matrix(x, z, p, N) for x, z, p in zip(*majorana_paulis(n))]
        for i in range(n):
            assert np.allclose(G[i], G[i].conj().T)
            for j in range(n):
                assert np.allclose(G[i] @ G[j] + G[j] @ G[i], 2 * (i == j) * np.eye(N))


def test_gamma_is_product_of_majoranas():
    n, k = 8, 4
    syk = SYK(n, k)
    N = syk.N
    G = [pauli_matrix(x, z, p, N) for x, z, p in zip(*majorana_paulis(n))]
    for d in [0, 7, 33, syk.D - 1]:
        prod = np.eye(N)
        for i in syk.S[d]:
            prod = prod @ G[i]
        assert np.allclose(1j ** (k * (k - 1) // 2) * prod, pauli_matrix(syk.x[d], syk.z[d], syk.p[d], N))


def test_blocks_match_dense():
    rng = np.random.default_rng(0)
    for n, k in [(6, 2), (8, 4), (10, 4), (10, 6)]:
        syk = SYK(n, k)
        g = rng.standard_normal(syk.D)
        Hd = syk.dense(g)
        assert np.allclose(Hd, Hd.conj().T)
        ev = np.sort(syk.eigvals(g[None])[0])
        assert np.allclose(ev, np.linalg.eigvalsh(Hd))
        assert np.isclose(np.mean(ev ** 2), np.sum(g ** 2))  # tr Gamma_S Gamma_T = delta_ST


def test_chi_rule():
    from scipy.special import hyp1f1
    for D in (28, 70, 495):
        r, w = chi_rule(D, 16)
        for t in (0.5, 2.0, 0.3 + 0.2j):
            exact = hyp1f1(D / 2, 0.5, (t * 2 / math.sqrt(D)) ** 2 / 2) if np.isreal(t) else None
            got = np.sum(w * np.cosh(t * 2 / math.sqrt(D) * r))
            if exact is not None:
                assert np.isclose(got, exact, rtol=1e-12)


def test_estimator_moments_and_eval():
    rng = np.random.default_rng(1)
    syk = SYK(8, 4)
    r_q, w_q = chi_rule(syk.D, 16)
    acc = BinnedSpectrum(0.01 * math.sqrt(syk.V), 14, 4)
    lam_all = []
    for it in range(8):
        g = rng.standard_normal((500, syk.D))
        a = syk.eigvals(g) / np.linalg.norm(g, axis=1)[:, None]
        lam = (syk.sigma * np.abs(a).ravel()[:, None] * r_q[None, :]).ravel()
        w = np.broadcast_to(w_q[None, :] / syk.N, (a.size, 16)).ravel()
        acc.add(it % 4, lam, w, 500)
        lam_all.append((lam, w))
    Z = AnnealedZ(acc.S, acc.count, acc.h)
    assert np.isclose(Z.moment(0), 1.0)
    assert np.isclose(Z.moment(2), syk.V, rtol=1e-12)  # second moment is exact under the radial rule
    lam = np.concatenate([l for l, _ in lam_all])
    w = np.concatenate([w for _, w in lam_all]) / 4000
    zs = np.array([0.3 + 2j, 1.5 + 7j, 4j, 2.0 - 1j]) / math.sqrt(syk.V)
    direct = np.array([np.sum(w * np.cosh(z * lam)) for z in zs])
    val, dval = Z(zs, deriv=True)
    assert np.allclose(val, direct, rtol=1e-11, atol=1e-13)
    ddirect = np.array([np.sum(w * lam * np.sinh(z * lam)) for z in zs])
    assert np.allclose(dval, ddirect, rtol=1e-10, atol=1e-12)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
