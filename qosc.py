"""Vacuum partition function of the q-deformed oscillator, L_q(z) = <0| exp(z (a_q + a_q^*)) |0> (paper 1, Thm 1.1).

T_q = a_q + a_q^* is the Jacobi matrix with off-diagonal sqrt([m]_q); truncating it at size K gives the Gauss rule
(theta_i, w_i) of the q-Gaussian measure, exact for the first 2K-1 moments, so L_q(z) = sum_i w_i exp(z theta_i).
Its zeros are +-i t_j(q); the paper's finite-n comparator is Q_n(z) = L_{q_n}(sqrt(V_n) z), q_n = exp(-2 k^2/n),
with zeros +-i t_j(q_n) / sqrt(V_n).
"""
import math

import numpy as np
from scipy.optimize import brentq

S_STAR = 1.19967864025773  # s tanh s = 1
ALPHA_STAR = 2 * S_STAR / math.cosh(S_STAR)  # ~ 1.3255 (Theorem 1.3)


def gauss_rule(q, K=600):
    m = np.arange(1, K)
    qint = m.astype(float) if q == 1 else (1 - q ** m) / (1 - q)
    off = np.sqrt(qint)
    theta, vec = np.linalg.eigh(np.diag(off, 1) + np.diag(off, -1))
    return theta, vec[0] ** 2


def L_q(z, q, K=600):
    theta, w = gauss_rule(q, K)
    z = np.asarray(z, complex)
    return (w * np.cosh(np.multiply.outer(z, theta))).sum(-1)


def t_zeros(q, tmax, K=600, n_grid=40000, floor=1e-11):
    """Zeros 0 < t_1 < t_2 < ... <= tmax of phi_q(t) = L_q(i t) = sum_i w_i cos(t theta_i). Sign changes where
    |phi| is below `floor` are double-precision noise (phi decays fast for q near 1) and are skipped."""
    theta, w = gauss_rule(q, K)
    phi = lambda t: np.cos(np.atleast_1d(t)[:, None] * theta) @ w
    ts = np.linspace(tmax / n_grid, tmax, n_grid)
    v = phi(ts)
    out = []
    for i in np.nonzero(v[:-1] * v[1:] < 0)[0]:
        if max(abs(v[i]), abs(v[i + 1])) < floor:
            break
        out.append(brentq(lambda t: phi(t)[0], ts[i], ts[i + 1], xtol=1e-14))
    return np.array(out)


def qn_paper(n, k):
    return math.exp(-2 * k * k / n)
