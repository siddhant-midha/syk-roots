"""Exact annealed partition function of SYK at n = 8, k = 4 via the traceless GOE(8).

For n = 8 the Majoranas have a real representation, so every Gamma_S (|S| = 4) is real symmetric, and in each
fermion-parity sector Gamma_{S^c} = +-Gamma_S. The 35 operators {Gamma_S : 1 in S} restricted to a sector are
orthonormal for the normalized trace and traceless, and 35 = dim{traceless real symmetric 8x8}. Hence H restricted to
each sector is sum_S (g_S +- g_{S^c}) Gamma_S, an orthogonally invariant Gaussian on traceless symmetric matrices, i.e.
a traceless GOE(8), and the two sectors are independent. So

    Z_{8,4}(z) = Z0(c z),   Z0(w) = E (1/8) tr exp(w H0),   H0 = H - (tr H / 8) I,   H ~ GOE(8) ~ exp(-tr H^2 / 2),

with c^2 = V_8 / E[(1/8) tr H0^2] = V_8 / ((N + 1)/2 - 1/N). For the GOE, H0 and tr H are independent, so
Z0(w) = Z_GOE(w) exp(-w^2 / (2N)). The GOE one-level density is known exactly (Mehta, N even),

    rho_N(x) = sum_{k<N} phi_k(x)^2 + sqrt(N/2) phi_{N-1}(x) int eps(x - t) phi_N(t) dt,   eps = sgn / 2,

(phi_k: Hermite functions), i.e. rho_N = exp(-x^2) A(x) + exp(-x^2/2) erf(x/sqrt 2) B(x) with polynomials A, B, whose
Laplace transforms follow from int x^m e^{-x^2 + w x} = d^m/dw^m sqrt(pi) e^{w^2/4} and
int x^m e^{-x^2/2 + w x} erf(x/sqrt 2) = d^m/dw^m sqrt(2 pi) e^{w^2/2} erf(w/2). For N = 8 (c^2 = 3/128):

    Z0(w) = e^{-w^2/16} / 92160 * [ e^{w^2/4} (w^10 + 90 w^8 + 6540 w^6 + 62520 w^4 + 159120 w^2 + 92160)
                                   + sqrt(pi) e^{w^2/2} erf(w/2) (1920 w^7 + 20160 w^5 + 50400 w^3 + 25200 w) ].

(A direct quadrature of rho_N fails at large Im w: Z0(i y) ~ e^{-3 y^2/16} poly(y) is below its rounding noise.)
"""
import math

import mpmath as mp
import numpy as np
import sympy as sp

N = 8


def density_parts(n=N):
    """Polynomials A, B with rho_n(x) = exp(-x^2) A(x) + exp(-x^2/2) erf(x/sqrt 2) B(x)."""
    x = sp.symbols("x", real=True)
    g = sp.exp(-x ** 2 / 2)
    phi = [g * sp.hermite(k, x) / sp.sqrt(2 ** k * sp.factorial(k) * sp.sqrt(sp.pi)) for k in range(n + 1)]
    # I_k = int_{-inf}^x phi_k = alpha_k * I_0 - (sum of phi's) for even k (recurrence from phi_k' relation)
    alpha = {0: sp.Integer(1)}
    rest = {0: sp.Integer(0)}
    for k in range(1, n, 2):
        alpha[k + 1] = sp.sqrt(sp.Rational(k, k + 1)) * alpha[k - 1]
        rest[k + 1] = sp.sqrt(sp.Rational(k, k + 1)) * rest[k - 1] + sp.sqrt(sp.Rational(2, k + 1)) * phi[k]
    # J = I_n - I_n(inf)/2 = alpha_n (I_0 - I_0(inf)/2) - rest_n,  I_0 - I_0(inf)/2 = pi^{-1/4} sqrt(pi/2) erf(x/sqrt2)
    pref = sp.sqrt(sp.Rational(n, 2)) * phi[n - 1]
    A = sp.expand(sp.simplify((sum(p ** 2 for p in phi[:n]) - pref * rest[n]) * sp.exp(x ** 2)))
    B = sp.expand(sp.simplify(pref * alpha[n] * sp.pi ** sp.Rational(-1, 4) * sp.sqrt(sp.pi / 2) * sp.exp(x ** 2 / 2)))
    return x, sp.Poly(A, x), sp.Poly(B, x)


def z0_expression(n=N):
    """Sympy expression of Z0(w) = E (1/n) tr exp(w H0) for the traceless GOE(n), and the symbol w."""
    x, A, B = density_parts(n)
    w = sp.symbols("w")
    G1 = sp.sqrt(sp.pi) * sp.exp(w ** 2 / 4)
    G2 = sp.sqrt(2 * sp.pi) * sp.exp(w ** 2 / 2) * sp.erf(w / 2)
    zg = sum(c * sp.diff(G1, w, m) for (m,), c in A.terms()) + sum(c * sp.diff(G2, w, m) for (m,), c in B.terms())
    return sp.simplify(zg * sp.exp(-w ** 2 / (2 * n)) / n), w


class ExactZ84:
    """Z_{8,4}(z) = E tr exp(z H_SYK) exactly (normalized trace, paper's normalization), for complex arrays z.
    Evaluated with numpy/scipy (Faddeeva-based complex erf); precise=True uses mpmath at `dps` digits instead."""

    def __init__(self, V, precise=False, dps=40):
        expr, w = z0_expression(N)
        mods = "mpmath" if precise else ["scipy", "numpy"]
        self.f = sp.lambdify(w, expr, mods)
        self.df = sp.lambdify(w, sp.diff(expr, w), mods)
        self.scale = math.sqrt(V / ((N + 1) / 2 - 1 / N))  # H_SYK ~ scale * H0 on each parity sector
        self.phase_rate = self.scale * 12.0  # |d arg Z / dz| bound for contour sampling (spectrum of H0 within ~ +-12)
        self.precise, self.dps = precise, dps

    def _eval(self, fun, z):
        z = np.asarray(z, complex)
        if not self.precise:
            return np.asarray(fun(z * self.scale), complex)
        with mp.workdps(self.dps):
            out = [complex(fun(mp.mpc(v.real, v.imag) * self.scale)) for v in z.ravel()]
        return np.array(out, complex).reshape(z.shape)

    def __call__(self, z, deriv=False):
        val = self._eval(self.f, z)
        if deriv:
            return val, self.scale * self._eval(self.df, z)
        return val
