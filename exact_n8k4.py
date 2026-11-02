"""Exact zeros of the n = 8, k = 4 annealed partition function (traceless GOE(8), see goe.py) and their comparison
with the ED/Monte Carlo zeros of analyze.py.

    python exact_n8k4.py results/n8_k4.json --out results/exact_n8_k4

On the imaginary axis Z0(i y) = e^{-3y^2/16} G(y) / 92160 with G(y) = P(-y^2) - 2 D(y/2) R(y) real (D = Dawson
function), P(-y^2) = -y^10 + 90 y^8 - 6540 y^6 + 62520 y^4 - 159120 y^2 + 92160 and
R(y) = -1920 y^7 + 20160 y^5 - 50400 y^3 + 25200 y. Since |D| < 0.55, -y^10 dominates for large y, so the sign
changes of G on (0, Y_MAX] are all the zeros on the positive imaginary axis.
"""
import argparse
import json
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import dawsn

from annealed import quadrant_zeros
from goe import N, ExactZ84

Y_MAX = 400.0


def axis_zeros_exact():
    P = lambda y: -y ** 10 + 90 * y ** 8 - 6540 * y ** 6 + 62520 * y ** 4 - 159120 * y ** 2 + 92160
    R = lambda y: -1920 * y ** 7 + 20160 * y ** 5 - 50400 * y ** 3 + 25200 * y
    G = lambda y: P(y) - 2 * dawsn(y / 2) * R(y)
    y = np.linspace(1e-6, Y_MAX, 4_000_001)
    g = G(y)
    idx = np.nonzero(np.sign(g[:-1]) * np.sign(g[1:]) < 0)[0]
    return np.array([brentq(G, y[i], y[i + 1], xtol=1e-15) for i in idx])  # in units of w = c z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ed_json")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    ed = json.load(open(args.ed_json))
    V = ed["V"]
    s = math.sqrt(V)
    ex = ExactZ84(V)
    t_per_w = math.sqrt((N + 1) / 2 - 1 / N)  # t = sqrt(V) z = w sqrt(E tr H0^2 / N)
    axis_t = axis_zeros_exact() * t_per_w
    X, Y = ed["window_t"]
    zs, on, _, _, _ = quadrant_zeros(ex, X / s, Y / s, 0.04 / s, 0.04 / s)
    zeros = [dict(t_re=z.real * s, t_im=z.imag * s, z_re=z.real, z_im=z.imag, on_axis=bool(a)) for z, a in zip(zs, on)]
    comp = []
    for e in ed["zeros"]:
        t = complex(e["t_re"], e["t_im"])
        d = np.abs(zs * s - t)
        i = int(d.argmin())
        err = math.hypot(e["err_re_t"], e["err_im_t"])
        comp.append(dict(ed_t=[t.real, t.imag], exact_t=[zs[i].real * s, zs[i].imag * s], dist_t=float(d[i]),
                         sigma=float(d[i] / max(err, 1e-300)), ed_resolved=e["resolved"]))
    res = dict(n=8, k=4, V=V, sqrt_V=s, c_squared=V / ((N + 1) / 2 - 1 / N), window_t=[X, Y],
               axis_zeros_t_all=axis_t.tolist(), zeros=zeros, ed_comparison=comp)
    with open(args.out + ".json", "w") as f:
        json.dump(res, f, indent=1)
    lines = ["### n = 8, k = 4: exact zeros (traceless GOE(8)) against ED", "",
             f"All zeros of the exact Z on the positive imaginary axis (t = sqrt(V) z): {len(axis_t)}, at "
             + ", ".join(f"{t:.4f}" for t in axis_t) + ". All other zeros are off the axis.", "",
             f"Zeros in the window |Re t| <= {X:g}, 0 < Im t <= {Y:g} (first quadrant):", "",
             "| exact t | exact z | nearest ED zero (t) | distance / jackknife error | ED resolved |",
             "|---|---|---|---|---|"]
    for z in zeros:
        t = complex(z["t_re"], z["t_im"])
        best = min(comp, key=lambda c: abs(complex(*c["ed_t"]) - t))
        near = abs(complex(*best["ed_t"]) - t) < 0.5
        lines.append(f"| {t.real:.4f} + {t.imag:.4f}i | {z['z_re']:.3f} + {z['z_im']:.3f}i | "
                     + (f"{best['ed_t'][0]:.4f} + {best['ed_t'][1]:.4f}i | {best['sigma']:.1f} | "
                        f"{'yes' if best['ed_resolved'] else 'no'} |" if near else "none | | |"))
    spurious = [c for c in comp if c["dist_t"] > 0.5]
    lines += ["", f"Every resolved ED zero lies within {max(c['sigma'] for c in comp if c['ed_resolved']):.1f} "
              f"jackknife errors of an exact zero. {len(spurious)} unresolved ED zeros have no exact counterpart "
              "(Monte Carlo noise; all were flagged unresolved)."
              if all(not c["ed_resolved"] for c in spurious) else "WARNING: a resolved ED zero has no exact counterpart."]
    with open(args.out + ".md", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
