"""Zeros of the annealed SYK partition function for one (n, k) from the sampler outputs.

    python analyze.py out/n8_k4/*.npz --out results/n8_k4

Works in the scaled variable t = z sqrt(V_n) (unit-variance Hamiltonian; the q-oscillator comparator Q_n has its zeros
at +-i t_j(q_n)). Steps:
  1. window: Im t <= Y = max(36, min(2.5 sqrt(2D), 120)), 0 <= Re t <= max(18, Y/2) (D = C(n,k); the radial chi_D
     smearing of the coupling scale damps the oscillations of Z(it) beyond |t| ~ sqrt(2D), where zeros leave the axis);
  2. jackknife noise map rel(z) = std_jk Z(z) / rms|Z| (rms over a small neighbourhood, so zeros do not blow it up);
  3. all zeros in the first quadrant of the window (grid winding + Newton, argument-principle checked);
  4. per zero: jackknife errors of Re and Im, local noise (std_jk Z / rms|Z| on a ring of radius dmin/2, dmin =
     distance to the nearest other zero, mirror images included), the zero recomputed from two disjoint halves of the
     samples, and a 'resolved' flag: noise < 0.4, jackknife error < dmin/5, halves differ by < dmin/5 and, off the
     axis, Re t > 5 err(Re t).
Writes <out>.json (zeros, flags, q-oscillator zeros, meta) and <out>.npz (log|Z| and arg Z on the plotting grid,
noise map).
"""
import argparse
import json
import math
import time

import numpy as np
from scipy.ndimage import uniform_filter

from annealed import AnnealedZ, jackknife_zeros, newton, quadrant_zeros, winding
from qosc import qn_paper, t_zeros

REL_MAX = 0.4
COMMON = (18.0, 36.0)  # common comparison window |Re t| <= 18, |Im t| <= 36 (always inside the analysis window)


def local_noise(Z, z, radius):
    ring = z + radius * np.exp(2j * np.pi * np.arange(8) / 8)
    val, err = Z.with_errors(ring)
    return float(np.mean(err) / np.sqrt(np.mean(np.abs(val) ** 2)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ymax", type=float, default=0, help="window Im t (default max(36, min(2.5 sqrt(2D), 120)))")
    ap.add_argument("--xmax", type=float, default=0, help="window Re t (default max(18, ymax/2))")
    ap.add_argument("--dt", type=float, default=0, help="grid step in t (default max(0.06, ymax/1500))")
    ap.add_argument("--plot-px", type=int, default=600, help="plotting grid rows per half-plane")
    args = ap.parse_args()

    t0 = time.time()
    Z0 = AnnealedZ.load(args.files, groups=16)
    m = Z0.meta
    n, k, D, V = m["n"], m["k"], m["D"], m["V"]
    s = math.sqrt(V)
    Y = args.ymax or max(COMMON[1], min(2.5 * math.sqrt(2 * D), 120.0))
    X = args.xmax or max(COMMON[0], Y / 2)
    dt = args.dt or max(0.06, Y / 1500)
    Z = Z0.for_window(math.hypot(X, Y) / s)
    print(f"n={n} k={k} D={D} V={V:.6g} samples={m['samples']} files={m['files']}  window Re t<={X:.1f} Im t<={Y:.1f} "
          f"dt={dt:.3f}  bins {Z0.NB}->{Z.NB}", flush=True)

    # noise map on a coarse grid (t units)
    nx, ny = 61, 121
    xs_n, ys_n = np.linspace(0, X, nx), np.linspace(0, Y, ny)
    val, err = Z.with_errors((xs_n[None, :] + 1j * ys_n[:, None]).ravel() / s)
    val, err = val.reshape(ny, nx), err.reshape(ny, nx)
    rel = err / np.sqrt(uniform_filter(np.abs(val) ** 2, 5, mode="nearest"))
    print(f"  noise map done {time.time() - t0:.0f}s", flush=True)

    # zeros (in z units internally)
    zs, axis, xs, ys, Zg = quadrant_zeros(Z, X / s, Y / s, dt / s, dt / s)
    nwind, _, _ = winding(Z, (xs[0], xs[-1], ys[0], ys[-1]))
    print(f"  {len(zs)} zeros ({axis.sum()} on axis) in quadrant; argument principle on [-dx/2, X] x [0, Y]: {nwind} "
          f"(mirror images with Re in (-dx/2, 0) included there)  {time.time() - t0:.0f}s", flush=True)
    ex, ey = jackknife_zeros(Z, zs)
    halves = [AnnealedZ(Z.Sb[sl], Z.counts[sl], Z.h) for sl in (slice(0, Z.K // 2), slice(Z.K // 2, Z.K))]
    zh = np.array([[newton(H, z, atol=1e-8 * abs(z))[0] for z in zs] for H in halves]).reshape(2, len(zs))
    zh = np.abs(zh.real) + 1j * zh.imag
    allz = np.concatenate([zs, -zs.conj()])
    rows = []
    for i, (z, on) in enumerate(zip(zs, axis)):
        others = np.abs(allz - z)
        others = others[others > 1e-12 * abs(z)]
        dmin = float(others.min()) if len(others) else float("inf")
        noise = local_noise(Z, z, 0.5 * dmin)  # ~ (zero displacement) / (dmin / 2)
        e = math.hypot(ex[i], ey[i])
        dh = abs(zh[0, i] - zh[1, i])
        resolved = bool(noise < REL_MAX and e < 0.2 * dmin and dh < 0.2 * dmin and (on or z.real > 5 * ex[i]))
        rows.append(dict(t_re=z.real * s, t_im=z.imag * s, z_re=z.real, z_im=z.imag, on_axis=bool(on),
                         err_re_t=ex[i] * s, err_im_t=ey[i] * s, noise=noise, nearest_t=dmin * s,
                         halves_t=[[w.real * s, w.imag * s] for w in zh[:, i]], resolved=resolved))
    off = [r for r in rows if not r["on_axis"] and r["resolved"]]
    onset = min((math.hypot(r["t_re"], r["t_im"]) for r in off), default=None)
    n_axis_before = sum(1 for r in rows if r["on_axis"] and (onset is None or r["t_im"] < onset))
    print(f"  resolved: {sum(r['resolved'] for r in rows)} / {len(rows)}; off-axis resolved {len(off)}; "
          f"first resolved off-axis |t| = {onset}  {time.time() - t0:.0f}s", flush=True)

    q_n = qn_paper(n, k)
    tq = t_zeros(q_n, Y)
    tq_exact = t_zeros(m["q_exact"], Y)

    # plotting grids (mirror-symmetric, z units): the full analysis window and the common comparison window
    grids = {}
    for tag, (gx, gy), py in [("", (X, Y), args.plot_px), ("c_", COMMON, args.plot_px // 2)]:
        px = max(2, int(round(py * gx / gy)))
        xp = np.linspace(-gx, gx, 2 * px) / s
        yp = np.linspace(-gy, gy, 2 * py) / s
        Zp = Z(xp[None, :] + 1j * yp[:, None])
        grids.update({tag + "x": xp, tag + "y": yp, tag + "logabs": np.log(np.abs(Zp) + 1e-300).astype(np.float32),
                      tag + "phase": np.angle(Zp).astype(np.float32)})
    np.savez_compressed(args.out + ".npz", rel=rel, rel_x=xs_n / s, rel_y=ys_n / s, **grids)
    res = dict(n=n, k=k, D=D, V=V, sqrt_V=s, sigma=m["sigma"], q_n=q_n, q_exact=m["q_exact"], alpha=k * k / n,
               samples=m["samples"], files=m["files"], window_t=[X, Y], dt=dt, sqrt_2D=math.sqrt(2 * D),
               argument_principle=nwind, onset_t=onset, n_axis_before_onset=n_axis_before,
               qosc_t=tq.tolist(), qosc_exact_t=tq_exact.tolist(), zeros=rows, wall=time.time() - t0)
    with open(args.out + ".json", "w") as f:
        json.dump(res, f, indent=1)
    print(f"  wrote {args.out}.json/.npz  {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
