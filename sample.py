"""Sample SYK disorder, diagonalize, and accumulate the annealed partition function (see annealed.py).

    python sample.py --n 8 --k 4 --samples 1000000 --out data/n8_k4_s1.npz --seed 1

Writes (atomically, also as periodic checkpoints) an npz with the binned moments S[batch, bin, p], per-batch
sample counts, and per-batch sums of the angular moments mu_l = (1/N) sum_j a_j^{2l}, l = 1..4, used to check
E tr H^4 against its exact value V^2 (2 + q_exact).
"""
import argparse
import os
import time

import numpy as np

from annealed import BinnedSpectrum, chi_rule
from sykmaj import SYK


def save(path, syk, acc, mom, args, wall):
    tmp = path + ".tmp.npz"
    np.savez(tmp, S=acc.S, count=acc.count, mom=mom, h=acc.h, P=acc.P, Q=args.Q, n=syk.n, k=syk.k, D=syk.D,
             N=syk.N, sigma=syk.sigma, V=syk.V, q_exact=syk.q_exact(), seed=args.seed, wall=wall)
    os.replace(tmp, path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--samples", type=int, required=True, help="number of disorder samples (upper bound if --time-limit)")
    ap.add_argument("--chunk", type=int, default=0, help="samples per ED batch (default: by size)")
    ap.add_argument("--batches", type=int, default=16, help="jackknife batches")
    ap.add_argument("--Q", type=int, default=16, help="radial chi_D Gauss nodes")
    ap.add_argument("--P", type=int, default=16, help="Taylor moments per bin")
    ap.add_argument("--h", type=float, default=0, help="bin width in units of sqrt(V_n); default min(0.01, 0.15/sqrt(2D)), "
                    "which keeps the Taylor error below 1e-20 for |z| sqrt(V_n) <= 4 sqrt(2D)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--time-limit", type=float, default=0, help="stop after this many seconds (0 = none)")
    ap.add_argument("--checkpoint", type=float, default=600, help="seconds between checkpoints")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    t0 = time.time()
    syk = SYK(args.n, args.k)
    r_q, w_q = chi_rule(syk.D, args.Q)
    chunk = args.chunk or max(1, min(20000, 2 ** 22 // syk.N ** 2))
    rng = np.random.default_rng([args.seed, args.n, args.k])
    h_rel = args.h or min(0.01, 0.15 / np.sqrt(2 * syk.D))
    acc = BinnedSpectrum(h_rel * np.sqrt(syk.V), args.P, args.batches)
    mom = np.zeros((args.batches, 4))
    print(f"n={syk.n} k={syk.k} D={syk.D} N={syk.N} sigma={syk.sigma:.6g} V={syk.V:.6g} q_exact={syk.q_exact():.6g} "
          f"chunk={chunk} setup={time.time() - t0:.1f}s", flush=True)
    done, it, last = 0, 0, time.time()
    while done < args.samples:
        m = min(chunk, args.samples - done)
        g = rng.standard_normal((m, syk.D))
        a = syk.eigvals(g) / np.linalg.norm(g, axis=1)[:, None]  # spectrum of sum_S ghat_S Gamma_S, |ghat| = 1
        a2 = a ** 2
        mom[it % args.batches] += [np.sum(np.mean(a2 ** l, axis=1)) for l in range(1, 5)]
        lam = (syk.sigma * np.abs(a).ravel()[:, None] * r_q[None, :]).ravel()
        w = np.broadcast_to(w_q[None, :] / syk.N, (a.size, args.Q)).ravel()
        acc.add(it % args.batches, lam, w, m)
        done += m
        it += 1
        now = time.time()
        if now - last > args.checkpoint:
            save(args.out, syk, acc, mom, args, now - t0)
            print(f"  {done} samples, {now - t0:.0f}s", flush=True)
            last = now
        if args.time_limit and now - t0 > args.time_limit:
            break
    save(args.out, syk, acc, mom, args, time.time() - t0)
    mu2 = mom[:, 1].sum() / done
    D = syk.D
    m4 = syk.sigma ** 4 * D * (D + 2) * mu2
    m4_exact = syk.V ** 2 * (2 + syk.q_exact())
    print(f"done {done} samples in {time.time() - t0:.1f}s; E tr H^2: {syk.sigma ** 2 * D * mom[:, 0].sum() / done:.8g} "
          f"(exact {syk.V:.8g}); E tr H^4: {m4:.8g} (exact {m4_exact:.8g}, rel diff {m4 / m4_exact - 1:+.2e})", flush=True)


if __name__ == "__main__":
    main()
