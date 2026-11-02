"""Domain-coloring plots of the annealed SYK partition function Z_n(z) = E tr exp(z H_SYK) and its zeros.

Domain coloring: hue = arg Z, brightness = |Z| (gamma-corrected log scale), zeros = white dots. Additionally: open white circles = zeros
+-i t_j(q_n)/sqrt(V_n) of the q-oscillator comparator Q_n(z) = L_{q_n}(sqrt(V_n) z) of paper 1 (Lemma 5.2);
grey veil = region where the Monte Carlo (jackknife) noise exceeds 40% of |Z|; small grey rings = zeros there.

    python plot_roots.py results/ --out figs/
"""
import argparse
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import hsv_to_rgb  # noqa: E402
from scipy.interpolate import RegularGridInterpolator  # noqa: E402

NOISE_VEIL = 0.4
plt.rcParams.update({"font.size": 9, "mathtext.fontset": "cm", "axes.linewidth": 0.6})


def load(prefix):
    with open(prefix + ".json") as f:
        res = json.load(f)
    return res, np.load(prefix + ".npz")


def domain_rgb(L, phase, gamma=0.6, lo_pct=2, hi_pct=99.5):
    hue = (phase / (2 * np.pi)) % 1.0
    lo, hi = np.percentile(L, [lo_pct, hi_pct])
    val = np.clip((L - lo) / (hi - lo), 0, 1) ** gamma
    return hsv_to_rgb(np.stack([hue, np.ones_like(hue), val], -1))


def noise_mask(grid, x, y):
    """rel noise (stored on the first quadrant, z units) evaluated on a plotting grid (z units) by symmetry."""
    f = RegularGridInterpolator((grid["rel_y"], grid["rel_x"]), grid["rel"], bounds_error=False, fill_value=None)
    X, Y = np.meshgrid(np.abs(x), np.abs(y))
    return f(np.stack([Y, X], -1))


def panel(ax, res, grid, ms=9, qosc=True, xlim=None, ylim=None, label=True, common=False):
    """Domain coloring of one (n, k). common=False: analysis window in z; common=True: common window in t = sqrt(V) z."""
    tag, sc = ("c_", res["sqrt_V"]) if common else ("", 1.0)
    x, y = grid[tag + "x"], grid[tag + "y"]
    rgb = domain_rgb(grid[tag + "logabs"].astype(float), grid[tag + "phase"].astype(float))
    rel = noise_mask(grid, x, y)
    veil = (rel > NOISE_VEIL)[..., None]
    rgb = np.where(veil, 0.35 * rgb + 0.65 * 0.55, rgb)
    x, y = x * sc, y * sc
    ax.imshow(rgb, origin="lower", extent=[x[0], x[-1], y[0], y[-1]], aspect="auto", interpolation="bilinear")
    zs = res["zeros"]
    for sgn_re in (1, -1):
        for sgn_im in (1, -1):
            good = [(sgn_re * z["z_re"] * sc, sgn_im * z["z_im"] * sc) for z in zs if z["resolved"]]
            bad = [(sgn_re * z["z_re"] * sc, sgn_im * z["z_im"] * sc) for z in zs if not z["resolved"]]
            if good:
                g = np.array(good)
                ax.scatter(g[:, 0], g[:, 1], s=ms, c="white", linewidths=0, zorder=4)
            if bad:
                b = np.array(bad)
                ax.scatter(b[:, 0], b[:, 1], s=ms * 0.8, facecolors="none", edgecolors="0.8", linewidths=0.5, zorder=3)
    if qosc:
        tq = np.array(res["qosc_t"]) / res["sqrt_V"] * sc
        tq = tq[tq <= y[-1]]
        ax.scatter(np.zeros(2 * len(tq)), np.concatenate([tq, -tq]), s=ms * 2.2, facecolors="none", edgecolors="white",
                   linewidths=0.35, alpha=0.7, zorder=3)
    ax.set_xlim(xlim or (x[0], x[-1]))
    ax.set_ylim(ylim or (y[0], y[-1]))
    if label:
        ax.set_xlabel(r"$\Re\, (\sqrt{V_n} z)$" if common else r"$\Re\, z$")
        ax.set_ylabel(r"$\Im\, (\sqrt{V_n} z)$" if common else r"$\Im\, z$")


def title(res):
    return rf"$n={res['n']},\ k={res['k']}$  ($k^2/n={res['alpha']:.3g}$)"


def fig_single(prefix, out, zoom=None):
    res, grid = load(prefix)
    ncol = 2 if zoom else 1
    fig, axs = plt.subplots(1, ncol, figsize=(3.4 * ncol + 0.2, 4.6), squeeze=False)
    panel(axs[0, 0], res, grid, ms=10)
    axs[0, 0].set_title(title(res))
    if zoom:
        panel(axs[0, 1], res, grid, ms=22, xlim=zoom[0], ylim=zoom[1])
        axs[0, 1].set_title("zoom")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def fig_remark_exact(ed_prefix, exact_json, out):
    """n = 8, k = 4: ED/Monte Carlo (left) next to the exact traceless-GOE(8) partition function (right), full window
    (top) and zoom on the first off-axis zero (bottom)."""
    from goe import ExactZ84
    res, grid = load(ed_prefix)
    exact = json.load(open(exact_json))
    x, y = grid["x"], grid["y"]
    Zx = ExactZ84(res["V"])(x[None, :] + 1j * y[:, None])
    grid_x = dict(x=x, y=y, logabs=np.log(np.abs(Zx) + 1e-300), phase=np.angle(Zx), rel=np.zeros_like(grid["rel"]),
                  rel_x=grid["rel_x"], rel_y=grid["rel_y"])
    res_x = dict(res, zeros=[dict(z, resolved=True) for z in exact["zeros"]])
    z0 = min((z for z in exact["zeros"] if not z["on_axis"]), key=lambda z: z["z_re"] ** 2 + z["z_im"] ** 2)
    zoom = ((-0.15 * z0["z_re"], 1.6 * z0["z_re"]), (0, 1.6 * z0["z_im"]))
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 9.6))
    mant, expo = f"{res['samples']:.1e}".split("e")
    for col, (r, g, name) in enumerate([(res, grid, rf"ED, ${mant}\times10^{{{int(expo)}}}$ disorder samples"),
                                        (res_x, grid_x, "exact: traceless GOE(8)")]):
        panel(axs[0, col], r, g, ms=10)
        axs[0, col].set_title(rf"$n=8,\ k=4$: {name}", fontsize=9)
        panel(axs[1, col], r, g, ms=24, xlim=zoom[0], ylim=zoom[1])
        axs[1, col].set_title("zoom", fontsize=9)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_grid(results, out, ks, ns):
    """All (n, k) in the common window of t = sqrt(V_n) z (unit-variance Hamiltonian), so panels compare directly."""
    fig, axs = plt.subplots(len(ks), len(ns), figsize=(1.45 * len(ns) + 0.5, 2.75 * len(ks) + 0.4), squeeze=False,
                            sharex=True, sharey=True)
    for i, k in enumerate(ks):
        for j, n in enumerate(ns):
            ax = axs[i, j]
            prefix = results.get((n, k))
            if prefix is None:
                ax.set_axis_off()
                continue
            res, grid = load(prefix)
            panel(ax, res, grid, ms=2.5, label=False, common=True)
            ax.set_title(rf"$n={n},\,k={k}$", fontsize=8, pad=2)
            ax.tick_params(labelsize=6, length=2, pad=1)
        first = next(j for j, n in enumerate(ns) if (n, k) in results)
        axs[i, first].set_ylabel(r"$\Im\,(\sqrt{V_n}\,z)$", fontsize=8, labelpad=1)
        axs[i, first].tick_params(labelleft=True)
    for ax in axs[-1]:
        ax.set_xlabel(r"$\Re\,(\sqrt{V_n}\,z)$", fontsize=8, labelpad=1)
    fig.tight_layout(pad=0.3, w_pad=0.15, h_pad=0.5)
    fig.savefig(out, dpi=200)
    plt.close(fig)


def fig_onset(results, out):
    """|t| = |z| sqrt(V_n) of the first resolved off-axis zero against sqrt(2D), D = C(n, k)."""
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    colors = {2: "#1f77b4", 4: "#d62728", 6: "#2ca02c", 8: "#9467bd"}
    xs = []
    for (n, k), prefix in sorted(results.items()):
        res = json.load(open(prefix + ".json"))
        x = res["sqrt_2D"]
        xs.append(x)
        c = colors.get(k, "k")
        if res["onset_t"] is not None:
            ax.plot(x, res["onset_t"], "o", color=c, ms=5)
            ax.annotate(str(n), (x, res["onset_t"]), textcoords="offset points", xytext=(4, -3), fontsize=6, color=c)
        else:
            ax_t = [z["t_im"] for z in res["zeros"] if z["on_axis"] and z["resolved"]]
            if ax_t:
                ax.plot(x, max(ax_t), "^", color=c, ms=5, mfc="none")
    for k, c in colors.items():
        if any(kk == k for _, kk in results):
            ax.plot([], [], "o", color=c, label=f"k={k}")
    ax.plot([], [], "^", color="k", mfc="none", label="all resolved zeros on axis\n(largest one shown)")
    lim = [0.8 * min(xs), 1.25 * max(xs)]
    ax.plot(lim, lim, "k:", lw=0.8, label=r"$|t| = \sqrt{2D}$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\sqrt{2D}$,  $D=\binom{n}{k}$")
    ax.set_ylabel(r"first off-axis zero  $|z|\sqrt{V_n}$")
    ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def axis_run(res):
    """Resolved imaginary-axis zeros below the first resolved off-axis zero (the on-axis regime)."""
    lim = res["onset_t"] if res["onset_t"] is not None else float("inf")
    return np.array([z["t_im"] for z in res["zeros"] if z["on_axis"] and z["resolved"] and z["t_im"] < lim])


def fig_axis_compare(results, out, ks):
    """On-axis zeros t_j of the annealed Z_n (t = sqrt(V_n) z) against the q-oscillator zeros t_j(q), for the
    paper's q_n = exp(-2k^2/n) (top) and the exact finite-n average commutation sign q~ (bottom)."""
    fig, axs = plt.subplots(2, len(ks), figsize=(2.6 * len(ks) + 0.5, 4.8), squeeze=False, sharex=True, sharey="row")
    cmap = plt.get_cmap("viridis")
    for col, k in enumerate(ks):
        ns = sorted(n for n, kk in results if kk == k)
        for i, n in enumerate(ns):
            res = json.load(open(results[(n, k)] + ".json"))
            t = axis_run(res)
            for row, key in enumerate(["qosc_t", "qosc_exact_t"]):
                tq = np.array(res[key])
                m = min(len(t), len(tq))
                if m:
                    axs[row, col].plot(np.arange(1, m + 1), t[:m] / tq[:m], ".-", color=cmap(i / max(1, len(ns) - 1)),
                                       ms=3, lw=0.8, label=f"n={n}")
        for row in range(2):
            axs[row, col].axhline(1, color="k", lw=0.6, ls=":")
            axs[row, col].set_xscale("log")
        axs[0, col].set_title(f"k={k}", fontsize=9)
        axs[1, col].set_xlabel("j")
        axs[0, col].legend(fontsize=6, frameon=False, ncol=2)
    axs[0, 0].set_ylabel(r"$t_j^{\rm ED} / t_j(q_n)$,  $q_n = e^{-2k^2/n}$")
    axs[1, 0].set_ylabel(r"$t_j^{\rm ED} / t_j(\tilde q)$,  $\tilde q = \mathbb{E}\,(-1)^{|S\cap T|}$")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("resdir")
    ap.add_argument("--out", default="figs")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    results = {}
    for path in glob.glob(os.path.join(args.resdir, "n*_k*.json")):
        res = json.load(open(path))
        results[(res["n"], res["k"])] = path[:-5]
    ns = sorted({n for n, _ in results})
    ks = sorted({k for _, k in results})
    for (n, k), prefix in sorted(results.items()):
        fig_single(prefix, os.path.join(args.out, f"zeros_n{n}_k{k}.png"))
    if (8, 4) in results:  # Remark 1.5: zoom on the first resolved off-axis zero
        res = json.load(open(results[(8, 4)] + ".json"))
        z0 = min((z for z in res["zeros"] if z["resolved"] and not z["on_axis"]), key=lambda z: z["t_re"] ** 2 + z["t_im"] ** 2)
        fig_single(results[(8, 4)], os.path.join(args.out, "remark_n8_k4.png"),
                   zoom=((-0.15 * z0["z_re"], 1.6 * z0["z_re"]), (0, 1.6 * z0["z_im"])))
        exact = os.path.join(args.resdir, "exact_n8_k4.json")
        if os.path.exists(exact):
            fig_remark_exact(results[(8, 4)], exact, os.path.join(args.out, "remark_n8_k4_ed_vs_exact.png"))
    fig_grid(results, os.path.join(args.out, "zeros_grid.png"), ks, ns)
    for k in ks:
        fig_grid(results, os.path.join(args.out, f"zeros_k{k}.png"), [k], ns)
    fig_onset(results, os.path.join(args.out, "onset.png"))
    fig_axis_compare(results, os.path.join(args.out, "axis_vs_qosc.png"), ks)


if __name__ == "__main__":
    main()
