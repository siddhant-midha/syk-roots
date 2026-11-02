"""Tables from the analyze.py outputs: results/summary.csv, results/summary.md (all (n, k)) and
results/zeros_n{n}_k{k}.md (resolved zeros with jackknife errors).

    python summarize.py results [--readme README.md]   (also refreshes the README figure gallery)
"""
import argparse
import csv
import glob
import json
import os


def zero_table(res):
    lines = [f"### n = {res['n']}, k = {res['k']}  ({res['samples']:.3g} disorder samples)", "",
             "Resolved zeros in the first quadrant (each comes with its images -z, +-conj z). "
             "t = sqrt(V_n) z; errors are jackknife (16 groups).", "",
             "| z | t = sqrt(V_n) z | err Re t | err Im t | on axis |", "|---|---|---|---|---|"]
    for z in res["zeros"]:
        if not z["resolved"]:
            continue
        lines.append(f"| {z['z_re']:.4f} + {z['z_im']:.4f}i | {z['t_re']:.4f} + {z['t_im']:.4f}i | "
                     f"{z['err_re_t']:.2g} | {z['err_im_t']:.2g} | {'yes' if z['on_axis'] else '**no**'} |")
    return "\n".join(lines) + "\n"


def readme_gallery(readme, figdir="figs"):
    """Thumbnail table of figs/zeros_n{n}_k{k}.png between the <!-- gallery --> markers of the README."""
    figs = {}
    for path in glob.glob(os.path.join(os.path.dirname(readme) or ".", figdir, "zeros_n*_k*.png")):
        n, k = os.path.basename(path)[7:-4].split("_k")
        figs[(int(n), int(k))] = f"{figdir}/{os.path.basename(path)}"
    ns = sorted({n for n, _ in figs})
    ks = sorted({k for _, k in figs})
    html = ["<table>", "<tr><th></th>" + "".join(f"<th>n={n}</th>" for n in ns) + "</tr>"]
    for k in ks:
        cells = "".join(f'<td><a href="{figs[(n, k)]}"><img src="{figs[(n, k)]}" width="105"></a></td>'
                        if (n, k) in figs else "<td></td>" for n in ns)
        html.append(f"<tr><th>k={k}</th>{cells}</tr>")
    html.append("</table>")
    text = open(readme).read()
    a, b = text.index("<!-- gallery -->"), text.index("<!-- /gallery -->")
    with open(readme, "w") as f:
        f.write(text[:a] + "<!-- gallery -->\n" + "\n".join(html) + "\n" + text[b:])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("resdir")
    ap.add_argument("--readme", default="")
    args = ap.parse_args()
    rows = []
    for path in sorted(glob.glob(os.path.join(args.resdir, "n*_k*.json"))):
        res = json.load(open(path))
        zs = res["zeros"]
        axis = [z for z in zs if z["on_axis"] and z["resolved"]]
        off = [z for z in zs if not z["on_axis"] and z["resolved"]]
        first_off = min(off, key=lambda z: z["t_re"] ** 2 + z["t_im"] ** 2) if off else None
        lim = res["onset_t"] if res["onset_t"] is not None else float("inf")
        rows.append(dict(
            n=res["n"], k=res["k"], D=res["D"], sqrt_2D=round(res["sqrt_2D"], 2), alpha=round(res["alpha"], 3),
            V_n=round(res["V"], 5), q_n=round(res["q_n"], 4), q_exact=round(res["q_exact"], 4),
            samples=res["samples"], window_re_t=round(res["window_t"][0], 1), window_im_t=round(res["window_t"][1], 1),
            zeros=len(zs), resolved=sum(z["resolved"] for z in zs), axis_resolved=len(axis), offaxis_resolved=len(off),
            axis_before_onset=sum(1 for z in axis if z["t_im"] < lim),
            onset_abs_t=round(res["onset_t"], 3) if res["onset_t"] is not None else "",
            first_offaxis_t=f"{first_off['t_re']:.3f}+{first_off['t_im']:.3f}i" if first_off else "",
            first_offaxis_err=f"({first_off['err_re_t']:.2g}, {first_off['err_im_t']:.2g})" if first_off else "",
            largest_resolved_axis_t=round(max(z["t_im"] for z in axis), 3) if axis else "",
            t1=round(axis[0]["t_im"], 4) if axis else "", t1_qn=round(res["qosc_t"][0], 4),
            t1_qexact=round(res["qosc_exact_t"][0], 4)))
        with open(os.path.join(args.resdir, f"zeros_n{res['n']}_k{res['k']}.md"), "w") as f:
            f.write(zero_table(res))
    rows.sort(key=lambda r: (r["k"], r["n"]))
    with open(os.path.join(args.resdir, "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    cols = [("n", "n"), ("k", "k"), ("D", "D"), ("sqrt_2D", "sqrt(2D)"), ("samples", "samples"),
            ("axis_before_onset", "axis zeros before onset"), ("onset_abs_t", "first off-axis abs(t)"),
            ("first_offaxis_t", "first off-axis t"), ("first_offaxis_err", "err (Re, Im)"),
            ("largest_resolved_axis_t", "largest resolved axis t"), ("window_im_t", "window Im t")]
    md = ["| " + " | ".join(c[1] for c in cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        md.append("| " + " | ".join(f"{r[c[0]]:.3g}" if c[0] == "samples" else str(r[c[0]]) for c in cols) + " |")
    with open(os.path.join(args.resdir, "summary.md"), "w") as f:
        f.write("\n".join(md) + "\n")
    print("\n".join(md))
    if args.readme:
        readme_gallery(args.readme)


if __name__ == "__main__":
    main()
