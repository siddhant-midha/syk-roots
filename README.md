# syk-roots

Zeros of the **annealed SYK partition function** on the complex plane, by exact diagonalization:

$$Z_n(z) = \mathbb{E}_{\boldsymbol g}\, \mathrm{tr}\, e^{z H_{\rm SYK}(\boldsymbol g)}, \qquad
H_{\rm SYK} = \sigma_{n,k} \sum_{|S|=k} g_S \Gamma_S,\quad \sigma_{n,k}^2 = \frac{(k-1)!}{2k\,n^{k-1}},$$

with $g_S$ i.i.d. $\mathcal N(0,1)$, $\Gamma_S = i^{k(k-1)/2}\gamma_{i_1}\cdots\gamma_{i_k}$ and $\mathrm{tr}$ the
normalized trace. These are the conventions of Basu–Midha, *Lee-Yang theorems for q-deformed oscillators and
doubly-scaled SYK* (Sec. 5), whose Theorem 1.4 says that all zeros of $Z_n$ go to the imaginary axis as
$n,k\to\infty$ with $k^2/n\to\alpha$, and whose **Remark 1.5** says that this fails at finite size, e.g. for
$n=8$, $k=4$ the annealed partition function has zeros with non-zero real part. This repo checks Remark 1.5
numerically and maps the zeros for a range of $(n,k)$. The plots use the style of Fig. 3 of *SYK thermal
expectations are classically easy at any temperature*.

Below, $t=\sqrt{V_n}\,z$ with $V_n=\mathbb E\,\mathrm{tr}H_{\rm SYK}^2=\sigma_{n,k}^2\binom nk$. This is the variable in
which the paper's comparator $Q_n(z)=\mathcal L_{q_n}(\sqrt{V_n}z)$ has its zeros at $\pm i\,t_j(q_n)$. Zeros come in
sets $\{\pm z,\pm\bar z\}$; only the first quadrant is listed.

## Results

### Remark 1.5 (n = 8, k = 4): confirmed

![remark](figs/remark_n8_k4.png)

$2.6\times10^9$ disorder samples ($\sqrt{V_8}=0.3202$). The first seven zero pairs lie exactly on the imaginary axis.
The plot shows a subsequent pair to be off it:

| t = √V z | z | jackknife error (Re t, Im t) |
|---|---|---|
| 1.9582 i | 6.115 i | (0, 3e-6) |
| 3.5648 i | 11.132 i | (0, 1e-5) |
| 5.1868 i | 16.198 i | (0, 2e-5) |
| 6.8728 i | 21.463 i | (0, 5e-5) |
| 8.6758 i | 27.093 i | (0, 1e-4) |
| **8.553 + 10.548 i** | **26.71 + 32.94 i** | **(0.017, 0.026)** |
| 10.6807 i | 33.355 i | (0, 1e-4) |
| 13.0868 i | 40.868 i | (0, 1e-4) |

So $Z_8$ has the quadruplet $z\approx\pm26.71\pm32.94\,i$, whose real part is about 500 standard errors from zero. Two
disjoint halves of the samples give the same zero. Beyond $|t|\approx13$ the remaining zeros fan out into the plane
(open grey rings: there the Monte Carlo noise is too large to pin them down individually). The zoom shows the
off-axis zero as a phase singularity, where all colours meet, away from the axis. Full table:
[results/zeros_n8_k4.md](results/zeros_n8_k4.md).

### Range of (n, k)

All panels use the common window $|\Re t|\le18$, $|\Im t|\le36$. White dots are resolved zeros, grey rings are
unresolved ones, open white circles are the zeros of the paper's $Q_n$, and the grey veil marks where the noise
exceeds 40% of $|Z|$. Single-$(n,k)$ figures in the paper's $z$ over the full analysis window follow in the
gallery below; tables are in `results/`.

![Zeros of the annealed SYK partition function for k = 2, 4, 6, 8 (rows) and n = 8 ... 18 (columns)](figs/zeros_grid.png)

| k | n | first resolved off-axis zero (t) | resolved zeros on the axis below it |
|---|---|---|---|
| 4 | 8 | 8.553 + 10.548 i | 7 |
| 4 | 10 | none resolved (axis zeros resolved up to t = 28.2) | 15 |
| 4 | 12 | 5.475 + 29.034 i | 15 |
| 4 | 14 | 1.81 + 62.65 i | 33 |
| 4 | 16 | none (all resolved zeros on the axis, up to t = 56.8) | 35 |
| 4 | 18 | none (up to t = 77.2) | 49 |
| 6 | 10 / 14 / 18 | none (up to t = 28.2 / 117.5 / 119.0) | 15 / 62 / 75 |
| 6 | 12 / 16 | 0.462 + 13.870 i / 0.086 + 27.864 i (a row of near-axis pairs) | 7 / 16 |
| 8 | 14 / 16 / 18 | none (up to t = 117.5 / 107.2 / 119.3) | 62 / 63 / 75 |
| 2 | 8 – 18 | between 0.25 + 7.96 i and 1.34 + 11.19 i at every n | 2 – 7 |

[results/summary.md](results/summary.md) has all pairs.

#### Each (n, k) in the paper's variable z (click to enlarge)

<!-- gallery -->
<table>
<tr><th></th><th>n=8</th><th>n=10</th><th>n=12</th><th>n=14</th><th>n=16</th><th>n=18</th></tr>
<tr><th>k=2</th><td><a href="figs/zeros_n8_k2.png"><img src="figs/zeros_n8_k2.png" width="105"></a></td><td><a href="figs/zeros_n10_k2.png"><img src="figs/zeros_n10_k2.png" width="105"></a></td><td><a href="figs/zeros_n12_k2.png"><img src="figs/zeros_n12_k2.png" width="105"></a></td><td><a href="figs/zeros_n14_k2.png"><img src="figs/zeros_n14_k2.png" width="105"></a></td><td><a href="figs/zeros_n16_k2.png"><img src="figs/zeros_n16_k2.png" width="105"></a></td><td><a href="figs/zeros_n18_k2.png"><img src="figs/zeros_n18_k2.png" width="105"></a></td></tr>
<tr><th>k=4</th><td><a href="figs/zeros_n8_k4.png"><img src="figs/zeros_n8_k4.png" width="105"></a></td><td><a href="figs/zeros_n10_k4.png"><img src="figs/zeros_n10_k4.png" width="105"></a></td><td><a href="figs/zeros_n12_k4.png"><img src="figs/zeros_n12_k4.png" width="105"></a></td><td><a href="figs/zeros_n14_k4.png"><img src="figs/zeros_n14_k4.png" width="105"></a></td><td><a href="figs/zeros_n16_k4.png"><img src="figs/zeros_n16_k4.png" width="105"></a></td><td><a href="figs/zeros_n18_k4.png"><img src="figs/zeros_n18_k4.png" width="105"></a></td></tr>
<tr><th>k=6</th><td><a href="figs/zeros_n8_k6.png"><img src="figs/zeros_n8_k6.png" width="105"></a></td><td><a href="figs/zeros_n10_k6.png"><img src="figs/zeros_n10_k6.png" width="105"></a></td><td><a href="figs/zeros_n12_k6.png"><img src="figs/zeros_n12_k6.png" width="105"></a></td><td><a href="figs/zeros_n14_k6.png"><img src="figs/zeros_n14_k6.png" width="105"></a></td><td><a href="figs/zeros_n16_k6.png"><img src="figs/zeros_n16_k6.png" width="105"></a></td><td><a href="figs/zeros_n18_k6.png"><img src="figs/zeros_n18_k6.png" width="105"></a></td></tr>
<tr><th>k=8</th><td></td><td><a href="figs/zeros_n10_k8.png"><img src="figs/zeros_n10_k8.png" width="105"></a></td><td><a href="figs/zeros_n12_k8.png"><img src="figs/zeros_n12_k8.png" width="105"></a></td><td><a href="figs/zeros_n14_k8.png"><img src="figs/zeros_n14_k8.png" width="105"></a></td><td><a href="figs/zeros_n16_k8.png"><img src="figs/zeros_n16_k8.png" width="105"></a></td><td><a href="figs/zeros_n18_k8.png"><img src="figs/zeros_n18_k8.png" width="105"></a></td></tr>
</table>
<!-- /gallery -->

#### Where the zeros leave the axis

<img src="figs/onset.png" width="480" alt="first off-axis zero against sqrt(2D)">

* **Approach to the axis.** For k = 4, 6 and 8 the off-axis zeros move outward as n grows. By n = 16–18 every
  zero resolved in windows up to $|t|\approx 120$ lies on the imaginary axis, as Theorem 1.4 predicts for
  $n\to\infty$. For k = 4 the first off-axis zero roughly tracks $|t|\sim\sqrt{2D}$ with $D=\binom nk$ (dotted
  line). Open triangles mark (n,k) with no resolved off-axis zero, plotted at the largest resolved axis zero.
  A heuristic for that scale: integrating out the overall coupling scale $|\boldsymbol g|\sim\chi_D$ damps the
  oscillations of $Z_n(it)$ beyond $|t|\sim\sqrt{2D}$.
* **k = 2 (free fermions)** has off-axis zeros at $|t|\approx 6$–11 at every n up to 18. Here $q_n=e^{-8/n}\to1$, so
  no Lee–Yang limit is expected.
* **k = 6 at n = 12, 16** (n ≡ 0 mod 4) has a regular row of near-axis pairs, Re t ≈ 0.1–2.4, robust to $10^{-3}$.
  At the same (n,k) the fixed-$|\boldsymbol g|$ spectrum has an anomaly at E = 0 on the scale of one level spacing:
  the density there is 0.15× (n = 12) and 1.7× (n = 16) the smooth value, versus ≈1 at n = 10, 14, 18. This is
  consistent with the particle-hole structure of q ≡ 2 (mod 4) SYK at n ≡ 0 (mod 4). It is a finite-size effect,
  whose scale shrinks with the level spacing.
* **Duality check.** $\Gamma_S\Gamma_{[n]}\propto\Gamma_{S^c}$, so (n,k) and (n,n−k) have the same $Z_n$ in $t$.
  Independent runs reproduce this: (12,4) and (12,8) give $5.475+29.034i$ and $5.474+29.034i$; (8,2)/(8,6) and
  (10,2)/(10,8) agree as well.

#### Axis zeros against the q-oscillator

![t_j(ED) / t_j(q) for the paper's q_n (top) and the exact q~ (bottom)](figs/axis_vs_qosc.png)

* **q-oscillator comparison.** On-axis zeros $t_j$ against $t_j(q)$, maximum relative deviation over the first ≤10:

  | | k=4, n=12 | n=14 | n=16 | n=18 | k=6, n=18 | k=8, n=18 |
  |---|---|---|---|---|---|---|
  | paper's $q_n=e^{-2k^2/n}$ | 8.9% | 9.9% | 9.2% | 8.6% | 0.7% | 0.0% |
  | exact $\tilde q=\mathbb E(-1)^{\lvert S\cap T\rvert}$ | 3.5% | 2.4% | 0.8% | 0.4% | 0.2% | 0.0% |

  With $\tilde q$ (the value that matches $\mathbb E\,\mathrm{tr}H^4=V_n^2(2+\tilde q)$), $\mathcal L_{\tilde q}(\sqrt{V_n}z)$
  locates the finite-n axis zeros to under 1% from n ≈ 16. Unlike $e^{-2k^2/n}$, $\tilde q$ also respects the
  k ↔ n−k duality.

n = 20, 22, 24 (ED on 2¹⁰–2¹² dimensional blocks) are being added.

## Method

* **ED.** $n$ Majoranas on $n/2$ qubits (Jordan–Wigner). Each $\Gamma_S$ is a Pauli string $i^p X^x Z^z$. For one
  disorder sample, $\sum_S g_S\Gamma_S$ is built in $O(N^2\log N)$ ($N=2^{n/2}$): scatter the couplings into
  $C[x,z]$, apply a Walsh–Hadamard transform over $z$, and read $H[b\oplus x,b]=W[x,b]$. Even $k$ conserves
  fermion parity, so the two $N/2$ blocks are diagonalized separately (`sykmaj.py`).
* **Annealed average with the radius integrated out exactly.** Write $\boldsymbol g = r\,\hat{\boldsymbol g}$ with
  $r\sim\chi_D$, $D=\binom nk$, independent of the direction. For each sampled $\hat{\boldsymbol g}$ with spectrum
  $a_j$, $\mathbb E_r \cosh(z\sigma r a_j) = {}_1F_1(D/2;1/2;z^2\sigma^2a_j^2/2)$ is computed with a 16-node
  Gauss rule for $\chi_D$ (machine precision). $Z_n$ is even because $\boldsymbol g\to-\boldsymbol g$ flips the
  spectrum. The second moment $\mathbb E\,\mathrm{tr}H^2=V_n$ is then exact, and the fourth moment is checked
  against $V_n^2(2+\tilde q)$, $\tilde q = \binom nk^{-1}\sum_j(-1)^j\binom kj\binom{n-k}{k-j}$.
* **Storage.** The weighted eigenvalue atoms $\lambda\ge 0$ are binned with 16 Taylor moments per bin, so
  $Z_n(z)=\sum_{b,p}S_{bp}\,z^p\cosh^{(p)}(zc_b)/p!$ can be evaluated anywhere afterwards. The truncation error is below
  $10^{-18}$ in every window used. Samples are split into 16 groups to get jackknife errors (`annealed.py`).
* **Zeros.** On a grid in the first quadrant whose cells straddle the imaginary axis, cells with non-zero phase
  winding are polished by Newton. The count is checked strip by strip with an adaptive argument-principle
  contour integral, and mismatched strips are re-searched by recursive subdivision. Because $Z_n$ is even and
  real on $\mathbb R$, zeros come as $\{\pm z,\pm\bar z\}$, and a simple zero on $i\mathbb R$ cannot leave the axis
  under Monte Carlo noise. Off-axis zeros are therefore a genuine signal once their real part exceeds its error.
* **Resolved zeros.** For each zero we compute jackknife errors of $\Re,\Im$ and the local noise
  $\mathrm{std}_{\rm jk}Z/|Z|$ on a ring of radius $d/2$ around it. A zero counts as *resolved* when the noise is $<0.4$,
  the error is $<1/5$ of the distance $d$ to the nearest other zero, the zeros recomputed from two disjoint halves of
  the samples differ by $<d/5$, and, for off-axis zeros, $\Re z > 5\,\mathrm{err}$. The halves check guards against
  a failure the jackknife cannot see: at large $\Re z$ the average is weighted by $e^{\Re z\,\lambda}$, so rare
  large-eigenvalue samples can dominate it.
* **q-oscillator comparison.** The paper's finite-$n$ comparator $Q_n(z)=\mathcal L_{q_n}(\sqrt{V_n}z)$,
  $q_n=e^{-2k^2/n}$, has zeros $\pm i\,t_j(q_n)/\sqrt{V_n}$. Here $t_j(q)$ are computed from a truncated Jacobi
  matrix of $a_q+a_q^*$, i.e. the Gauss rule of the $q$-Gaussian (`qosc.py`). They are drawn as open circles.

## Reproduce

```
python tests/test_syk.py                                   # Majorana algebra, builder vs dense, estimator checks
python sample.py --n 8 --k 4 --samples 1000000 --seed 1 --out out/n8_k4/s1.npz
python analyze.py out/n8_k4/*.npz --out results/n8_k4      # zeros + jackknife errors -> json, plotting grid -> npz
python plot_roots.py results --out figs
python summarize.py results --readme README.md             # tables + refresh the figure gallery above
```

`della/` holds the Slurm scripts used on Princeton's Della cluster.
