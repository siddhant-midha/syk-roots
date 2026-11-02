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

RESULTS_PLACEHOLDER

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
```

`della/` holds the Slurm scripts used on Princeton's Della cluster.
