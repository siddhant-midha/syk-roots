"""Majorana operators, the SYK Hamiltonian and its exact diagonalization.

n Majoranas on m = n/2 qubits (Jordan-Wigner): gamma_{2j} = Z_0..Z_{j-1} X_j and gamma_{2j+1} = Z_0..Z_{j-1} Y_j.
A Pauli string is stored as (x, z, p), meaning i^p X^x Z^z (bit q of x / z acts on qubit q), so that

    i^p X^x Z^z |b> = i^p (-1)^{|z & b|} |b ^ x>.

The Hamiltonian follows paper 1 (Basu-Midha, Sec. 5):

    H_SYK = sigma_{n,k} sum_{|S|=k} g_S Gamma_S,   Gamma_S = i^{k(k-1)/2} gamma_{i_1} ... gamma_{i_k},
    sigma_{n,k}^2 = (k-1)! / (2 k n^{k-1}),

with g_S i.i.d. N(0,1) and tr the normalized trace, so V_n = E tr H^2 = sigma^2 C(n,k).
For even k every Gamma_S commutes with the fermion parity Z^{all}, so H is block diagonal in the two parity sectors.
"""
import itertools
import math

import numpy as np


def popcount(a):
    return np.bitwise_count(np.asarray(a, dtype=np.uint64)).astype(np.int64)


def majorana_paulis(n):
    """Pauli strings (x, z, p) of gamma_0 ... gamma_{n-1}."""
    j = np.arange(n // 2, dtype=np.int64)
    x = np.repeat(1 << j, 2)
    zlow = (1 << j) - 1
    z = np.stack([zlow, zlow | (1 << j)], axis=1).ravel()
    p = np.tile(np.array([0, 1], dtype=np.int64), n // 2)  # Y = i X Z
    return x, z, p


def gamma_strings(n, k):
    """All S in lexicographic order and the Pauli strings (x, z, p) of Gamma_S."""
    gx, gz, gp = majorana_paulis(n)
    S = np.array(list(itertools.combinations(range(n), k)), dtype=np.int64).reshape(-1, k)
    x = np.zeros(len(S), np.int64)
    z = np.zeros(len(S), np.int64)
    p = np.full(len(S), k * (k - 1) // 2, np.int64)
    for c in range(k):
        xi, zi, pi = gx[S[:, c]], gz[S[:, c]], gp[S[:, c]]
        p += pi + 2 * popcount(z & xi)  # (X^x Z^z)(X^xi Z^zi) = (-1)^{|z & xi|} X^{x^xi} Z^{z^zi}
        x ^= xi
        z ^= zi
    p %= 4
    assert np.all((p - popcount(x & z)) % 2 == 0), "Gamma_S is not Hermitian"
    return S, x, z, p


def pauli_matrix(x, z, p, N):
    """Dense matrix of i^p X^x Z^z on N = 2^m basis states (for tests)."""
    b = np.arange(N)
    M = np.zeros((N, N), complex)
    M[b ^ x, b] = (1j) ** p * (-1.0) ** popcount(z & b)
    return M


def wht(a):
    """Unnormalized Walsh-Hadamard transform over the last axis, in place: a[..., b] <- sum_z a[..., z] (-1)^{|z & b|}."""
    shape = a.shape
    N = shape[-1]
    h = 1
    while h < N:
        v = a.reshape(-1, N // (2 * h), 2, h)
        lo, hi = v[:, :, 0], v[:, :, 1]
        diff = lo - hi
        lo += hi
        hi[...] = diff
        h *= 2
    return a.reshape(shape)


class SYK:
    """Even-k SYK model on n Majoranas: builds the two parity blocks of sum_S g_S Gamma_S for batches of g."""

    def __init__(self, n, k):
        assert n % 2 == 0 and k % 2 == 0 and 2 <= k <= n
        self.n, self.k = n, k
        self.N = 1 << (n // 2)
        self.S, self.x, self.z, self.p = gamma_strings(n, k)
        self.D = len(self.S)
        self.phase = (1j) ** self.p
        self.sigma = math.sqrt(math.factorial(k - 1) / (2 * k * n ** (k - 1)))
        self.V = self.sigma ** 2 * self.D
        b = np.arange(self.N)
        par = popcount(b) % 2
        # H[b ^ x, b] = W[x, b] with W = WHT over z of the coefficients C[x, z]. Only even-weight x occur (k even);
        # row(x) indexes them. Block s: H_s[i, j] = W[row(B_i ^ B_j), B_j].
        row = np.full(self.N, -1, np.int64)
        row[par == 0] = np.arange(self.N // 2)
        self.row_x = row[self.x]
        self.flat_idx = []
        for s in (0, 1):
            B = b[par == s]
            self.flat_idx.append(row[B[:, None] ^ B[None, :]] * self.N + B[None, :])

    def blocks(self, g):
        """Parity blocks of sum_S g_S Gamma_S for g of shape (M, D): list of two (M, N/2, N/2) Hermitian arrays."""
        g = np.atleast_2d(g)
        C = np.zeros((g.shape[0], self.N // 2, self.N), complex)
        C[:, self.row_x, self.z] = g * self.phase
        W = wht(C).reshape(g.shape[0], -1)
        return [W[:, idx] for idx in self.flat_idx]

    def dense(self, g):
        """Full 2^{n/2} x 2^{n/2} matrix of sum_S g_S Gamma_S (for tests, small n only)."""
        return sum(gs * pauli_matrix(x, z, p, self.N) for gs, x, z, p in zip(g, self.x, self.z, self.p))

    def eigvals(self, g):
        """Eigenvalues of sum_S g_S Gamma_S (sigma not included) for g of shape (M, D): array (M, N)."""
        return np.concatenate([np.linalg.eigvalsh(Hs) for Hs in self.blocks(g)], axis=-1)

    def q_exact(self):
        """Average commutation sign E_{S,T} (-1)^{|S cap T|} at finite n (E tr H^4 = V^2 (2 + q_exact))."""
        n, k = self.n, self.k
        return sum((-1) ** j * math.comb(k, j) * math.comb(n - k, k - j) for j in range(k + 1)) / math.comb(n, k)
