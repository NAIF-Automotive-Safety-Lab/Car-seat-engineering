"""Wilson-Dirac fermion operator on a 3-D periodic lattice.

Physics summary
---------------
Naive (central-difference) lattice fermions suffer from the *doubling problem*:
in d space dimensions the spectrum of the free Dirac operator has 2^d
degenerate zero-modes at the 2^d corners of the Brillouin zone.  In the
Jackiw-Rebbi computation this means 8 spurious near-zero eigenvalues appear
alongside the single physical zero-mode, making the result uninterpretable.

Wilson's cure (1975) is to add the irrelevant dimension-5 operator

    S_W  =  −r/(2a) ψ̄ ∇² ψ         (r = Wilson parameter, 0 < r ≤ 1)

to the lattice action.  This generates an additive mass

    m_doubler(p)  =  r/a · Σ_k (1 − cos p_k a)

for every lattice mode.  At the zone centre p ≈ 0 this vanishes, leaving the
physical mode alone.  At zone corners, every cos p_k a = −1, so each doubler
picks up mass ~2rd/a → ∞ in the continuum limit a → 0.

The Wilson-Dirac operator in d=3 spatial dimensions is

    D_W  =  m(x)·I_4  +  Σ_{k=1}^{3} γ^k · [ψ_{n+k} − ψ_{n−k}]/(2a)
            −  r/(2a) · Σ_{k=1}^{3} [ψ_{n+k} − 2ψ_n + ψ_{n−k}]/a

with m(x) = M tanh(x/ξ) the kink background.

In matrix form for the full spinor vector (4·N³ complex components):

    D_W = m_diag·I  + γ^k ⊗ T_k^fwd  −  r/(2a) I_4 ⊗ L_k

where T_k^fwd is the off-diagonal lattice difference matrix along axis k and
L_k is the second-difference (Laplacian) matrix along axis k.

We build D_W as a complex sparse CSR matrix and then form H² = D_W^† D_W,
which is Hermitian positive-semidefinite.  Its lowest eigenvalue λ_0 satisfies
λ_0 = E_0² where E_0 is the fermion ground-state energy.  For a true
Jackiw-Rebbi zero-mode E_0 = 0, so λ_0 → 0 as a → 0 (N → ∞) at fixed L.

Conventions
-----------
- Natural units  ħ = c = 1
- Lattice spacing  a = L / N
- Spinor index α ∈ {0,1,2,3}
- Lattice site  n = (i,j,k),  flat index  i·N²+j·N+k
- Full DOF index  α·N³ + flat(n)
- Periodic boundary conditions (torus topology)
- Wilson parameter r = 1  (standard choice)

Gamma matrices (Dirac-Pauli representation)
--------------------------------------------
    γ^0 = diag(I₂, −I₂)
    γ^1 = [[0, σ_x], [−σ_x, 0]]
    γ^2 = [[0, σ_y], [−σ_y, 0]]
    γ^3 = [[0, σ_z], [−σ_z, 0]]

The spatial γ^k are anti-Hermitian (iγ^k is Hermitian), so the kinetic term
gives an anti-Hermitian contribution and D_W is not Hermitian—we must work
with H² = D_W^† D_W.

Performance notes
-----------------
For N = 50: matrix size = 4 × 125,000 = 500,000.  NNZ ≈ 14 × 500,000 = 7 M.
For N = 100: matrix size = 4,000,000.  NNZ ≈ 56 M.  Still fits in ~4 GB RAM.
For N = 200: matrix size = 32,000,000.  NNZ ≈ 450 M.  ~36 GB RAM; use on
  machines with at least 64 GB, or restrict eigsh to k=1 (only need lowest).

We use scipy.sparse.linalg.eigsh with which='SM' (smallest magnitude) on H²
and sigma=0 shift-invert for reliable near-zero eigenvalue extraction.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import NamedTuple

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


# ---------------------------------------------------------------------------
# Gamma matrices (Dirac-Pauli / standard representation)
# ---------------------------------------------------------------------------

def _gamma_matrices() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (γ^0, γ^1, γ^2, γ^3) as 4×4 complex arrays.

    Dirac-Pauli representation:
        γ^0 = diag(I, -I)
        γ^k = [[0, σ_k], [-σ_k, 0]]  for k = 1,2,3

    Returns:
        Tuple of four (4,4) complex128 arrays.
    """
    I2 = np.eye(2, dtype=np.complex128)
    Z2 = np.zeros((2, 2), dtype=np.complex128)
    sx = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sy = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    sz = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    g0 = np.block([[I2, Z2], [Z2, -I2]])
    g1 = np.block([[Z2, sx], [-sx, Z2]])
    g2 = np.block([[Z2, sy], [-sy, Z2]])
    g3 = np.block([[Z2, sz], [-sz, Z2]])
    return g0, g1, g2, g3


GAMMA: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] = _gamma_matrices()


# ---------------------------------------------------------------------------
# Lattice parameters
# ---------------------------------------------------------------------------

@dataclass
class WilsonParams:
    """Parameters controlling the Wilson-Dirac lattice computation.

    Attributes:
        N:    Lattice side length.  Total sites = N³.
        L:    Physical box length (in units where ξ = 1).
        M:    Kink mass parameter (asymptotic |m(x → ±∞)|).
        xi:   Kink width (tanh profile m(x) = M tanh(x/xi)).
        r:    Wilson parameter (standard value r = 1).
        kink_axis: Which axis carries the kink (0=x,1=y,2=z).  z is natural.
        boundary: Boundary condition for the spinor in the kink axis only.
                  'pbc'        -> ψ(x+L) = +ψ(x)  (periodic; kink+antikink)
                  'apbc'       -> ψ(x+L) = −ψ(x)  (anti-periodic; flips bdry hop)
                  'dirichlet'  -> ψ(0) = ψ(N−1) = 0 in the kink axis
                                  (open; ZERO boundary hops; SINGLE isolated kink)
                  Transverse axes are always PBC.

    Boundary-condition physics:

      PBC: Both the spinor and the periodic identification carry the
        m(x)=M tanh(x/ξ) profile around the box, but tanh is not periodic,
        so the boundary acts as an antikink.  Density carries weight at
        BOTH the centre kink (x=0) and the boundary (x=±L/2).  This is
        the standard textbook lattice setup; the price is doubled zero
        modes.

      APBC: A spinor sign-flip at the boundary shifts allowed momenta
        by π/L (removes the zero-momentum doubler) and gives a slightly
        different λ_0 than PBC, but the mass field is unchanged so the
        boundary antikink remains.  APBC is useful for matching finite-
        temperature fermion sums but does NOT isolate a single kink on
        its own.

      DIRICHLET: The boundary hop coefficients in the kink axis are zeroed
        out entirely — there is no wraparound and no boundary antikink.
        With m(x)=M tanh(x/ξ) the lattice carries exactly ONE Jackiw-Rebbi
        zero-mode at x=0.  This is the cleanest setup for testing the
        cosh⁻²(x/ξ) profile shape and the 1/N approach to λ_0=0; the
        cost is broken translation invariance (open boundary).
    """
    N: int = 50
    L: float = 16.0
    M: float = 1.0
    xi: float = 1.0
    r: float = 1.0
    kink_axis: int = 2   # kink along z-axis
    boundary: str = "pbc"   # 'pbc', 'apbc', or 'dirichlet'

    def __post_init__(self) -> None:
        if self.boundary not in ("pbc", "apbc", "dirichlet"):
            raise ValueError(
                f"WilsonParams.boundary must be 'pbc', 'apbc', or 'dirichlet'; "
                f"got {self.boundary!r}"
            )

    @property
    def a(self) -> float:
        """Lattice spacing a = L / N."""
        return self.L / self.N

    @property
    def n3(self) -> int:
        """Number of spatial lattice sites."""
        return self.N ** 3

    @property
    def ndof(self) -> int:
        """Total degrees of freedom = 4 spinor components × N³ sites."""
        return 4 * self.n3


# ---------------------------------------------------------------------------
# Kink background  m(x) = M tanh(x/ξ)
# ---------------------------------------------------------------------------

def kink_mass_field(p: WilsonParams) -> np.ndarray:
    """Construct the kink background mass profile m(x) = M tanh(x_k / ξ).

    The kink is centred at x_k = 0 (half-box, exploiting periodic boundary).
    We shift coordinates so that the kink centre lies at the box midpoint and
    the field runs from −M (left) to +M (right) along the kink axis.

    Args:
        p: Wilson lattice parameters.

    Returns:
        Real array of shape (N, N, N) with values in (−M, +M).
    """
    a = p.a
    coords = (np.arange(p.N) - p.N / 2.0) * a   # centred: [−L/2, L/2)

    # Build 3D mass field: kink along axis p.kink_axis, trivial transverse
    shape = [1, 1, 1]
    shape[p.kink_axis] = p.N
    m1d = (p.M * np.tanh(coords / p.xi)).reshape(shape)
    return np.broadcast_to(m1d, (p.N, p.N, p.N)).copy()


# ---------------------------------------------------------------------------
# Sparse Wilson-Dirac matrix construction
# ---------------------------------------------------------------------------

def _flat_index(i: int | np.ndarray, j: int | np.ndarray, k: int | np.ndarray, N: int) -> np.ndarray:
    """Convert (i,j,k) triple to flat site index with periodic wrapping.

    Args:
        i, j, k: Lattice coordinate arrays.
        N:        Lattice size.

    Returns:
        Flat index array of same shape.
    """
    i = np.asarray(i) % N
    j = np.asarray(j) % N
    k = np.asarray(k) % N
    return i * N * N + j * N + k


def build_wilson_dirac(p: WilsonParams) -> sp.csr_matrix:
    """Construct the Wilson-Dirac operator D_W as a sparse CSR matrix.

    The operator acts on the full spinor-lattice vector v of length 4·N³,
    ordered as:  v[α·N³ + flat(n)]  for spinor component α, site n.

    The Wilson-Dirac operator is:

        D_W = m(x) I₄  +  Σ_k γ^k (∂_k)_lattice  −  r/(2a) Σ_k (∇²_k)_lattice I₄

    where:
        (∂_k ψ)_n  =  [ψ_{n+k} − ψ_{n−k}] / (2a)       (central difference)
        (∇²_k ψ)_n =  [ψ_{n+k} − 2ψ_n + ψ_{n−k}] / a²  (second difference)

    The Wilson term  −r/(2a) ∇²_k  gives a positive mass ~r·d/a to doublers
    while leaving the zone-centre mode (physical) untouched.

    Matrix structure (block form, 4 spinor × N³ sites):
        Row α·N³+n, column β·N³+m:

        Diagonal (n=m, α=β):
            +m(x_n)                                  [mass term]
            +r/a · d                                  [Wilson diagonal, d=3 dims]

        Off-diagonal hopping (m = n±k̂):
            −r/(2a²) I_{αβ}                           [Wilson hopping, same spin]
            ±(γ^k)_{αβ} / (2a)                       [Dirac hopping, mix spin]

    Total NNZ ≈ (1 + 6 + 6) × 4·N³ = 52·N³ entries.

    Args:
        p: Wilson lattice parameters (includes kink background).

    Returns:
        Complex CSR matrix of shape (4·N³, 4·N³).
    """
    N = p.N
    n3 = p.n3
    ndof = p.ndof
    a = p.a
    r = p.r

    m_field = kink_mass_field(p)   # shape (N,N,N)

    g0, g1, g2, g3 = GAMMA
    gammas_spatial = [g1, g2, g3]  # γ^1, γ^2, γ^3

    # Pre-compute all (i,j,k) site indices as flat arrays for vectorised construction
    ii, jj, kk = np.meshgrid(
        np.arange(N), np.arange(N), np.arange(N), indexing='ij'
    )
    ii = ii.ravel()   # length N³
    jj = jj.ravel()
    kk = kk.ravel()
    site_flat = ii * N * N + jj * N + kk   # [0, N³)

    # Accumulate COO data
    rows_list: list[np.ndarray] = []
    cols_list: list[np.ndarray] = []
    data_list: list[np.ndarray] = []

    # ------------------------------------------------------------------
    # 1. Diagonal: mass term + Wilson diagonal  r·d/a · I
    # ------------------------------------------------------------------
    # Entry for (α, n) diagonal: m(x_n) + r·d/a
    m_vals = m_field[ii, jj, kk]   # shape (N³,)
    diag_val = m_vals + r * 3.0 / a   # scalar + Wilson diagonal contribution

    for alpha in range(4):
        row_idx = alpha * n3 + site_flat
        col_idx = alpha * n3 + site_flat
        rows_list.append(row_idx)
        cols_list.append(col_idx)
        data_list.append(diag_val.astype(np.complex128))

    # ------------------------------------------------------------------
    # 2. Hopping terms along each spatial axis k = 0,1,2
    # ------------------------------------------------------------------
    # For each spatial direction k, the neighbouring sites are:
    #   n+k̂:  displacement (+1) along axis k
    #   n−k̂:  displacement (−1) along axis k
    #
    # The contribution to D_W from axis k is:
    #   Dirac forward:   +γ^k_{αβ} / (2a)   × ψ_β(n+k̂)
    #   Dirac backward:  −γ^k_{αβ} / (2a)   × ψ_β(n−k̂)
    #   Wilson forward:  −r/(2a²) δ_{αβ}    × ψ_β(n+k̂)
    #   Wilson backward: −r/(2a²) δ_{αβ}    × ψ_β(n−k̂)

    axes_shifts: list[tuple[np.ndarray, np.ndarray, np.ndarray]] = []
    # For each axis, precompute the ±1 neighbour flat indices
    ip = (ii + 1) % N; im = (ii - 1) % N   # x-axis (axis 0)
    jp = (jj + 1) % N; jm = (jj - 1) % N   # y-axis (axis 1)
    kp = (kk + 1) % N; km = (kk - 1) % N   # z-axis (axis 2)

    # site_plus[k] = flat index of n+k̂
    site_plus = [
        ip * N * N + jj * N + kk,    # n + x̂
        ii * N * N + jp * N + kk,    # n + ŷ
        ii * N * N + jj * N + kp,    # n + ẑ
    ]
    site_minus = [
        im * N * N + jj * N + kk,    # n − x̂
        ii * N * N + jm * N + kk,    # n − ŷ
        ii * N * N + jj * N + km,    # n − ẑ
    ]

    # ------------------------------------------------------------------
    # Boundary-condition sign masks (per-axis).
    #
    # For PBC the hop  ψ(N−1) → ψ(0)  carries coefficient +1.
    # For APBC in the kink axis the same hop carries −1, while transverse
    # axes stay PBC.  A boundary-wrap is detected by:
    #     forward  hop wrapping when index value == N−1  (next-site is 0)
    #     backward hop wrapping when index value == 0    (prev-site is N−1)
    # The per-site arrays  fwd_sign[k_axis]  and  bwd_sign[k_axis]  hold the
    # multiplier; they are +1 everywhere except at the wrap sites for APBC.
    # ------------------------------------------------------------------
    fwd_sign = [np.ones(n3, dtype=np.complex128) for _ in range(3)]
    bwd_sign = [np.ones(n3, dtype=np.complex128) for _ in range(3)]
    if p.boundary in ("apbc", "dirichlet"):
        idx_arrs = [ii, jj, kk]
        idx_at_kink_axis = idx_arrs[p.kink_axis]
        # APBC: sign-flip the boundary-wrapping hops in the kink axis.
        # DIRICHLET: zero them entirely — no wraparound, no antikink.
        bdry_value = -1.0 if p.boundary == "apbc" else 0.0
        # Sites where the FORWARD hop wraps to the other side
        fwd_sign[p.kink_axis] = np.where(
            idx_at_kink_axis == N - 1, bdry_value, 1.0
        ).astype(np.complex128)
        # Sites where the BACKWARD hop wraps to the other side
        bwd_sign[p.kink_axis] = np.where(
            idx_at_kink_axis == 0, bdry_value, 1.0
        ).astype(np.complex128)

    wilson_hop = -r / (2.0 * a * a)           # Wilson hopping coefficient
    dirac_coeff_fwd = +1.0 / (2.0 * a)        # +γ^k / (2a)
    dirac_coeff_bwd = -1.0 / (2.0 * a)        # −γ^k / (2a)

    for k_axis, gk in enumerate(gammas_spatial):
        s_plus = site_plus[k_axis]
        s_minus = site_minus[k_axis]
        sign_fwd = fwd_sign[k_axis]
        sign_bwd = bwd_sign[k_axis]

        for alpha in range(4):
            row_alpha = alpha * n3 + site_flat   # row indices for this (α,n)

            for beta in range(4):
                # Wilson hopping: −r/(2a²) δ_{αβ}  (only diagonal in spinor)
                if alpha == beta:
                    # Forward neighbour (sign-flipped at wrap if APBC)
                    rows_list.append(row_alpha)
                    cols_list.append(beta * n3 + s_plus)
                    data_list.append(wilson_hop * sign_fwd)
                    # Backward neighbour
                    rows_list.append(row_alpha)
                    cols_list.append(beta * n3 + s_minus)
                    data_list.append(wilson_hop * sign_bwd)

                # Dirac hopping: ±γ^k_{αβ} / (2a)
                gk_ab = gk[alpha, beta]
                if abs(gk_ab) > 1e-14:
                    # Forward: +γ^k_{αβ} / (2a)
                    rows_list.append(row_alpha)
                    cols_list.append(beta * n3 + s_plus)
                    data_list.append(dirac_coeff_fwd * gk_ab * sign_fwd)
                    # Backward: −γ^k_{αβ} / (2a)
                    rows_list.append(row_alpha)
                    cols_list.append(beta * n3 + s_minus)
                    data_list.append(dirac_coeff_bwd * gk_ab * sign_bwd)

    # Assemble COO → CSR
    rows_arr = np.concatenate(rows_list)
    cols_arr = np.concatenate(cols_list)
    data_arr = np.concatenate(data_list)

    D_W = sp.coo_matrix(
        (data_arr, (rows_arr, cols_arr)),
        shape=(ndof, ndof),
        dtype=np.complex128,
    ).tocsr()

    return D_W


def build_H_squared(D_W: sp.csr_matrix) -> sp.csr_matrix:
    """Compute H² = D_W^† D_W (Hermitian positive-semidefinite).

    The Wilson-Dirac operator D_W is not Hermitian.  To get a Hermitian
    eigenproblem we form H² = D_W† D_W.  The eigenvalues of H² are E²,
    where E are the singular values of D_W (equivalently the eigenvalues of
    the Hermitian H = √(D_W† D_W)).  A zero-mode satisfies D_W ψ = 0, so
    H² ψ = 0 as well.

    Args:
        D_W: Wilson-Dirac operator as CSR matrix.

    Returns:
        H² = D_W^H D_W as CSR matrix (Hermitian by construction).
    """
    D_W_dag = D_W.conj().T
    return (D_W_dag @ D_W).tocsr()


# ---------------------------------------------------------------------------
# Eigenvalue solver
# ---------------------------------------------------------------------------

class ZeroModeResult(NamedTuple):
    """Results from the zero-mode eigenvalue extraction.

    Attributes:
        N:            Lattice size used.
        a:            Lattice spacing.
        eigenvalue:   Lowest eigenvalue of H² (= E_0² for the zero-mode).
        energy:       E_0 = sqrt(max(eigenvalue, 0)) — fermion ground energy.
        eigenvector:  ψ_0 as full (4·N³,) complex vector.
        density:      |ψ_0|² summed over spinor components, shape (N,N,N).
        kink_profile: m(x) array, shape (N,N,N).
        params:       WilsonParams used.
        wall_time_s:  Computation time in seconds.
    """
    N: int
    a: float
    eigenvalue: float
    energy: float
    eigenvector: np.ndarray
    density: np.ndarray
    kink_profile: np.ndarray
    params: WilsonParams
    wall_time_s: float


def find_zero_mode(
    p: WilsonParams,
    *,
    verbose: bool = True,
) -> ZeroModeResult:
    """Find the Jackiw-Rebbi zero-mode via sparse Hermitian eigenvalue solve.

    Algorithm:
        1. Build D_W (Wilson-Dirac sparse matrix, shape 4N³ × 4N³).
        2. Form H² = D_W† D_W (Hermitian, positive-semidefinite).
        3. Call eigsh(H², k=4, which='SM') to find the 4 smallest eigenvalues.
           We request 4 to check there is a single isolated near-zero mode.
        4. Extract eigenvector ψ_0 for the smallest eigenvalue.
        5. Reshape to (4, N, N, N) and compute density |ψ_0|² = Σ_α |ψ_α|².

    Args:
        p:       Wilson lattice parameters.
        verbose: Print progress messages.

    Returns:
        ZeroModeResult containing eigenvalue, eigenvector, and density profile.
    """
    t0 = time.perf_counter()

    if verbose:
        print(f"  Building Wilson-Dirac matrix  N={p.N}  a={p.a:.4f}  ndof={p.ndof:,}")

    D_W = build_wilson_dirac(p)

    if verbose:
        print(f"  D_W built: {D_W.nnz:,} non-zeros  ({time.perf_counter()-t0:.1f}s)")

    t1 = time.perf_counter()
    H2 = build_H_squared(D_W)

    if verbose:
        print(f"  H² built: {H2.nnz:,} non-zeros  ({time.perf_counter()-t1:.1f}s)")

    # Request k=4 smallest eigenvalues.  For large N we only strictly need
    # k=1 but 4 lets us verify the spectral gap is clean.
    k_eigs = min(6, p.ndof - 2)

    t2 = time.perf_counter()

    # MEMORY-SAFE EIGENVALUE STRATEGY (LOBPCG):
    #
    #   shift-invert (sigma=0) is the textbook method for smallest
    #   eigenvalues of a sparse PSD matrix, but it requires an LU
    #   decomposition of (H² - σI) which suffers from catastrophic
    #   fill-in on a 3D lattice (the LU factor is dense banded with
    #   bandwidth ~ N², so the factor takes O(N⁴) memory — for N=20
    #   that's 25 GB even though H² itself is only ~30 MB sparse).
    #
    #   The Lanczos-on-(μI−H²) trick (deflated power iteration) was
    #   tried first but missed the true minima by ~100× because the
    #   spectrum is densely clustered near zero (huge degeneracy from
    #   doubler/zero-mode interaction) and ARPACK could not separate
    #   them with a small Krylov subspace.
    #
    #   LOBPCG (Locally Optimal Block Preconditioned Conjugate Gradient)
    #   is the right tool: it works directly on H² (matrix-vector only),
    #   is provably convergent for PSD operators, and naturally handles
    #   degenerate clusters via the block size k.  Memory stays at
    #   O(sparse-NNZ + k·N³).
    if verbose:
        print(f"  Running LOBPCG (k={k_eigs}, smallest of H²) ...")

    # Initial guess: random orthonormal block.  Match dtype of H² so
    # LOBPCG can stay in complex arithmetic (necessary because the γ
    # matrices give H² complex entries even though it is Hermitian).
    rng = np.random.default_rng(42)
    is_complex = np.iscomplexobj(H2.data)
    if is_complex:
        X0 = (rng.standard_normal((p.ndof, k_eigs))
              + 1j * rng.standard_normal((p.ndof, k_eigs))).astype(np.complex128)
    else:
        X0 = rng.standard_normal((p.ndof, k_eigs)).astype(np.float64)
    # Cheap orthonormalisation
    X0, _ = np.linalg.qr(X0)

    try:
        vals_raw, vecs_raw = spla.lobpcg(
            H2,
            X0,
            tol=1e-9,
            maxiter=2_000,
            largest=False,
            verbosityLevel=0,
        )
    except Exception as exc:
        if verbose:
            print(f"  LOBPCG raised {type(exc).__name__}: {exc}; falling back to deflated ARPACK")
        # Fall back to ARPACK matrix-free deflation if LOBPCG fails.
        rng2 = np.random.default_rng(123)
        v = rng2.standard_normal(p.ndof).astype(np.float64)
        v /= np.linalg.norm(v)
        for _ in range(50):
            v = H2 @ v
            nv = np.linalg.norm(v)
            if nv == 0:
                break
            v /= nv
        mu = 1.10 * float(np.real(v @ (H2 @ v))) + 1.0
        n = p.ndof

        def _mv(x: np.ndarray) -> np.ndarray:
            return mu * x - H2 @ x
        A_shift = spla.LinearOperator(
            (n, n), matvec=_mv, dtype=H2.dtype,
        )
        vals_A, vecs_raw = spla.eigsh(
            A_shift, k=k_eigs, which='LM', tol=1e-9, maxiter=20_000,
        )
        vals_raw = mu - vals_A

    # LOBPCG returns real eigenvalues for any Hermitian operator;
    # promote to complex for the rest of the pipeline.
    vals = np.asarray(vals_raw, dtype=np.complex128)
    vecs = np.asarray(vecs_raw, dtype=np.complex128)

    # Sort by ascending eigenvalue (should already be sorted, but be safe)
    order = np.argsort(np.abs(vals))
    vals = vals[order]
    vecs = vecs[:, order]

    eig0 = float(np.real(vals[0]))
    E0 = float(np.sqrt(max(eig0, 0.0)))

    # ----------------------------------------------------------------------
    # PBC creates a kink-antikink pair (m(x) = M tanh(x/ξ) is not periodic,
    # so going around the box requires one anti-kink to flip the sign back).
    # The Jackiw-Rebbi zero-mode therefore has 4-fold degeneracy:
    #     2 (kink + antikink) × 2 (Dirac spin doublet)
    # ARPACK returns these 4 in an arbitrary basis, so picking vecs[:,0]
    # gives a *random* superposition whose density depends on the seed.
    #
    # The basis-INDEPENDENT density is the trace over the degenerate
    # subspace:   ρ̄(x) = Σ_α |ψ_α(x)|² / d
    # where the sum runs over all eigenvectors with eigenvalue within
    # tol of eig0 (the ones returned by ARPACK with k≥4).
    # This recovers the symmetric kink-antikink density profile.
    # ----------------------------------------------------------------------
    deg_tol = max(10.0 * abs(eig0), 1e-9)
    n_deg = int(np.sum(np.abs(np.real(vals) - eig0) < deg_tol))
    n_deg = max(n_deg, 1)

    psi_4d_subspace = np.stack(
        [vecs[:, i].reshape(4, p.N, p.N, p.N) for i in range(n_deg)],
        axis=0,
    )   # shape (n_deg, 4, N, N, N)

    # ----------------------------------------------------------------------
    # JR-chirality projector.
    #
    # The Jackiw-Rebbi zero-mode in the kink direction (axis k) satisfies
    #     iγ^k · ψ = +ψ      (spinors aligned with the kink slope)
    # Project the degenerate subspace onto this (1+iγ^k)/2 eigenspace
    # before computing the density.
    #
    # In the Dirac-Pauli representation γ^k for k=1,2,3 has γ_k² = -I_4,
    # so iγ^k has γ_k² = +I_4 with eigenvalues ±1.  The +1 eigenspace is
    # the JR bound state.
    #
    # HONEST CAVEAT: this projector helps when LOBPCG actually converges
    # on the JR zero-mode, but Wilson fermions with Dirichlet BCs admit
    # surface-localised states (a sort of "topological-insulator surface
    # mode") whose eigenvalues compete with the bulk kink zero-mode.
    # When LOBPCG settles on a surface mode rather than the bulk kink,
    # the JR projector cannot rescue the centred profile because the
    # surface mode is in a *different* chirality sector.  Resolving this
    # cleanly requires either Domain-Wall Fermions or Overlap Fermions
    # (which have exact lattice chiral symmetry); both are substantial
    # rewrites flagged as future work.  For now the projector is wired
    # in so that *when* the JR mode is found, its density is correctly
    # centred, but the surface-state competition is not eliminated.
    # ----------------------------------------------------------------------
    g_kink = GAMMA[1 + p.kink_axis]   # γ^1, γ^2, or γ^3 along kink direction
    P_JR = 0.5 * (np.eye(4, dtype=np.complex128) + 1j * g_kink)
    # Apply spinor projector to each mode in the subspace
    psi_proj = np.einsum("ab,nbijk->naijk", P_JR, psi_4d_subspace)

    # Density from JR-projected subspace (sum over modes and spinors)
    density_proj = np.sum(np.abs(psi_proj) ** 2, axis=(0, 1))

    # Sanity: if the projection annihilates everything (numerical zero),
    # fall back to the unprojected density rather than producing zeros.
    proj_weight = float(np.sum(density_proj) * p.a ** 3)
    if proj_weight > 1e-12:
        density = density_proj
    else:
        # Try the orthogonal projector (1 − iγ^k)/2; one of the two should
        # carry weight under Dirichlet, both under PBC.
        P_JR_minus = 0.5 * (np.eye(4, dtype=np.complex128) - 1j * g_kink)
        psi_proj_m = np.einsum("ab,nbijk->naijk", P_JR_minus, psi_4d_subspace)
        density_m = np.sum(np.abs(psi_proj_m) ** 2, axis=(0, 1))
        if float(np.sum(density_m) * p.a ** 3) > 1e-12:
            density = density_m
        else:
            # Both projectors annihilate — eigenvector lies entirely outside
            # the JR chirality space.  Fall back to the unprojected trace.
            density = np.sum(np.abs(psi_4d_subspace) ** 2, axis=(0, 1))

    if verbose:
        print(f"  eigsh done  ({time.perf_counter()-t2:.1f}s)")
        print(f"  4 smallest eigenvalues of H²: {np.real(vals[:min(4,len(vals))])}")
        print(f"  Zero-mode eigenvalue λ_0 = {eig0:.6e},  E_0 = {E0:.6e}")
        print(f"  Degenerate subspace: averaging over {n_deg} mode(s)")

    # Keep psi0 for backward compatibility (used by callers expecting one vec)
    psi0 = vecs[:, 0]
    psi_4d = psi0.reshape(4, p.N, p.N, p.N)

    # Normalise density so integral = 1 over the box volume
    norm = float(np.sum(density) * p.a ** 3)
    if norm > 1e-30:
        density /= norm
        psi_4d /= np.sqrt(max(norm / n_deg, 1e-30))

    kink_profile = kink_mass_field(p)
    wall_time = time.perf_counter() - t0

    if verbose:
        print(f"  Total wall time: {wall_time:.1f}s")

    return ZeroModeResult(
        N=p.N,
        a=p.a,
        eigenvalue=eig0,
        energy=E0,
        eigenvector=psi_4d.ravel(),
        density=density,
        kink_profile=kink_profile,
        params=p,
        wall_time_s=wall_time,
    )


# ---------------------------------------------------------------------------
# Jackiw-Rebbi analytic zero-mode profile
# ---------------------------------------------------------------------------

def jackiw_rebbi_analytic(
    p: WilsonParams,
    normalise: bool = True,
) -> np.ndarray:
    """Compute the analytic Jackiw-Rebbi zero-mode density profile.

    For the 1+1D Jackiw-Rebbi model with m(x) = M tanh(x/ξ), the exact
    zero-mode wavefunction is:

        ψ_0(x) ∝ exp(−M·ξ · ln cosh(x/ξ))
                = cosh(x/ξ)^(−Mξ)     [in units where coupling g=1]

    In 3D the kink is planar (infinite transverse extent), so the density
    depends only on the kink-axis coordinate:

        ρ(x_k) ∝ cosh(x_k / ξ)^(−2Mξ)

    Lattice corrections modify this at O(a/ξ) but the overall shape and
    the localisation length are predicted correctly.

    Args:
        p:          Wilson lattice parameters.
        normalise:  Divide by ∫ρ dV so ∫ρ dV = 1.

    Returns:
        Density array of shape (N,N,N), peaked at kink position.
    """
    a = p.a
    coords = (np.arange(p.N) - p.N / 2.0) * a   # centred coordinates

    # 1D profile: cosh^(−2Mξ) at the centre kink.
    exponent = -2.0 * p.M * p.xi
    profile_centre = np.cosh(coords / p.xi) ** exponent   # shape (N,)

    if p.boundary == "dirichlet":
        # Boundary hops are zeroed -> NO antikink at the box edge.
        # The lattice carries a single kink at x=0; analytic template
        # is a single cosh⁻²(x/ξ) peak.
        profile_1d = profile_centre
    else:
        # PBC and APBC: the mass field m(x) = M tanh(x/ξ) keeps its
        # +M → −M jump at the boundary, so the lattice still hosts an
        # antikink there regardless of the spinor BC.  Two-peak template.
        L = p.N * a
        dist_from_boundary = (L / 2.0) - np.abs(coords)
        profile_boundary = np.cosh(dist_from_boundary / p.xi) ** exponent
        profile_1d = profile_centre + profile_boundary

    # Broadcast to 3D transverse to kink axis
    shape = [1, 1, 1]
    shape[p.kink_axis] = p.N
    profile_3d = np.broadcast_to(profile_1d.reshape(shape), (p.N, p.N, p.N)).copy()

    if normalise:
        norm = float(np.sum(profile_3d) * a ** 3)
        if norm > 1e-30:
            profile_3d /= norm

    return profile_3d


# ---------------------------------------------------------------------------
# Coulomb interaction estimate
# ---------------------------------------------------------------------------

def coulomb_interaction_energy(
    density: np.ndarray,
    p: WilsonParams,
    e_charge: float = 1.0,
) -> float:
    """Estimate the Coulomb self-energy of the zero-mode charge distribution.

    Uses the Fourier-space Poisson solver:
        E_Coulomb = (e²/2) ∫∫ ρ(x) G(x−y) ρ(y) dx dy
                  = (e²/2) Σ_k |ρ̃(k)|² / k²    (in 3D, k≠0)

    where G(x) = 1/(4π|x|) is the 3D Coulomb Green's function.

    The zero-mode density ρ is delocalised in the transverse directions
    (it doesn't decay transversely in the planar-kink case), so we work
    with the 1D projection along the kink axis to get a finite estimate:

        ρ_1D(x_k) = ∫ ρ(x_k, x_⊥) d²x_⊥  ×  a²  (transverse integral)

    and use the 1D Coulomb Green's function G_1D(x) = −|x|/2 for the
    charge-line interaction.

    Args:
        density: Normalised density |ψ_0|², shape (N,N,N).
        p:       Wilson lattice parameters.
        e_charge: Fermion electric charge (in lattice units).

    Returns:
        Coulomb interaction energy estimate (dimensionless lattice units).
    """
    a = p.a
    N = p.N

    # Project density onto kink axis
    axes = [0, 1, 2]
    transverse_axes = tuple(ax for ax in axes if ax != p.kink_axis)
    rho_1d = np.sum(density, axis=transverse_axes) * a ** 2   # shape (N,)

    # Fourier-space solve: V(k) = ρ̃(k) / k²  (1D Poisson: d²V/dx² = −ρ)
    rho_k = np.fft.rfft(rho_1d)
    k_vals = np.fft.rfftfreq(N, d=a) * 2.0 * np.pi   # physical wavenumbers

    # Avoid k=0 singularity (charge-neutral system has vanishing k=0 mode)
    k_safe = np.where(np.abs(k_vals) < 1e-12, 1.0, k_vals)
    V_k = np.where(np.abs(k_vals) < 1e-12, 0.0, rho_k / k_safe ** 2)

    E_coulomb = float(
        0.5 * e_charge ** 2
        * np.real(np.sum(np.conj(rho_k) * V_k))
        * a   # 1D integral weight
    )
    return E_coulomb


def alpha_effective(
    E_coulomb: float,
    E_fermion: float,
    pi: float = np.pi,
) -> float:
    """Estimate effective fine structure constant from Coulomb vs fermion energy.

    In natural units the Coulomb energy of a bound fermion scales as:
        E_Coulomb ~ e² / r_0 ~ α · (ħc / r_0)

    Comparing to the fermion binding energy E_fermion = M (mass gap), we get:
        α_eff ≈ E_Coulomb / E_fermion   (very rough dimensional estimate)

    This is a purely lattice-derived estimate, not a renormalised coupling.

    Args:
        E_coulomb: Coulomb interaction energy (same units as E_fermion).
        E_fermion: Fermion mass / binding energy scale.
        pi:        π (default numpy.pi).

    Returns:
        Dimensionless α_eff estimate.
    """
    if abs(E_fermion) < 1e-30:
        return 0.0
    return abs(E_coulomb) / abs(E_fermion)


# ---------------------------------------------------------------------------
# Scaling analysis helper
# ---------------------------------------------------------------------------

def fit_power_law(
    ns: list[int],
    vals: list[float],
) -> tuple[float, float]:
    """Fit eigenvalue ~ C / N^p via log-log linear regression.

    Args:
        ns:   List of lattice sizes.
        vals: Corresponding eigenvalues.

    Returns:
        (C, p) such that val ≈ C · N^(-p).
    """
    log_n = np.log(np.array(ns, dtype=float))
    log_v = np.log(np.array(vals, dtype=float))
    # Linear fit: log(val) = log(C) − p · log(N)
    coeffs = np.polyfit(log_n, log_v, 1)
    p_fit = -coeffs[0]    # exponent (positive for decay)
    C_fit = np.exp(coeffs[1])
    return C_fit, p_fit
