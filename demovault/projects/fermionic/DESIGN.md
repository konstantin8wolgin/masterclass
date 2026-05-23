---
title: "Fermionic one-site sampling via Gaussian decomposition + PIK"
date: 2026-05-22
status: design
project: projects/fermionic
---

# Context

The other vault projects (`phi4imaginary`, `linearphi4`) ship as numbered Julia scripts and saved JLD2 checkpoints. For the fermionic project we want a single executable Jupyter notebook (Julia kernel via IJulia) that runs every step of the sampling procedure side-by-side, with explicit verification gates between phases.

The physics target is the **one-site Hubbard model with the spin-1/2 fermionic action**:

$$
S(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + e^{\beta(i\phi + \mu + U/2)}\right]
$$

$S$ is complex on $\phi \in \mathbb{R}$, so $e^{-S}/Z$ is not a probability measure — this is the sign problem at one site. The notebook constructs a sign-problem-free sampler by:

1. Expanding the Boltzmann weight algebraically as **three** complex shifted Gaussians (binomial square — exact, not Taylor).
2. Identifying the physics-informed kernel $\dot\phi = iU$ that connects each component to a single real-Gaussian base $\mathcal{N}(0, U/\beta)$.
3. Drawing $M$ real-Gaussian samples once and evaluating every observable as a weighted sum at the three deterministic imaginary shifts.

No MCMC, no thimble integration, no neural-network training is needed at one site. The KL-mixture learner for the multi-site / unknown-decomposition case is recorded in the agent's memory (`memory/kl_gaussian_mixture_decomposition.md`) for future reuse and is **out of scope** for this notebook.

# Theory recap (concise)

Algebraic identity (this is what every cell rests on):

$$
e^{-S(\phi)} = \sum_{j=0}^{2} w_j\, e^{-\frac{1}{2\sigma^2}(\phi - i m_j)^2},
\qquad
\sigma^2 = \frac{U}{\beta},
$$

with

| $j$ | center $m_j$ | weight $w_j$ |
|---|---|---|
| 0 | $0$ | $1$ |
| 1 | $U$ | $2\,e^{\beta\mu}$ |
| 2 | $2U$ | $e^{\beta(2\mu - U)}$ |

The integrated partition function reproduces the Fock trace exactly: $Z = 1 + 2 e^{\beta\mu} + e^{\beta(2\mu - U)}$.

Observable estimator (Fukushima 2026 eq. 20), with $\phi_0^{(i)} \sim \mathcal{N}(0, \sigma^2)$:

$$
\langle\mathcal{O}\rangle \;=\; \frac{\sum_{j=0}^{2} w_j\,\frac{1}{M}\sum_{i=1}^{M} \mathcal{O}(\phi_0^{(i)} + i m_j)}{\sum_{j=0}^{2} w_j}.
$$

**Critical: `O(φ)` must be holomorphic on the strip $\mathrm{Im}\,\phi \in [0, 2U]$.** The naive per-φ trace-ratio expressions $n(\phi) = 2\sigma(\beta(i\phi+\mu+U/2))$ etc. have poles at $\mathrm{Im}\,\phi = \mu + U/2$, which lies *inside* the contour-shift strip whenever $\mu < 3U/2$ — i.e. 7 of the 8 grid cases. The deformed integral then picks up residues and the estimator returns the wrong answer (verified empirically: $\langle n^2\rangle$ via the sigmoid form fails 8/8 cases by 100σ to 10⁶σ; even $\langle n\rangle$ fails 2/8 cases catastrophically).

The fix is the **holomorphic representation from Fukushima eq. 21**: the auxiliary field *is* the density operator, $\langle n\rangle = -i\langle\phi\rangle/U$. Derivation by direct Gaussian integration of $e^{-S}\cdot\phi^k$ against $(1+\xi)^2$:

| observable | PIK operator $O(\phi)$ |
|---|---|
| $n$ | $-i\phi/U$ |
| $n^2$ | $1/(\beta U) - \phi^2/U^2$ |
| $n_\uparrow n_\downarrow$ | $(O_{n^2}(\phi) - O_n(\phi))/2$ |

All three are polynomials in $\phi$ — analytic everywhere — so every contour shift is exact. **A side-effect:** at one site $\langle n\rangle$ is exactly machine-precision under PIK (every per-sample contribution has the same real part $\mathrm{Im}(\phi_i)/U = 0$ plus the constant Fock value), so its sampling variance is zero. The 1/√M convergence study therefore uses $\langle n^2\rangle$, which inherits genuine variance from $\phi^2$.

# Notebook layout

Single notebook: `demovault/projects/fermionic/fermionic_one_site.ipynb`. Five top-level sections, each ending with an explicit `@assert` cell so the user pauses and confirms before moving on.

| § | Title | Purpose | Verification gate |
|---|---|---|---|
| 0 | Setup & action | Imports, parameters, `S(φ)` and `e^{-S(φ)}` as complex Julia functions. Plot `Re S`, `Im S`, `|e^{-S}|` over a real-φ grid for the headline case. | Visual inspection only. |
| 1 | Analytic 3-Gaussian decomposition | Compute `w_j, m_j, σ²` from the formulas above. Define mixture density `p_mix(φ)`. Plot `p_mix` vs `e^{-S}` (both complex; show Re and Im) for the headline case. | `maximum(abs.(p_mix.(grid) .- exp.(-S.(grid)))) < 1e-12` over `grid = -L:dφ:L`. |
| 2 | Partition function & analytic observables | Closed-form integrals (Gaussian formula) → analytic `Z`, `⟨n⟩`, `⟨n²⟩`, `⟨n↑n↓⟩` for all 8 grid cases. Cross-check against direct Fock-trace formulas. | Per case: `abs(Z_gauss - Z_trace) < 1e-12` and same for `⟨n⟩`. |
| 3 | PIK sampler | For each of the 8 cases: draw `φ₀ ∼ 𝒩(0, U/β)` at `M = 10⁵`; build the three deterministic complex images `φ₀ + ij·U`; compute the estimator using the *holomorphic* PIK operators (see above) for `⟨n⟩`, `⟨n²⟩`, `⟨n↑n↓⟩`. Convergence study at headline uses `⟨n²⟩` at `M ∈ {10³, 10⁴, 10⁵, 10⁶}` — `⟨n⟩` is zero-variance at one site. | `⟨n⟩` exact to `1e-10` (zero-variance gate); `⟨n²⟩` and `⟨n↑n↓⟩` within `5·SE` of analytic for every case; convergence SE scales clean `1/√M`. |
| 4 | Distribution comparison | Reproduce a `04_t1_comparison`-style figure for the headline case: histogram of weighted samples vs analytic mixture density. Also: 2×4 grid of subplots (one per case) showing `⟨n⟩` analytic vs sampled with error bars. | Visual + KS-style numeric agreement on the histogram. |

# Parameter grid (8 cases)

Cartesian product:
- $U \in \{1, 4\}$
- $\beta \in \{2, 10\}$
- $\mu \in \{0, 1\}$

Indexed `case_1 ... case_8`. **Headline case** for the convergence study and the marquee histogram: $(U, \beta, \mu) = (1, 10, 1)$ — same as Ihssen 2026 Fig. 6.

# File layout

```
demovault/projects/fermionic/
├── DESIGN.md                          this file
├── Project.toml                       IJulia, Plots, Distributions, QuadGK, JLD2, Random
├── Manifest.toml
├── fermionic_one_site.ipynb           the notebook
├── data/
│   ├── mixture_params.jld2            (w_j, m_j, σ²) per case
│   ├── samples_headline.jld2          φ₀ samples + shifted images at M=10⁵
│   ├── observables.jld2               table of ⟨n⟩, ⟨n²⟩, ⟨n↑n↓⟩, Z per case
│   ├── 01_action_check.png            Re/Im S, |e^{-S}|
│   ├── 01_mixture_vs_target.png       p_mix vs e^{-S} pointwise
│   ├── 03_convergence.png             1/√M scaling, headline case
│   └── 04_t1_comparison.png           histogram vs theory, headline case
└── 04_t1_comparison.png               existing reference (kept untouched as the visual target)
```

# Implementation notes

- Reuse vault conventions: scripts/notebook activate the local environment via `using Pkg; Pkg.activate("."); Pkg.instantiate()`; data saved as JLD2; PNGs in `data/`.
- Complex arithmetic throughout (mirrors `phi4imaginary`); no neural networks anywhere.
- The PIK is hardcoded as the constant shift $\dot\phi = iU$ — derivable in two lines from the Wegner equation on each component but not re-derived in the notebook.
- Random seed fixed (`Random.seed!(20260522)`) so the headline figure is reproducible.

# Verification end-to-end

After running the notebook in order:

1. Phase 1 assert passes (algebra correct).
2. Phase 2 assert passes for all 8 cases (analytic integration correct).
3. Phase 3 headline convergence figure shows clean `1/√M` slope.
4. Phase 3 8-case table: sampled `⟨n⟩` within `5σ_stat` of analytic for every case.
5. Phase 4 headline histogram visually matches `04_t1_comparison.png` style and analytic density.

If any gate fails the notebook stops via `@assert`; the user re-runs and inspects.

# Out of scope (intentional)

- KL-fit Gaussian-mixture learner for actions without a closed-form decomposition (recorded in agent memory for future use).
- Multi-site Hubbard.
- Lefschetz thimble construction via PIK (covered in `phi4imaginary`).
- PIKfold ripping diagnostics (one-site, μ ≤ 1 stays safely below the pole strip; not needed here).
