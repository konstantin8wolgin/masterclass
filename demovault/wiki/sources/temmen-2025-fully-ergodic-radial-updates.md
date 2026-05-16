---
title: "Fully ergodic simulations using radial updates"
authors: [Finn L. Temmen, Evan Berkowitz, Anthony Kennedy, Thomas Luu, Johann Ostmeyer, Xinhao Yu]
year: 2025
type: paper
tags: [hmc, ergodicity, radial-updates, hubbard-model, metropolis-hastings]
source_file: raw/papers/2504.01575v2.pdf
---

## Summary

Proposes **radial updates** — a multiplicative Metropolis-Hastings step interleaved with HMC — to resolve ergodicity violations from infinite potential barriers in non-compact field spaces. The proposal $\phi' = e^\gamma\phi$ rescales the entire field configuration and jumps over barriers without computing gradients. Demonstrated on the Hubbard model and the C$_{20}$H$_{12}$ perylene molecule.

**When to use**: whenever HMC is trapped by manifolds of vanishing fermion determinant, topological freezing, or any codimension-1 surface where $S \to \infty$. Radial updates fix ergodicity at comparably low cost.

---

## 1. HMC review

Compute expectation values:
$$\langle\mathcal{O}\rangle = Z^{-1}\int\mathcal{D}\phi\,\mathcal{O}[\phi]\,e^{-S[\phi]}, \quad Z = \int\mathcal{D}\phi\,e^{-S[\phi]}. \tag{1}$$

HMC introduces conjugate momenta $\pi_i$ and the artificial Hamiltonian:
$$H[\phi,\pi] = \frac{1}{2}\sum_{i=1}^d\pi_i^2 + S[\phi]. \tag{5}$$

Evolve $(\phi,\pi)$ by Hamilton's equations for fictitious time $\tau$:
$$\frac{d\pi_i}{d\tau} = -\frac{\partial H}{\partial\phi_i} = -\frac{\partial S}{\partial\phi_i}, \qquad \frac{d\phi_i}{d\tau} = \frac{\partial H}{\partial\pi_i} = \pi_i. \tag{6}$$

The MD trajectory is integrated numerically (leapfrog) for length $T$, divided into $N_\text{MD}$ steps of size $\epsilon = T/N_\text{MD}$. Accept/reject with:
$$\alpha_\text{HMC} = \min\!\left(1, e^{-\Delta H}\right), \quad \Delta H = H[\tilde\phi,\tilde\pi] - H[\phi,\pi]. \tag{7}$$

**Ergodicity problem**: if $S$ diverges on a codimension-1 manifold, the force $F_i = -\partial S/\partial\phi_i$ diverges there and the MD trajectory can never cross it, partitioning configuration space.

---

## 2. Radial update — formulation

A **radial update** is a multiplicative MH proposal on non-compact field $\phi = (\phi_1,\ldots,\phi_d)$.

**Algorithm** (Section III.A):

1. Draw $\gamma \sim \mathcal{N}(0, \sigma_R^2)$ (zero-mean Gaussian with standard deviation $\sigma_R$)
2. Propose $\phi' = e^\gamma\phi$, i.e. $\phi'_i = e^\gamma\phi_i$ for all $i$. This rescales the radius:
   $$R' = e^\gamma R, \qquad R = \sqrt{\sum_i\phi_i^2}. \tag{8/9}$$
3. Accept with probability:
   $$\alpha_R = \min\!\left(1,\, e^{-\Delta S + d\gamma}\right), \qquad \Delta S = S[\phi'] - S[\phi]. \tag{10}$$
   The factor $e^{d\gamma}$ is the **Jacobian** of the multiplicative map ($d$-dimensional volume element transforms as $d\phi' = e^{d\gamma}d\phi$).

The combined HMC + radial update algorithm satisfies **detailed balance**:
$$p[\phi]\,\Omega(\phi\to\phi') = p[\phi']\,\Omega(\phi'\to\phi). \tag{3}$$

**Why it works**: the proposal $\phi' = e^\gamma\phi$ does not require computing $\partial S/\partial\phi$, so it is agnostic to the diverging force at barriers and can freely cross them.

---

## 3. Tuning $\sigma_R$

The integrated autocorrelation time as a function of $\sigma_R$ follows the ansatz:
$$\tau_\text{int}(\sigma_R) = a\sigma_R^{-2} + b + c\sigma_R. \tag{25}$$

- Small $\sigma_R$: diffusive regime, $\tau_\text{int} \propto \sigma_R^{-2}$
- Large $\sigma_R$: acceptance rate drops, $\tau_\text{int} \propto \sigma_R$
- Optimum at $\sigma_R^\text{(min)} = (2a/c)^{1/3}$

**Scaling with dimensionality**:
$$\sigma_R^\text{(min)}(d) \propto d^{-1/2}. \tag{26}$$

In practice, aim for a radial acceptance rate of ~30% to be near the optimal $\sigma_R$.

**Practical implementation**: interleave one radial update before each HMC step (or after — order does not matter for detailed balance). Both steps are individually MH-correct.

---

## 4. Toy model benchmark

Distribution $p[x] \propto \prod_{i=1}^d\cos^2(x_i)\,e^{-\beta x_i^2}$ with infinite barriers at $x_i = \frac{2k+1}{2}\pi$. For $d=2$, $\beta=0.125$, $N_\text{MD}=12$, $T=1$: radial update with $\sigma_R = 1.75$ resolves ergodicity violations completely.

---

## 5. Hubbard model

**Hamiltonian** (particle/hole basis):
$$H = H_K + H_U = -\kappa\sum_{\langle x,y\rangle}(a_x^\dagger a_y - b_x^\dagger b_y) + \frac{U}{2}\sum_x(a_x^\dagger a_x - b_x^\dagger b_x)^2. \tag{12}$$

Parameters: hopping $\kappa$, on-site interaction $U$.

**Path integral via Hubbard-Stratonovich**: apply second-order Suzuki-Trotter with $N_t$ time slices ($\Delta_t = \beta/N_t$), then continuous HS transformation:
$$\exp\!\left\{-\frac{\Delta_t U}{2}\sum_x(a_x^\dagger a_x - b_x^\dagger b_x)^2\right\} \propto \int_{-\infty}^\infty\prod_x d\phi_{tx}\,\exp\!\left\{-\frac{1}{2U\Delta_t}\phi_{tx}^2 - i\phi_{tx}(a_x^\dagger a_x - b_x^\dagger b_x)\right\}. \tag{14}$$

**Effective (exponential discretisation) action**:
$$S[\phi] = \frac{1}{2U\Delta_t}\sum_{t,x}\phi_{tx}^2 - \log\!\left(\det M[\phi|\kappa]\,\det M[-\phi|-\kappa]\right). \tag{15}$$

**Fermion matrix** (exponential discretisation):
$$M[\kappa]_{tx,t'y} = \delta_{tt'}\delta_{xy} - (e^{\kappa h})_{xy}\,e^{i\phi_{tx}}\,\mathcal{B}_{t'}\,\delta_{t',t+1}, \tag{16}$$
where $h$ is the hopping matrix, $\mathcal{B}_t = +1$ for $0 < t < N_t$ and $\mathcal{B}_0 = -1$ (anti-periodic boundary conditions).

**Ergodicity problem in Hubbard model**: $\det M[\phi|\kappa]\det M[-\phi|-\kappa] = 0$ on a manifold of codimension 1, creating infinite potential barriers that partition the $(\phi_1,\phi_2)$-plane into diagonal bands.

**Optimal trajectory length** (eliminates harmonic autocorrelations):
$$T = \frac{\pi}{2}\sqrt{U\beta/N_t}. \tag{17}$$

**Benchmark parameters** (two-site model): $U=18$, $\kappa=1$, $\beta=1$, $N_t=1$ (then scaled up). $N_\text{MD}=60$ for $>99\%$ acceptance. Radial $\sigma_R = 1.8$ for $N_t=1$; $\sigma_R = 0.6$ for $N_t=8$. Record $10^4$ configurations saving every 10th HMC step.

**Observable**: spatial correlator $C_{xy}(\tau) = \langle a_x(\tau)a_y^\dagger(0)\rangle = \langle\sum_t(M^{-1})_{tx,(t+\tau)y}\rangle$, projected to irrep $C_-(\tau)$ via parity vector $\mathcal{A}^- = \frac{1}{\sqrt{2}}(1,-1)^\top$.

---

## 6. Algorithm (Julia implementation outline)

```julia
# Radial update
function radial_update!(ϕ, S, σ_R)
    d = length(ϕ)
    γ = σ_R * randn()          # draw proposal
    ϕ_new = exp(γ) .* ϕ        # multiplicative rescaling
    ΔS = S(ϕ_new) - S(ϕ)
    α = min(1.0, exp(-ΔS + d*γ))  # Jacobian factor: e^{dγ}
    if rand() < α
        ϕ .= ϕ_new             # accept
    end
end

# Combined HMC + radial update
function combined_step!(ϕ, π, S, ∇S, T, N_MD, σ_R)
    radial_update!(ϕ, S, σ_R)   # radial first
    hmc_step!(ϕ, π, S, ∇S, T, N_MD)  # then HMC leapfrog
end
```

**Key Julia notes**:
- `det M[ϕ|κ]`: use `LinearAlgebra.det` or `logdet` for numerical stability
- For large $N_t$: $M$ is block-tridiagonal → use sparse representation
- Tune $\sigma_R$ by scanning $\tau_\text{int}$ and fitting ansatz (25)
- $\sigma_R^\text{opt} \approx \alpha d^{-0.5}$ for initial guess

---

## Notation summary

| Symbol | Meaning |
|--------|---------|
| $R = \sqrt{\sum_i\phi_i^2}$ | radial coordinate in field space |
| $\gamma \sim \mathcal{N}(0, \sigma_R^2)$ | radial proposal variable |
| $\phi' = e^\gamma\phi$ | proposed configuration |
| $\alpha_R = \min(1, e^{-\Delta S + d\gamma})$ | radial acceptance probability |
| $\sigma_R$ | proposal standard deviation (key tuning parameter) |
| $M[\phi|\kappa]$ | fermion matrix (exponential discretisation) |
| $\tau_\text{int}$ | integrated autocorrelation time |
| $\tau_\text{int}(\sigma_R) = a\sigma_R^{-2} + b + c\sigma_R$ | autocorrelation ansatz |

---

## Open questions

- Optimal ratio of HMC steps to radial updates as function of system size
- Extension to compact fields ($U(1)$, $SU(N)$) — multiplicative map is specific to non-compact fields
- Combination with normalising flows or parallel tempering
- Theoretical bound on autocorrelation improvement from radial updates

---

## Links

- [[hybrid-monte-carlo]]
- [[radial-updates]]
- [[ergodicity-violation]]
- [[radial-update-acceptance]]
- [[hubbard-action]]
