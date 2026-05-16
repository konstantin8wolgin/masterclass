# Overview

*This page is the evolving synthesis of how topics in this wiki connect. The LLM updates it after significant ingests or when a new cross-cutting theme emerges.*

## Domain

This wiki covers **lattice field theory** and **physics-informed kernels** — two areas that intersect around discretization of continuous field theories and machine-learning methods constrained by physical laws.

## Current understanding

### The sign problem and its severity

All papers centre on the same root obstruction: the Boltzmann weight $e^{-S[\phi]}$ (or $\det M(U;\mu)e^{-S_W}$ in QCD) becomes complex or oscillatory, making direct importance sampling impossible. The standard diagnostic is the average phase $\langle e^{i\varphi}\rangle_{pq}$, whose signal-to-noise ratio decays exponentially as $\text{SNR} \sim e^{-(F-F_{pq})/T}$ with volume [Aarts & Sexty]. The silver blaze phenomenon [Balassa] is a sharp manifestation: below the onset chemical potential $\mu_\text{onset}$, all observables are $\mu$-independent despite $\mu\neq 0$, so any inaccuracy in the phase cancellation is immediately visible.

### Shared mathematical structure: flowing distributions

Three distinct approaches in this vault share the same skeleton — they all define a one-parameter family of distributions $p_t$ or integration contours $\mathcal{M}_t$ that interpolates between a tractable reference ($t=0$) and the target ($t=1$):

| Method | Flow equation | Reference → Target |
|--------|--------------|-------------------|
| PIK / PIRG [Ihssen 2025, 2026] | $\frac{dp_t}{dt} = -\frac{\partial}{\partial\phi_i}[\dot\phi_{t,i}\,p_t]$ (Wegner) | real Gaussian → interacting / complex |
| Holomorphic gradient flow [Aarts] | $\frac{dA_\mu}{dt} = -(\partial S/\partial A_\mu)^*$ | original contour $\mathcal{M}_0$ → deformed $\mathcal{M}_t$ |
| Lefschetz thimbles [Aarts] | $\partial_t\tilde\phi = -\overline{(\partial S/\partial\tilde\phi)}$ (same ODE, $t\to\infty$) | original integration cycle → stable thimbles $\mathcal{J}_k$ |

The PIK approach is unique in that the flow equation is inverted: given an interpolating action $S_t$, the kernel $\dot\phi_t$ is solved from the **linear** first-order PDE (Wegner equation)
$$\frac{dS_t(\phi)}{dt} + \dot\phi_t\,\frac{\partial S_t}{\partial\phi} = \frac{\partial\dot\phi_t}{\partial\phi}$$
making each layer an independent linear system $A_t k_t = b_t$ — no backpropagation through time needed.

### Weight-preservation and sign-problem removal

The PIK complex extension [Ihssen 2026] adds a second requirement: the flow must map the complex measure $p(\varphi)d\mu$ to a **real, positive** measure on the deformed manifold $\mathcal{M}_1$ (a PIKfold). The weight-preservation identity
$$\int_{\mathcal{M}_t}d\mu_t\,p_t(\phi)\,\mathcal{O}(\phi) = \int_\mathcal{M}d\mu\,p(\phi)\,\mathcal{O}(\varphi(\phi))$$
guarantees that expectation values computed on the PIKfold agree with those of the original theory. For the harmonic oscillator the PIK is analytic: $\dot\phi_t = -\frac{1}{2}M_t^{-1}\frac{dM_t}{dt}\phi$ where the coupling $\kappa_t = e^{-i\pi t/2}$ rotates the action from Euclidean to Minkowski signature.

### Subtraction as a unifying language

Lawrence & Yamauchi show that **every** infinitesimal contour deformation is a subtraction of the form
$$e^{-S_u} = e^{-S}\bigl(1 + u\cdot\nabla S - \nabla\cdot u\bigr)$$
where the divergence theorem guarantees $\int \nabla\cdot(e^{-S}u)\,d\phi = 0$. The converse fails (not every subtraction is a contour deformation), so subtractions are strictly more general. Neural networks parametrising $u(\phi)$ extend the method to nonperturbative regimes.

### Algebraic bypass: RBF path integrals

Balassa's approach avoids flow entirely. The interaction factor $e^{-S_1}$ is approximated by $K$ Gaussian kernels:
$$e^{-S_{1,x}(m,\lambda)} \approx \sum_{k=1}^K a_k\,e^{-A\sum_a(\phi_{a,x}-c_{a,k})^2}$$
After a unitary transformation $U$ that diagonalises $W = M(\mu) - 2\bar{A}$ in momentum space, the partition function factorises as $Z\approx\prod_p\sum_k Q_k[\tilde\Phi_{a,p},\mu,A,c_{a,k},a_k]$, reducing complexity from $\mathcal{O}(K^{N^4})$ to $\mathcal{O}(KN^4)$.

### Ergodicity as an independent obstruction

Even when the sign problem is handled, MCMC can fail due to ergodicity violations at infinite potential barriers (e.g. zeros of $\det M$ in the Hubbard model with exponential discretisation). Temmen et al. fix this with multiplicative proposals $\phi' = e^\gamma\phi$, $\gamma\sim\mathcal{N}(0,\sigma_R^2)$, which jump over barriers without computing gradients. The acceptance criterion $\alpha_R = \min(1,e^{-\Delta S + d\gamma})$ includes the Jacobian $e^{d\gamma}$ from the volume-changing map, preserving detailed balance.

## Open threads

- **PIK existence criteria**: a PIKfold $\mathcal{M}_1$ exists only when the Wegner equation has a real solution at $t=1$; systematic criteria for general gauge theories are unknown.
- **Fermionic extensions**: all current PIK demonstrations are bosonic. Finite-density QCD requires handling a complex fermion determinant $\det M(U;\mu)$ — an open problem for both PIK and RBF approaches.
- **Subtraction + contour combination**: Lawrence & Yamauchi show subtractions are strictly more general than contour deformations; optimal combinations with holomorphic flows for severe sign problems are unexplored.
- **Scaling of RBF kernels**: the required number $K$ of Gaussian kernels with lattice volume and strong coupling $\lambda$ is not theoretically bounded.
- **Radial updates on compact fields**: the multiplicative proposal $\phi' = e^\gamma\phi$ is specific to non-compact fields; extension to $U(1)$ or $SU(N)$ gauge fields is open.

## Key tensions

- **Complex Langevin vs PIK**: the PIK complex architecture correctly handles the real-time harmonic oscillator where Complex Langevin (CL) fails [Ihssen 2026]; however, it is not known whether PIK replaces CL in all cases where CL fails, or only a subset. CL is simpler to implement for gauge theories.
- **Lefschetz thimbles vs PIKfolds**: both deform the integration manifold into $\mathbb{C}^N$; thimbles satisfy $\text{Im}\,S = \text{const}$ and in principle completely eliminate the intra-thimble sign problem, but suffer from a residual inter-thimble global sign problem and require computing Jacobians. PIKfolds are defined by the Wegner equation condition and do not automatically satisfy $\text{Im}\,S = \text{const}$; their advantage is the constructive analytic flow.
- **NP-hardness vs method optimism**: Troyer-Wiese 2005 [cited in Aarts] proves the sign problem is NP-hard in general; yet all specific methods in this vault achieve correct results on their test models. The resolution is that NP-hardness is a worst-case statement — structured physical theories may lie outside the hard instances.
