# Wiki Index

A catalog of all pages in this wiki. Updated by the LLM on every ingest.

## Sources

| Page | Summary | Date Added |
|------|---------|------------|
| [[lawrence-2023-deep-learning-fermion-sign]] | $e^{-S_u}=e^{-S}(1+u\cdot\nabla S-\nabla\cdot u)$; every contour deformation is a subtraction; NN trains $u$ by gradient $\partial C/\partial\lambda = -\langle\text{Re}\,S\rangle_Q$; Thirring $S=\sum\frac{2}{g^2}(1-\cos A_\nu)-\log\det K[A]$ | 2026-05-05 |
| [[temmen-2025-fully-ergodic-radial-updates]] | $\phi'=e^\gamma\phi$, $\alpha_R=\min(1,e^{-\Delta S+d\gamma})$; Hubbard $S[\phi]=\frac{1}{2U\Delta_t}\sum\phi_{tx}^2-\log\det M[\phi|\kappa]\det M[-\phi|-\kappa]$; opt. trajectory $T=\frac\pi2\sqrt{U\beta/N_t}$; $\sigma_R^\text{opt}\propto d^{-1/2}$ | 2026-05-05 |
| [[balassa-2026-sign-problem-rbf-networks]] | $e^{-S_{1,x}}\approx\sum_k a_k e^{-A\sum_a(\phi_{a,x}-c_{a,k})^2}$; $Z\to\prod_p\sum_k Q_k$ after diagonalising $W=M(\mu)-2\bar{A}$; $\mathcal{A}_p=-2A+4\sum\sin^2(p_i/2)+2[1-\cosh\mu\cos p_4]$; $\mathcal{O}(KN^4)$; $\mu_c=1.17\pm0.018$ | 2026-05-05 |
| [[ihssen-2025-generative-sampling-piks]] | Wegner eq. $\frac{dS_t}{dt}+\dot\phi_t\frac{\partial S_t}{\partial\phi}=\frac{\partial\dot\phi_t}{\partial\phi}$ gives linear system $A_t k_t=b_t$ per layer; basis $K_{1,t}=\phi$, $K_{j,t}=\sin(\omega_j\phi)$; RK4 $\Delta t=1/25$; correction via Wendland RBF collocation | 2026-05-05 |
| [[ihssen-2026-solving-sign-problems-piks]] | Weight-preservation: $\int_{\mathcal{M}_t}d\mu_t\,p_t\mathcal{O}=\int d\mu\,p\,\mathcal{O}\circ\varphi$; PIKfold ODE $d\phi_I/d\phi_R=-\text{Im}[p_t]/\text{Re}[p_t]$; harmonic osc. PIK $\dot\phi_t=-\frac12M_t^{-1}\dot M_t\phi$; $\kappa_t=e^{-i\pi t/2}$, $N_t\geq512$ steps | 2026-05-05 |
| [[aarts-2026-lattice-sign-problem-review]] | $Z=\int DU\det M(U;\mu)e^{-S_W}$; SNR$\sim e^{-(F-F_{pq})/T}$; thimble flow $\partial_t\tilde\phi=-\overline{\partial_{\tilde\phi}S}$; holomorphic flow $\partial_t A_\mu=-(\partial_A S)^*$; CL $dx/dt=-\partial_x S+\eta$; review | 2026-05-05 |

## Concepts

| Page | Summary | Sources | Date Added |
|------|---------|---------|------------|
| [[sign-problem]] | Complex Boltzmann weight makes importance sampling impossible; SNR decays exponentially with volume | 4 | 2026-05-05 |
| [[physics-informed-kernels]] | PIK: layer-independent generative flow solving Wegner equation; resolves OOD and sign problems | 2 | 2026-05-05 |
| [[physics-informed-rg-flows]] | PIRG flows in distribution space; Wegner equation as governing PDE; basis for PIK architecture | 2 | 2026-05-05 |
| [[ergodicity-violation]] | MCMC trapped by infinite potential barriers; biased observables; fixed by radial updates | 2 | 2026-05-05 |
| [[silver-blaze-problem]] | Observables independent of μ for μ < μc at T=0; requires exact sign cancellations | 2 | 2026-05-05 |
| [[average-phase]] | Standard diagnostic for sign problem severity; target of subtraction and contour deformation methods | 2 | 2026-05-05 |
| [[overlap-problem]] | Poor overlap between phase-quenched and full theory; exponentially bad reweighting | 1 | 2026-05-05 |
| [[contour-deformation]] | Deform integration contour into C^N; reduce phase fluctuations; every deformation is a subtraction | 2 | 2026-05-05 |
| [[lefschetz-thimbles]] | Contours where Im S = const; reduce intra-thimble sign problem; two residual sign problems remain | 2 | 2026-05-05 |
| [[complex-langevin-dynamics]] | Stochastic process in complexified field space; convergence not guaranteed; fails for some real-time cases | 2 | 2026-05-05 |

## Methods

| Page | Summary | Sources | Date Added |
|------|---------|---------|------------|
| [[hybrid-monte-carlo]] | MD + Metropolis algorithm for sampling e^{-S}; susceptible to ergodicity violations | 2 | 2026-05-05 |
| [[radial-updates]] | Multiplicative MH proposal φ' = e^γφ; overcomes infinite potential barriers; detailed balance satisfied | 1 | 2026-05-05 |
| [[subtraction-method]] | Subtract g = ∇·v from Boltzmann factor; physical results unchanged; neural network parametrisation of v | 1 | 2026-05-05 |
| [[rbf-path-integral]] | Expand e^{-S1} as sum of Gaussians; momentum-space diagonalisation; O(KN^4) complexity | 1 | 2026-05-05 |
| [[holomorphic-gradient-flow]] | Flow into C^N along −(∂S/∂A)*; deformed contours; approaches Lefschetz thimbles at t→∞ | 2 | 2026-05-05 |

## Equations

| Page | Summary | Sources | Date Added |
|------|---------|---------|------------|
| [[partition-function-lattice]] | Z = ∫Dφ e^{-S}; expectation values; phase-quenched theory; SNR decay | 3 | 2026-05-05 |
| [[wegner-equation]] | dS/dt + φ̇·∂S/∂φ = ∂φ̇/∂φ; linear first-order PDE for PIK given S_t | 2 | 2026-05-05 |
| [[radial-update-acceptance]] | αR = min(1, e^{-ΔS + dγ}); Jacobian from multiplicative map; detailed balance | 1 | 2026-05-05 |
| [[thirring-model-action]] | S = Σ 2/g²(1-cos Aν) − log det K[A]; staggered Thirring in 1+1D; μ-dependent sign problem | 1 | 2026-05-05 |
| [[hubbard-action]] | S[φ] = φ²/2UΔt − log det M[φ|κ]det M[−φ|−κ]; exponential discretisation; infinite barriers | 1 | 2026-05-05 |



## Creative

| Page | Summary | Date Added |
|------|---------|------------|
| *(none yet)* | | |

## Problems

| Page | Status | Date Added |
|------|--------|------------|
| *(none yet)* | | |
