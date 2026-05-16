# Operation Log

Append-only record of all wiki operations. Each entry begins with a parseable prefix.

Format:
- `## [YYYY-MM-DD] ingest | <source title> | <N pages updated>`
- `## [YYYY-MM-DD] query  | <question summary> | <pages consulted>`
- `## [YYYY-MM-DD] lint   | <summary of issues found/fixed>`

---

## [2026-05-05] init | Wiki created | 3 pages created

## [2026-05-05] ingest | Lawrence & Yamauchi 2023 — Deep Learning of Fermion Sign Fluctuations | 1 source + sign-problem, average-phase, contour-deformation, subtraction-method, holomorphic-gradient-flow, thirring-model-action
Key math: subtracted Boltzmann factor $e^{-S_u} = e^{-S}(1 + u\cdot\nabla S - \nabla\cdot u)$; any contour deformation is a subtraction; Thirring action $S = \sum \frac{2}{g^2}(1-\cos A_\nu) - \log\det K[A]$.

## [2026-05-05] ingest | Temmen et al. 2025 — Fully ergodic simulations using radial updates | 1 source + ergodicity-violation, hybrid-monte-carlo, radial-updates, radial-update-acceptance, hubbard-action
Key math: proposal $\phi' = e^\gamma\phi$, acceptance $\alpha_R = \min(1, e^{-\Delta S + d\gamma})$; Hubbard action $S[\phi] = \frac{1}{2U\Delta_t}\sum\phi_{tx}^2 - \log\det M[\phi|\kappa]\det M[-\phi|-\kappa]$; $\sigma_R\propto 1/\sqrt{d}$ scaling.

## [2026-05-05] ingest | Balassa 2026 — Addressing the sign-problem with RBF neural networks | 1 source + silver-blaze-problem, rbf-path-integral, partition-function-lattice
Key math: $e^{-S_{1,x}} \approx \sum_k a_k e^{-A\sum_a(\phi_{a,x}-c_{a,k})^2}$; $Z\to\prod_p\sum_k Q_k[...]$ after diagonalising $W = M(\mu)-2\bar{A}$; complexity $\mathcal{O}(KN^4)$; $\mu_c = 1.17\pm 0.018$.

## [2026-05-05] ingest | Ihssen, Kapust, Pawlowski 2025 — Generative sampling with physics-informed kernels | 1 source + physics-informed-kernels, physics-informed-rg-flows, wegner-equation
Key math: Wegner equation $\frac{dS_t}{dt} + \dot\phi_t\frac{\partial S_t}{\partial\phi} = \frac{\partial\dot\phi_t}{\partial\phi}$; continuity equation $\frac{dp_t}{dt} = -\frac{\partial}{\partial\phi_i}[\dot\phi_{t,i}p_t]$; linear system $A_t k_t = b_t$ per layer.

## [2026-05-05] ingest | Ihssen, Kapust, Pawlowski 2026 — Solving sign problems with physics-informed kernels | 1 source + lefschetz-thimbles, complex-langevin-dynamics (updated)
Key math: weight-preservation $\int_{\mathcal{M}_t}d\mu_t\,p_t\mathcal{O} = \int_\mathcal{M}d\mu\,p\,\mathcal{O}\circ\varphi$; harmonic oscillator PIK $\dot\phi_t = -\frac{1}{2}M_t^{-1}\dot M_t\phi$; $\kappa_t = e^{-i\pi t/2}$ rotates Euclidean to Minkowski.

## [2026-05-05] ingest | Aarts & Sexty 2026 — Lattice field theories with a sign problem | 1 source + overlap-problem; sign-problem, contour-deformation, lefschetz-thimbles, complex-langevin updated
Key math: $Z = \int DU\det M(U;\mu)e^{-S_W}$; SNR $\sim e^{-(F-F_{pq})/T}$; thimble flow $\partial_t\tilde\phi = -\overline{\partial_{\tilde\phi}S}$; holomorphic flow $\partial_t A_\mu = -(\partial_A S)^*$; Langevin $dx/dt = -\partial_x S + \eta$.

## [2026-05-05] ingest | Batch ingest of 6 papers | 6 sources, 10 concepts, 5 methods, 5 equations created; index and log updated
