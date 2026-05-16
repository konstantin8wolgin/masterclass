---
title: "Physics-informed kernels (PIKs)"
tags: [physics-informed-kernels, generative-network, renormalisation-group, sign-problem]
sources: [
  "[[ihssen-2025-generative-sampling-piks]]",
  "[[ihssen-2026-solving-sign-problems-piks]]"
]
---

## Definition

A physics-informed kernel (PIK) is a map $\dot{\phi}_t: I \times \mathcal{C} \to \mathcal{C}$ that characterises the infinitesimal transformation of a field $\phi$ at RG-time $t \in [0,1]$, satisfying the Wegner equation for a given action $S_t(\phi)$. The PIK architecture uses a sequence of such kernels to build a generative flow from a simple base distribution $p_0(\varphi)$ to a target distribution $p_1(\phi) = p(\phi) \propto e^{-S(\phi)}$. The key property is that the kernel at each time $t$ is determined by a simple linear first-order PDE that is independent of all other times.

## Intuition

Standard generative networks (normalising flows) learn a global map between two distributions by training a deep neural network end-to-end. The layers are entangled: changing the kernel at time $t_1$ affects the solution at $t_2 > t_1$. PIKs break this entanglement by using the physics — specifically, the action $S_t(\phi)$ is analytically given, and the Wegner equation converts the global map-finding problem into a sequence of independent local PDE solutions. This means: (1) no backpropagation through time, (2) OOD corrections are analytically accessible, (3) the distribution is correct by construction rather than by fitting.

## Key results

- The Wegner equation for the PIK: $\frac{dS_t}{dt} + \dot{\phi}_t \frac{\partial S_t}{\partial\phi} = \frac{\partial \dot{\phi}_t}{\partial\phi}$ — a linear first-order PDE for $\dot{\phi}_t$ given $S_t$.
- Layer independence: kernel parameters $\theta_t$ at different times are not entangled (no coupling through the loss function).
- Kernel expansion: $\dot{\phi}_t(\phi) = \sum_j k_{j,t} K_{j,t}(\phi)$; coefficients from independent linear systems $A_t k_t = b_t$.
- The complex PIK-architecture extends this to complex $S_t$, yielding a weight-preserving deformation that maps complex distributions to sign-problem-free manifolds ([[ihssen-2026-solving-sign-problems-piks]]).
- For the harmonic oscillator, the PIK is purely analytic: $\dot{\phi}_t = -\frac{1}{2}M_t^{-1}\frac{dM_t}{dt}\phi$.
- Weight-preserving property: $\int_{\mathcal{M}_1} d\mu_1\, p_1(\phi)\,\mathcal{O}(\phi) = \int_\mathcal{M} d\mu\, p(\phi)\,\mathcal{O}(\varphi)$ holds exactly.

## Connections

- [[physics-informed-rg-flows]] — the PIRG framework provides the differential equations that PIKs solve.
- [[wegner-equation]] — the key PDE that determines the kernel $\dot{\phi}_t$.
- [[sign-problem]] — complex PIK-architecture solves sign problems by constructing PIKfolds.
- [[contour-deformation]] — PIKfold is a specific, analytically constructed deformation of the integration manifold.
- [[lefschetz-thimbles]] — PIKfolds generalise thimbles; exist in cases where thimbles fail.

## Open questions

- Efficient parametrisation of $S_t$ for strongly interacting theories.
- Global optimisation of the PIRG pair $(F_t, \dot{\phi}_t)$ beyond local layer-wise optimisation.
- Extension to gauge theories (non-Abelian link fields, compact configuration space).
- Scaling with lattice volume.
