---
title: "Physics-informed renormalisation group flows (PIRGs)"
tags: [pirg, renormalisation-group, wegner-equation, generative-network]
sources: ["[[ihssen-2025-generative-sampling-piks]]", "[[ihssen-2026-solving-sign-problems-piks]]"]
---

## Definition

Physics-informed renormalisation group (PIRG) flows are RG flows in the space of distributions governed by the equation $\frac{d p_t(\phi)}{dt} = -\frac{\partial}{\partial\phi}[\Psi_t(\phi)p_t(\phi)]$ with $\Psi_t(\phi) \propto \partial S_t(\phi)/\partial\phi$. For the specific choice $\Psi_t = \dot{\phi}_t$ (the PIK), this becomes the Wegner equation. PIRGs describe scale-dependent transformations of field configurations that are used in the PIK architecture to build generative flows from a simple base distribution to a complex target.

## Intuition

Standard RG flows evolve the action as the UV cutoff is integrated out. PIRG flows instead evolve the field coordinates $\phi$ themselves, interpreting the RG-time $t$ as a parameter interpolating between a simple Gaussian theory (at $t=0$) and the full interacting theory (at $t=1$). The flow is invertible (the PIK $\dot{\phi}_t$ generates a diffeomorphism), so the partition function $Z$ is conserved along the flow. The physics — encoded in the action $S_t(\phi)$ — is fully specified analytically at every $t$, which is what makes the kernels computable without training.

## Key results

- The PIRG pair $(S_t, \dot{\phi}_t)$ with $p_t(\phi) = e^{-S_t(\phi)}$ satisfies $\frac{dS_t}{dt} + \dot{\phi}_t\frac{\partial S_t}{\partial\phi} = \frac{\partial\dot{\phi}_t}{\partial\phi}$ (Wegner equation).
- Total derivative structure: $\frac{dZ}{dt} = 0$ (partition function conserved).
- PIRG flows reduce to Fokker-Planck dynamics (Langevin) when $\Psi_t \propto \partial S_t/\partial\phi$ with RG-time-dependent classical action.
- The Polchinski equation for the Wilsonian effective action is recovered as a special case.
- PIRG flows provide a general framework for invertible field transformations; the PIK is the map that the flow generates.

## Connections

- [[physics-informed-kernels]] — PIKs are the infinitesimal generators of PIRG flows.
- [[wegner-equation]] — the specific differential equation governing PIRG flows.
- [[sign-problem]] — PIRG-based PIK architecture solves the sign problem for complex distributions.

## Open questions

- Whether PIRG flows with other choices of $\Psi_t$ (beyond the PIK/Wegner choice) offer advantages.
- Relationship to functional renormalisation group (fRG) and the Wetterich equation.
