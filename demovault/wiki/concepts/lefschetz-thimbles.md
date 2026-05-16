---
title: "Lefschetz thimbles"
tags: [lefschetz-thimbles, sign-problem, complex-analysis, holomorphic]
sources: ["[[aarts-2026-lattice-sign-problem-review]]", "[[ihssen-2026-solving-sign-problems-piks]]"]
---

## Definition

Lefschetz thimbles are specific integration contours in the complexified field space $\mathbb{C}^N$ attached to critical points (saddle points) of the action. On a stable thimble $\mathcal{J}_k$, the imaginary part of the action $\text{Im}\, S$ is constant, so the oscillatory phase $e^{-i\text{Im}\, S}$ does not contribute to cancellations along the contour. The original path integral can be rewritten as a sum over thimbles by Picard-Lefschetz theory.

## Intuition

The sign problem comes from rapid phase oscillations of $e^{-S}$. Lefschetz thimbles are the steepest-descent contours for $\text{Re}\, S$; along them, $\text{Re}\, S$ is maximal at the critical point and decreases away from it (the integrand is peaked), while $\text{Im}\, S$ stays constant (no oscillation). Deforming the integration to a thimble removes the intra-thimble sign problem. However, if multiple thimbles contribute, their relative phases give a global sign problem, and the Jacobian of the thimble parametrisation gives a residual sign problem.

## Key results

- Stable thimbles $\mathcal{J}_k$ are defined by the flow equations $\partial_t \tilde{\phi}(x,t) = -\overline{(\partial S/\partial \tilde{\phi})}$ with $\lim_{t\to\infty}\tilde{\phi} = \tilde{\phi}_0^{(k)}$ (critical point).
- The original path integral equals $Z = \sum_k n_k e^{-i\text{Im}\, S(\tilde{\phi}_0^{(k)})} \int_{\mathcal{J}_k} \mathcal{D}\tilde{\phi}\, e^{-\text{Re}\, S[\tilde{\phi}]}$ where $n_k$ is the intersection number with the unstable thimble.
- Two sign problems remain: (i) global — relative phases $e^{-i\text{Im}\, S(\tilde{\phi}_0^{(k)})}$ between contributing thimbles; (ii) residual — Jacobian $J(\phi) = \det(\partial\tilde{\phi}/\partial\phi)$ is complex.
- Holomorphic gradient flow (same equations, finite flow time) gives a practical family of contours without the need to enumerate all critical points.
- PIKfolds ([[physics-informed-kernels]]) generalise thimbles and exist in cases where thimbles fail (e.g. the Minkowski oscillator).

## Connections

- [[contour-deformation]] — thimbles are a specific type of contour deformation.
- [[holomorphic-gradient-flow]] — the flow that defines thimbles (at $t \to \infty$) and practical deformed contours (at finite $t$).
- [[sign-problem]] — thimbles reduce but do not fully solve the sign problem.
- [[physics-informed-kernels]] — PIKfolds as an alternative that avoids thimble limitations.
- [[complex-langevin-dynamics]] — the other main holomorphic extension approach.

## Open questions

- Efficient algorithms to find all relevant critical points and compute intersection numbers for large lattice volumes.
- Whether machine learning can find better approximations to thimbles than the holomorphic gradient flow.
