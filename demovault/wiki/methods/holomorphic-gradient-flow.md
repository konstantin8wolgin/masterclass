---
title: "Holomorphic gradient flow"
tags: [holomorphic-flow, contour-deformation, lefschetz-thimbles, sign-problem]
sources: ["[[lawrence-2023-deep-learning-fermion-sign]]", "[[aarts-2026-lattice-sign-problem-review]]"]
---

## Definition

The holomorphic gradient flow evolves a field configuration $A_\mu(x) \in \mathbb{R}$ into the complexification according to $\frac{dA_\mu(x)}{dt} = -(\partial S/\partial A_\mu(x))^*$, generating a one-parameter family of integration contours $\mathcal{M}_t$ in the complexified field space. At finite flow time $t$, the contour $\mathcal{M}_t$ is a deformed integration surface with improved sign problem properties; as $t \to \infty$, it approaches Lefschetz thimbles.

## Intuition

The gradient flow moves field configurations in the direction that most rapidly decreases $\text{Re}\, S$ (the real part of the action) in the complexified space. Along the flow, $\text{Im}\, S$ is conserved: $\frac{d}{dt}\text{Im}\, S[\phi(t)] = 0$. This means the contour at any finite flow time has the same $\text{Im}\, S$ as the starting configuration. As $t \to \infty$, the flow converges to Lefschetz thimbles where $\text{Im}\, S = \text{const}$ exactly. At finite $t$, the contour is a smooth deformation of the real plane, easier to work with than the full thimble.

## Algorithm

1. Start at $A_\mu(x) \in \mathbb{R}$ (real configuration).
2. Integrate $\frac{dA_\mu}{dt} = -(\partial S/\partial A_\mu)^*$ numerically (e.g. Euler or Runge-Kutta).
3. At chosen flow time $t$, the flowed configuration $A_\mu(t) \in \mathbb{C}$ defines a point on $\mathcal{M}_t$.
4. Monte Carlo sampling along $\mathcal{M}_t$ uses the effective action $S_\text{eff} = S[\tilde{\phi}(\phi)] - \log J(\phi)$ where $J = \det(\partial\tilde{\phi}/\partial\phi)$.

For the Thirring model in $1+1$D, the holomorphic contours studied earliest were defined by specific Fourier series representations of $\text{Im}\, A_0(x)$.

## Key results

- Every contour deformation parameterised by the real plane is a subtraction: $g_\text{sub} = e^{-S(\tilde{\phi})}\det(\partial\tilde{\phi}/\partial\phi) - e^{-S(\phi)}$ ([[lawrence-2023-deep-learning-fermion-sign]], Eq. 16).
- At small $t$: holomorphic flow matches the gradient descent subtraction to first order in $u$.
- Long flow times found empirically effective for alleviating the sign problem in Thirring model [Alexandru et al.].
- For $U(1)$ gauge fields, optimal contours of the form $\text{Im}\, A_0(x) = \sum_n c_n \cos(n\,\text{Re}\, A_0(x))$ outperform gradient flow at the same computational cost.

## Connections

- [[lefschetz-thimbles]] — limit of holomorphic gradient flow at $t\to\infty$.
- [[contour-deformation]] — the general framework; gradient flow is one construction.
- [[subtraction-method]] — every contour deformation corresponds to a subtraction.
- [[sign-problem]] — the target.

## Open questions

- Optimal flow time: too short gives little improvement; too long approaches thimbles with their multi-thimble sign problem.
- Extension to non-Abelian gauge theories (SU(3) links).
