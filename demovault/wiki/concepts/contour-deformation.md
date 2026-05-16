---
title: "Contour deformation"
tags: [contour-deformation, sign-problem, holomorphic, cauchy]
sources: ["[[lawrence-2023-deep-learning-fermion-sign]]", "[[aarts-2026-lattice-sign-problem-review]]"]
---

## Definition

A contour deformation replaces the original integration contour (real field space $\mathbb{R}^N$) with a deformed contour $\mathcal{M} \subset \mathbb{C}^N$ in the complexification. By Cauchy's theorem, the value of the path integral is unchanged (if the integrand is holomorphic and boundary terms vanish), but the integrand on $\mathcal{M}$ may have smaller phase fluctuations and hence a milder sign problem.

## Intuition

The original integration is performed along the real axis; the Boltzmann factor $e^{-S(\phi)}$ oscillates wildly there. By deforming into the complex plane — to a contour where $\text{Im}\, S$ is approximately constant — the oscillations are reduced. At the extreme, on a Lefschetz thimble, $\text{Im}\, S = \text{const}$ exactly. The price is that the deformed contour must be found, sampled along, and the Jacobian of the map must be computed.

## Key results

- Every contour deformation parameterised by the real plane ($\tilde{\phi}: \mathbb{R}^N \to \mathbb{C}^N$) corresponds to a subtraction: $g_\text{sub}(\phi) = e^{-S(\tilde{\phi}(\phi))}\det\frac{\partial\tilde{\phi}}{\partial\phi} - e^{-S(\phi)}$ ([[lawrence-2023-deep-learning-fermion-sign]], Eq. 16).
- Converse is false: not all subtractions correspond to contour deformations (e.g. $S = \cos\theta + \epsilon$ admits a perfect subtraction but no contour deformation that entirely removes the sign problem).
- Holomorphic gradient flow defines a one-parameter family of contours; at small flow times it improves the average phase.
- ML-based contour deformations train a neural network to find a good deformation map, and have been applied to scalar theories and QCD.
- Contour deformations work for all theories (unlike thimble methods which face the multi-thimble sign problem for certain theories).

## Connections

- [[lefschetz-thimbles]] — Lefschetz thimbles are the ideal contours where $\text{Im}\, S = \text{const}$.
- [[holomorphic-gradient-flow]] — constructs contours by following the gradient flow of $S$ in the complex plane.
- [[subtraction-method]] — every contour deformation is a subtraction; not vice versa.
- [[sign-problem]] — the target: reduce phase fluctuations.
- [[physics-informed-kernels]] — PIKfolds are analytically constructed integration manifolds; related concept.

## Open questions

- Extension of ML contour deformation from scalar field theories to full lattice QCD (gauge fields, fermion determinant).
- Computational cost of evaluating the Jacobian $\det(\partial\tilde{\phi}/\partial\phi)$ for large lattice volumes.
