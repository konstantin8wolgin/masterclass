---
title: "Subtraction method for the sign problem"
tags: [subtraction, sign-problem, vector-field, neural-network]
sources: ["[[lawrence-2023-deep-learning-fermion-sign]]"]
---

## Definition

The subtraction method alleviates the fermion sign problem by adding a function $g(\phi)$ that integrates to zero — a subtraction — to the Boltzmann factor: $e^{-S_\text{sub}} = e^{-S} - g$. Physical expectation values are unchanged (since $\int g = 0$), but the phase fluctuations of the modified Boltzmann factor can be much smaller. The standard construction uses $g = \nabla \cdot v$ for a vector field $v$, guaranteed to integrate to zero by the divergence theorem.

## Intuition

The sign problem arises because $e^{-S}$ has complex values that cancel. Adding a zero-integral function doesn't change the total, but it can redistribute the cancellations: we can try to choose $g$ such that $e^{-S} - g$ has a nearly constant phase. The divergence construction gives a large family of such functions. A neural network can parametrise the vector field $v$ and be optimised to maximise the average phase of $e^{-S_\text{sub}}$.

## Algorithm

**Setup:** Write $v = e^{-S}u$ where $u(\phi)$ is a vector field of order 1. The subtracted action is:

$$e^{-S_u} = e^{-S}(1 + u \cdot \nabla S - \nabla \cdot u)$$

**Analytic subtractions (no ML):**

- Heavy-dense expansion: at large $\mu$, expand $\det K$ via polymer representation. Dominant term gives an explicit analytic $g_\text{hd}$.
- Gradient descent: set $u \propto \frac{\delta Z_{Q,u}}{\delta u} = -i\nabla\text{Im}\, S$; take one gradient step.
- Heat equation: evolve $e^{-S_t}$ under $\partial_t e^{-S_t} = -\sum_i \partial^2_{\phi_i^2}e^{-S_t}$; long-time limit removes all phase fluctuations (but computationally intractable).

**Neural network subtraction:**

1. Parametrise $u(\phi)$ as a deep neural network.
2. Train to maximise the average phase $\langle e^{i\text{Im}\, S_u} \rangle$ of the subtracted ensemble.
3. At inference: sample configurations from $e^{-S_u}$ (importance sampling), compute observables using the estimator:
$$\langle\mathcal{O}\rangle = \left\langle \mathcal{O} - \frac{e^{-S_u}}{e^{-S_v}} u \cdot \nabla\mathcal{O} \right\rangle_{S_u}$$

## Key results

- Perfect subtraction always exists: $g_\text{exact}(\phi) = e^{-S(\phi)} - \frac{\int \mathcal{D}\phi'\, e^{-S(\phi')}}{\int \mathcal{D}\phi'}$ entirely removes the sign problem.
- Every infinitesimal contour deformation is a subtraction (to first order); not every subtraction is a contour deformation.
- The scaling trick ($v = e^{-S}u$) breaks the guarantee that a good $u$ gives a perfect subtraction when $e^{-S}$ has zeros.
- Estimator Eq. (15) of [[lawrence-2023-deep-learning-fermion-sign]] empirically avoids severe sign-to-noise without the naive ratio's exponential variance.
- Demonstrated on staggered Thirring model in $1+1$D: neural network subtraction significantly improves average phase compared to analytic subtractions at large $\mu$.

## Connections

- [[sign-problem]] — the problem being addressed.
- [[average-phase]] — the diagnostic optimised during training.
- [[contour-deformation]] — every contour deformation is a subtraction.
- [[holomorphic-gradient-flow]] — generates subtractions at small flow time.
- [[thirring-model-action]] — testbed model.

## Open questions

- Why does estimator Eq. (15) empirically avoid sign-to-noise — is there a rigorous argument?
- Generalisation to lattice QCD where $u$ must respect gauge symmetry.
