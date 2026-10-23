---
title: "PIKfold ripping"
tags: [pikfold, sign-problem, physics-informed-kernels, singularity, rg-flow]
sources: [
  "[[fukushima-2026-one-site-hubbard-piks]]",
  "[[ihssen-2026-fermionic-sign-piks]]"
]
---

## Definition

PIKfold ripping is the failure mode of a PIKfold $\mathcal{M}_t$ that occurs when a pole $\phi_\star(t)$ of the complexified action $S_t(\phi)$ crosses the sampling manifold $\mathcal{M}_t$ during RG-time evolution $t \in [0,1]$. At the crossing point, the imaginary part of $p_t = e^{-S_t}$ diverges on $\mathcal{M}_t$, making the PIKfold ODE $d\phi_I/d\phi_R = -\mathrm{Im}[p_t]/\mathrm{Re}[p_t]$ ill-defined. The manifold can no longer be continuously extended — it is said to have ripped open.

## Intuition

The PIKfold is a real curve $\phi_I(\phi_R)$ in the complex $\phi$-plane satisfying a weight-preservation condition ([[ihssen-2026-solving-sign-problems-piks]], eq. 16). This curve must avoid the poles of the action, which function as branch points in the complexified field space. During RG-time evolution the action changes, and so do the pole positions $\phi_\star(t)$. If the trajectory of a pole in $\mathbb{C}$ intersects the trajectory of the manifold in $(t, \mathbb{C})$-space, the ODE for the manifold encounters a singularity and the manifold tears. This is analogous to a branch cut crossing in analytic continuation or the topological obstruction encountered when Lefschetz thimbles collide.

## Key results

For the one-site Hubbard model, poles sit at ([[fukushima-2026-one-site-hubbard-piks]], eq. 8):
$$\phi_\star^{(n)} = \frac{(2n+1)\pi}{\beta} + i\!\left(\mu + \frac{U}{2}\right), \qquad n \in \mathbb{Z}$$

Two choices of RG-time path give different pole dynamics:

**Path 1** (eq. 9a of [[fukushima-2026-one-site-hubbard-piks]]):
$$S_t(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + e^{t\beta(i\phi+\mu+U/2)}\right]$$
Poles: $\mathrm{Re}(\phi_\star) = (2n+1)\pi/(t\beta)$, $\mathrm{Im}(\phi_\star) = \mu+U/2$ — poles begin at $|\mathrm{Re}(\phi_\star)| = \infty$ at $t=0$ and move inward. Ripping risk arises for $\mu \gtrsim 0.1$ (near half-filling) where the inward-moving poles reach the manifold.

**Path 2** (eq. 10a of [[fukushima-2026-one-site-hubbard-piks]]):
$$S_t(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + e^{\beta(i\phi+\mu+tU/2)}\right]$$
Poles: $\mathrm{Re}(\phi_\star) = (2n+1)\pi/\beta$, $\mathrm{Im}(\phi_\star) = \mu+tU/2$ — poles move upward in the imaginary direction as $t$ increases. This path performs better for $\mu < -U/2$ (where the manifold stays safely below the poles) but rips for $\mu \approx 0.2$.

**Remedies:**
1. **Two-step PIK-pair resolution** ([[fukushima-2026-one-site-hubbard-piks]], Sec. VI): find an intermediate $p_t$ via $\mathcal{T}_t = \log\sum_i c_{i,t}\rho_{i,t}(\phi)$ (log-sum-exp Gaussian mixture) that keeps $\mathcal{M}_t$ away from all poles at every $t$.
2. **Complex initial manifold** $\mathcal{M}_0 \neq \mathbb{R}$: start on a manifold already shifted above the relevant poles.
3. **PIKfold Mixture** ([[ihssen-2026-fermionic-sign-piks]]): distribute $p$ across Gaussian components, apply a separate PIK per component that avoids the poles relevant to that component.

**Naive diagnostic for ripping:** large imaginary-to-real ratios $|\mathrm{Im}[p_t(\phi^{(i)})]|/|\mathrm{Re}[p_t(\phi^{(i)})]| \to \infty$ at the peak positions $\phi^{(i)}$ of the distribution as $t \to t_\mathrm{rip}$.

## Connections

- [[physics-informed-kernels]] — PIKfold ripping is a failure mode specific to the complex PIK architecture.
- [[physics-informed-rg-flows]] — the two-step resolution framework addresses ripping within the PIRG context.
- [[lefschetz-thimbles]] — ripping is the PIK analogue of the inter-thimble phase problem; LT encounters topological obstructions when Stokes phenomena occur.
- [[one-site-hubbard-fermionic-action]] — provides the concrete pole positions that cause ripping in the Hubbard model.
- [[sign-problem]] — ripping is a structural obstacle to PIKfold-based resolution of the fermionic sign problem.
- [[contour-deformation]] — any contour deformation must avoid poles; ripping signals that the chosen deformation path is topologically incompatible.

## Open questions

- General topological criteria for predicting ripping from the pole structure of $S(\phi)$ without explicit RG-time integration.
- Optimal choice of RG-time path $S_t$ to maximize the distance between poles and the sampling manifold across the $(\beta, \mu, U)$ parameter space.
- Extension to multi-site models where the pole structure involves a fermion determinant vanishing on a codimension-1 manifold in field space.
