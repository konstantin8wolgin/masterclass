---
title: "RBF expansion of the path integral"
tags: [rbf, radial-basis-function, path-integral, sign-problem, finite-density]
sources: ["[[balassa-2026-sign-problem-rbf-networks]]"]
---

## Definition

The RBF (radial basis function) path integral method approximates the nonlinear interaction part of the Euclidean action's Boltzmann factor $e^{-S_{1,x}(m,\lambda)}$ as a sum of $K$ Gaussian kernels centred at positions $c_k$ in field space, making the total path integral analytically tractable by diagonalisation in momentum space.

## Intuition

The path integral $Z = \int \mathcal{D}\phi\, e^{-S}$ is hard to evaluate when $S$ contains self-interaction terms (e.g. $\phi^4$) because the resulting integral does not factorise. If the interaction Boltzmann factor can be written as a sum of Gaussians, each term becomes a Gaussian integral — exactly solvable in momentum space. The momentum-space diagonalisation then reduces the complexity from $\mathcal{O}(K^{N^4})$ (all kernel combinations across all sites) to $\mathcal{O}(KN^4)$ by expressing the partition function as a product over momenta.

## Algorithm

1. **Split the action:** $S = S_{0,x}(\mu) + S_{1,x}(m,\lambda)$ where $S_0$ contains $\mu$-dependent and kinetic terms (quadratic in $\phi$); $S_1$ contains $\phi^4$ self-interaction.

2. **RBF expand the interaction:**
$$e^{-S_{1,x}(m,\lambda)} \approx \sum_{k=1}^K a_k e^{-A\sum_a (\phi_{a,x} - c_{a,k})^2}$$
where $a_k$ are weights determined by least squares against training field configurations; $c_{a,k}$ are pre-chosen centers; $A$ is the width.

3. **Write the partition function** as a sum over all $K^{N^4}$ combinations of kernel choices at each site, each of which is Gaussian.

4. **Transform to momentum space** using unitary $U$ that diagonalises $W = [M(\mu) - 2\bar{A}]$:
$$Z \to \prod_p \sum_{k=1}^K Q_k[\tilde{\Phi}_{a,p}, \mu, A, c_{a,k}, a_k]$$
Complexity: $\mathcal{O}(KN^4)$.

5. **Compute observables** (e.g. $\langle n\rangle = \frac{1}{V}\frac{\partial\ln Z^\text{rbf}}{\partial\mu}$) analytically from the approximated $Z^\text{rbf}$.

## Key results

- Relative error $E_R = \frac{|\ln Z_U - \ln Z_0|}{|\ln Z_U|}$ stays below a few percent for symmetric center distributions around 0.
- Fixed center grid ($K=25$, $A=1$, $\Delta c = 0.5$) achieves smaller error variance than random centers, though both are below 1%.
- Phase transition at $\mu_c = 1.17\pm 0.018$ extracted from the effective potential.
- Silver blaze: $\langle n\rangle \approx 0$ for $\mu < \mu_c$ correctly reproduced.

## Connections

- [[sign-problem]] — overcome by making the integral analytically tractable.
- [[silver-blaze-problem]] — correctly reproduced by the method.
- [[partition-function-lattice]] — the quantity approximated.
- [[complex-langevin-dynamics]] — compared against for validation.

## Open questions

- Extension to fermionic theories.
- Scaling of required $K$ with coupling strength $\lambda$ and lattice volume.
- Error bounds when $e^{-S_1}$ is not well approximated by Gaussians (e.g. strong coupling).
