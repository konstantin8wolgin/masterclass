---
title: "Gaussian mixture decomposition of fermionic partition function"
tags: [gaussian-mixture, hubbard-model, fermionic, sign-problem, analytic-integration]
sources: [
  "[[ihssen-2026-fermionic-sign-piks]]",
  "[[fukushima-2026-one-site-hubbard-piks]]"
]
---

## Definition

The Gaussian mixture decomposition is the explicit expansion of a fermionic partition function (after [[hubbard-stratonovich-transformation]]) as a finite sum of shifted Gaussians in the auxiliary field $\phi$. Each term has the form $e^{-a\phi^2 + b\phi + c}$ with $b \in i\mathbb{R}$ (purely imaginary shift), and is integrable in closed form via:
$$\int_{-\infty}^{+\infty} dx\, e^{-ax^2 + bx + c} = \sqrt{\frac{\pi}{a}}\, e^{\frac{b^2}{4a} + c}$$

The decomposition arises from the finiteness of the fermionic Fock space: for $N_s$ spin states, the trace over fermions at one site yields $N_s + 1$ Gaussian terms (here $N_s = 2$ gives three terms).

## Intuition

The one-site Hubbard model has only four Fock states: $|0\rangle$, $|\uparrow\rangle$, $|\downarrow\rangle$, $|\uparrow\downarrow\rangle$. The trace $\mathrm{Tr}[e^{-\beta H}]$ is therefore a sum of four terms. After the Hubbard-Stratonovich transformation introduces the auxiliary field $\phi$, the factor $(1 + e^{\beta(i\phi+\mu+U/2)})^2$ arises from the spin-up and spin-down fermion traces; expanding this square yields exactly three distinct Gaussian weights (since $|\uparrow\rangle$ and $|\downarrow\rangle$ are degenerate). The imaginary shifts $j\cdot iU$ in the Gaussian centers ($j = 0, 1, 2$) directly encode the fermion occupation number. In a multi-site model the fermion determinant replaces this with a product structure, exponentially increasing the number of contributing Gaussians.

## Key results

For the one-site Hubbard model, expanding the squared bracket in the path integral gives ([[fukushima-2026-one-site-hubbard-piks]], eq. 6):

$$Z = Z_0 + Z_1 + Z_2$$

- $Z_0 \propto \int d\phi\, e^{-\frac{\beta}{2U}\phi^2}$: free Gaussian, weight 1 (empty state)
- $Z_1 \propto \int d\phi\, e^{-\frac{\beta}{2U}\phi^2 + i\beta\phi + \beta(\mu+U/2)}$: Gaussian shifted by $+i\beta$ in exponent, weight $2e^{\beta(\mu+U/2)}$ (two single-occupation spin states)
- $Z_2 \propto \int d\phi\, e^{-\frac{\beta}{2U}\phi^2 + 2i\beta\phi + 2\beta(\mu+U/2)}$: Gaussian shifted by $+2i\beta$, weight $e^{2\beta(\mu+U/2)}$ (double occupation)

The same decomposition appears in [[ihssen-2026-fermionic-sign-piks]] (eq. 4):
$$\mathcal{N} = \sqrt{\frac{\beta}{2\pi U}}\int d\phi\left[e^{-\frac{\beta}{2U}\phi^2} + 2e^{\beta(\mu+U/2+i\phi)}e^{-\frac{\beta}{2U}\phi^2} + e^{2\beta(\mu+U/2+i\phi)}e^{-\frac{\beta}{2U}\phi^2}\right]$$

**Occupation statistics from the weights:**
- Relative weight of $Z_1$ to $Z_0$: $2e^{\beta\mu}\cdot e^{\beta U/2}$ (chemical potential + exchange energy)
- Relative weight of $Z_2$ to $Z_0$: $e^{\beta(2\mu-U)}$ (double occupation penalized by $-U$ in the exponent)
- At half-filling $\mu = -U/2$: $Z_0 = Z_1/2 = Z_2$ (equal weights; no preference for occupation)

**Source of the sign problem:** The imaginary shifts $j \cdot iU$ in the Gaussian exponents make the integrand complex for real $\phi$. The distribution is thus a mixture of complex Gaussians, not a positive measure.

**CDF approach:** Because each component is a Gaussian, its CDF is an error function, which can be analytically continued to $\mathbb{C}$. The full CDF decomposes as $F_t = \sum_i w_i F_{t,i}$, enabling the analytic PIK velocity formula $\dot\phi_t = -(1/p_t)(d/dt)F_t$ ([[ihssen-2026-fermionic-sign-piks]], eq. 14).

## Connections

- [[gaussian-composition-effective-action]] — the $t$-dependent version of this decomposition; provides the analytic PIK $\dot\phi = iU$.
- [[hubbard-stratonovich-transformation]] — the transformation that introduces $\phi$ and produces the Gaussian structure.
- [[one-site-hubbard-fermionic-action]] — the action whose expansion yields these Gaussians.
- [[physics-informed-kernels]] — the PIKfold Mixture approach applies one PIK per Gaussian component; CDF reformulation gives analytic kernel.
- [[sign-problem]] — the imaginary Gaussian shifts are the direct source of the sign problem.
- [[average-phase]] — $Q_{sp} = Z/Z_{pq}$ is computed from these weights; its drop signals severe sign problem.

## Open questions

- Whether a tractable finite-dimensional Gaussian decomposition survives in multi-site models after partial integration over some sites.
- Whether exact results for sign problem severity as a function of $\beta, \mu, U$ can be derived purely from the mixture weights.
