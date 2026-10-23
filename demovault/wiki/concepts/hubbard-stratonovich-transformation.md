---
title: "Hubbard-Stratonovich transformation"
tags: [hubbard-stratonovich, path-integral, auxiliary-field, sign-problem, fermionic]
sources: [
  "[[fukushima-2026-one-site-hubbard-piks]]",
  "[[ihssen-2026-fermionic-sign-piks]]",
  "[[temmen-2025-fully-ergodic-radial-updates]]"
]
---

## Definition

The Hubbard-Stratonovich (HS) transformation is an exact Gaussian integral identity that decouples a quartic fermion interaction into a quadratic fermionic coupling to a new real bosonic auxiliary field $\phi$. The identity is:
$$e^{-\frac{\Delta_t U}{2}(\hat{n}_\uparrow - \hat{n}_\downarrow)^2} \propto \int_{-\infty}^{\infty} d\phi\, e^{-\frac{1}{2U\Delta_t}\phi^2 - i\phi(\hat{n}_\uparrow - \hat{n}_\downarrow)}$$

After integration over $\phi$ on the right-hand side, the original quartic operator is exactly reproduced. Before integration, the fermionic operator appears linearly in $\phi$, making the fermion trace (determinant) tractable.

## Intuition

The HS transformation trades a hard many-body problem (quartic fermion-fermion interaction in the Hubbard Hamiltonian) for a fluctuating-field problem: sampling over the auxiliary bosonic field $\phi$. The price is that $\phi$ couples to fermions via the factor $e^{i\phi\hat{n}}$, which is complex — this is the direct cause of the sign problem in the resulting bosonic path integral. The transformation is exact; no approximation is involved. The choice of decoupling channel (spin, density, or pairing) affects the form of $\phi$ and the severity of the resulting sign problem.

## Key results

**Zero-dimensional (one-site, one-time-slice) limit**: applied to $H = U\hat{n}_\uparrow\hat{n}_\downarrow - \mu\hat{n}$ with $\Delta_t = \beta$ and $\kappa = 0$ (no hopping), the HS transformation gives ([[ihssen-2026-fermionic-sign-piks]], eq. 3; [[fukushima-2026-one-site-hubbard-piks]], eq. 7):
$$Z = \sqrt{\frac{\beta}{2\pi U}}\int_{-\infty}^\infty d\phi\, e^{-S(\phi)}, \qquad S(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + e^{\beta(\mu + i\phi + U/2)}\right]$$

The logarithm term arises from $\mathrm{Tr}_\uparrow[1 + e^{\beta(\mu+i\phi+U/2)}]\cdot\mathrm{Tr}_\downarrow[\ldots]$ where each trace evaluates to $1 + e^{(\ldots)}$. The imaginary coupling $i\phi$ in the exponent makes $S(\phi) \in \mathbb{C}$ for $\phi \in \mathbb{R}$. See [[one-site-hubbard-fermionic-action]].

**Lattice (multi-site, finite-temperature) case**: applied via Suzuki-Trotter decomposition to the full Hubbard Hamiltonian on a lattice with $N_t$ time slices and hopping $\kappa$ ([[temmen-2025-fully-ergodic-radial-updates]], [[hubbard-action]]):
$$S[\phi] = \frac{1}{2U\Delta_t}\sum_{t,x}\phi_{tx}^2 - \log\det M[\phi|\kappa]\det M[-\phi|-\kappa]$$

Here $\det M[\phi|\kappa]$ is the fermion determinant encoding hopping. At $\mu = 0$ (half-filling), particle-hole symmetry guarantees $\det M[\phi|\kappa]\det M[-\phi|-\kappa] > 0$ — no sign problem. At $\mu \neq 0$ the determinant becomes complex.

## Algorithm section: standard implementation steps

1. Write the Hubbard interaction $U\hat{n}_\uparrow\hat{n}_\downarrow$ using the identity $\hat{n}_\uparrow\hat{n}_\downarrow = \frac{1}{4}(\hat{n}_\uparrow + \hat{n}_\downarrow)^2 - \frac{1}{4}(\hat{n}_\uparrow - \hat{n}_\downarrow)^2$ to choose a decoupling channel.
2. Apply the Suzuki-Trotter approximation: $e^{-\beta H} \approx (e^{-\Delta_t H})^{N_t}$ with $\Delta_t = \beta/N_t$.
3. Apply the HS identity to the interaction part of $e^{-\Delta_t H_\mathrm{int}}$ at each time slice $t$.
4. Integrate out the fermion fields $c, c^\dagger$ analytically — this produces the fermion determinant $\det M[\phi]$.
5. The resulting bosonic path integral $Z = \int \mathcal{D}\phi\, e^{-S[\phi]}$ has a complex action for $\mu \neq 0$ → sign problem.

## Connections

- [[one-site-hubbard-fermionic-action]] — the action produced by the HS transformation in the zero-dimensional limit.
- [[hubbard-action]] — the full multi-site action arising from the HS transformation on the lattice.
- [[sign-problem]] — the imaginary coupling $i\phi\hat{n}$ introduced by the HS transformation is the direct cause of the sign problem.
- [[gaussian-mixture-decomposition]] — the HS-transformed partition function decomposes as a Gaussian mixture in $\phi$.
- [[physics-informed-kernels]] — PIKs operate on the bosonic field $\phi$ introduced by the HS transformation.

## Open questions

- Whether discrete HS transformations (as opposed to continuous) offer advantages for the sign problem structure, e.g. reducing the number of poles.
- Connection between the choice of HS decoupling channel (spin vs. density vs. pairing) and the severity and structure of the resulting sign problem.
- Extension to non-Abelian gauge theories where the HS field is a matrix-valued field.
