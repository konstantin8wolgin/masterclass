---
title: "Partition function (lattice field theory)"
tags: [partition-function, path-integral, lattice-qft]
sources: [
  "[[temmen-2025-fully-ergodic-radial-updates]]",
  "[[aarts-2026-lattice-sign-problem-review]]",
  "[[ihssen-2025-generative-sampling-piks]]"
]
---

## Statement

The partition function of a lattice field theory is the path integral

$$Z = \int \mathcal{D}\phi\, e^{-S[\phi]} = Z^{-1}\int \mathcal{D}\phi\, \mathcal{O}[\phi]\, e^{-S[\phi]}$$

Expectation values of observables are:

$$\langle\mathcal{O}\rangle = \frac{\int\mathcal{D}\phi\, \mathcal{O}[\phi]\, e^{-S[\phi]}}{\int\mathcal{D}\phi\, e^{-S[\phi]}}$$

In HMC, this is approximated by sampling $\phi^{(i)} \sim e^{-S}/Z$ and averaging:

$$\langle\mathcal{O}\rangle \approx \frac{1}{N_\text{conf}}\sum_{i=1}^{N_\text{conf}}\mathcal{O}[\phi^{(i)}]$$

## Variables

- $\phi$: the field (scalar, gauge, fermion auxiliary, etc.)
- $\mathcal{D}\phi = \prod_x d\phi(x)$: path integral measure (product over lattice sites)
- $S[\phi]$: the Euclidean action (real for standard importance sampling; complex at non-zero chemical potential)
- $\mathcal{O}[\phi]$: observable
- $N_\text{conf}$: number of Monte Carlo configurations

## Derivation sketch

The Euclidean path integral arises from Wick-rotating the Minkowski time integral $t \to -i\tau$ with $\tau \in [0, \beta]$ (where $\beta = 1/T$ is the inverse temperature). The resulting integral over Euclidean field configurations has a real, positive weight $e^{-S_E}$ for theories without a chemical potential. Discretising the spacetime on a lattice with spacing $a$ gives the lattice path integral.

## Limits and special cases

- At $\mu = 0$ (zero chemical potential): $S$ is real and $e^{-S} \geq 0$; standard importance sampling applies.
- At $\mu \neq 0$: for fermionic theories, $\det M(U;\mu) \in \mathbb{C}$ makes $e^{-S}$ complex — the sign problem.
- Large volume limit ($V \to \infty$): the partition function is dominated by the saddle point of $S$; the sign problem SNR decays as $e^{-V\Delta f}$ where $\Delta f = f - f_{pq}$ is the free energy density difference.
- Phase-quenched: replace $e^{-S}$ by $|e^{-S}|$; tractable but gives wrong physics at $\mu \neq 0$.

## Sources

- [[temmen-2025-fully-ergodic-radial-updates]] (Eq. 1–2: HMC partition function)
- [[aarts-2026-lattice-sign-problem-review]] (Eq. 1–2: lattice QCD partition function)
- [[ihssen-2025-generative-sampling-piks]] (Eq. 1: general normalised distribution)
