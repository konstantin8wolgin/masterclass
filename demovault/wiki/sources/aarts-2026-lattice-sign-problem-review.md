---
title: "Lattice field theories with a sign problem"
authors: [Gert Aarts, Dénes Sexty]
year: 2026
type: paper
tags: [review, sign-problem, lattice-qcd, lefschetz-thimbles, complex-langevin, holomorphic-flow, contour-deformation, machine-learning]
source_file: raw/papers/2604.24290v1.pdf
---

## Summary

Comprehensive review of the sign problem in lattice field theory, focussing on QCD at non-zero baryon chemical potential and methods which may allow exploration of the full QCD phase diagram. Covers the standard formulation, the overlap problem and silver blaze phenomenon, importance-sampling based methods (Taylor expansion, imaginary chemical potential, reweighting), and the main class of approaches: holomorphic extensions — Lefschetz thimbles, holomorphic flow, contour deformations, and complex Langevin dynamics. Machine learning contributions to each approach are highlighted throughout. The review also covers dual variables, tensor networks, and the introduction of new degrees of freedom.

## Key claims

- The fermion determinant $\det M(U;\mu)$ is complex at non-zero baryon chemical potential $\mu$; its complex conjugate satisfies $[\det M(U;\mu)]^* = \det M(U;-\mu)$.
- The signal-to-noise ratio for the phase $e^{i\varphi}$ decays exponentially: $\text{SNR}(e^{i\varphi}) \to \exp[-(F - F_{pq})/T]$, where $F - F_{pq}$ is the free energy difference between full and phase-quenched theories.
- The overlap problem is a manifestation of the sign problem: phase-quenched and full QCD occupy different regions of configuration space; the reweighting factor has exponentially bad SNR.
- Silver blaze: below the onset chemical potential $\mu_\text{onset} \approx m_N/3$, the baryon number density is zero despite $\mu \neq 0$, due to exact cancellations.
- Lefschetz thimbles deform integration to complex manifolds where $\text{Im}\, S = \text{const}$; global sign problem remains (relative phases between thimbles) plus residual sign problem on each thimble (Jacobian).
- Holomorphic gradient flow evolves $A$ according to $\frac{dA_\mu(x)}{dt} = -(\partial S/\partial A_\mu(x))^*$; improves average phase at small flow times.
- Complex Langevin dynamics extends the real stochastic process to the complexified field space; convergence is not guaranteed and correctness criteria must be checked.
- Contour deformations work for all theories (unlike some thimble methods); recent ML approaches find good contours efficiently.
- The sign problem is NP-hard in general [Troyer-Wiese 2005].

## Methods used

- QCD partition function: $Z = \int DU \det M(U;\mu) e^{-S_W[U]}$; $S_W$ is the Wilson action.
- Phase-quenched theory: replace $\det M(U;\mu)$ by $|\det M(U;\mu)|$; reweighting by phase gives full theory results.
- Taylor expansion of observables in $\mu_B/T$; valid for small $\mu_B/T \lesssim 3$.
- Imaginary chemical potential: $\mu = i\mu_I$ makes the determinant real; analytic continuation to real $\mu$.
- Lefschetz thimbles: stable thimbles $\mathcal{J}_k$ defined by $\partial_t \tilde{\phi}(x,t) = -\overline{(\partial S/\partial\tilde{\phi})}$ with $\tilde{\phi}|_{t=\infty} = \tilde{\phi}_0^{(k)}$ (critical point).
- Holomorphic gradient flow: same equation as thimble flow, defines integration contour $\mathcal{M}_t$.
- Complex Langevin: extends Langevin dynamics $dx/dt = -\partial S/\partial x + \eta$ to complex $x \in \mathbb{C}$.
- Benchmark models: Bose gas at non-zero $\mu$ (complex scalar $\phi^4$ in $d=4$); Thirring model at non-zero $\mu$.

## Notation introduced

- $\det M(U;\mu) = |\det M(U;\mu)|e^{i\varphi}$: fermion determinant decomposed into magnitude and phase.
- $\langle e^{i\varphi} \rangle_{pq}$: average phase in the phase-quenched theory (measures sign problem severity).
- $Z_{pq} = \exp[-(F_{pq}/T)]$, $Z = \exp[-(F/T)]$: phase-quenched and full partition functions; SNR $\sim e^{-(F-F_{pq})/T}$.
- $\mathcal{J}_k$: Lefschetz thimble attached to critical point $\tilde{\phi}_0^{(k)}$.
- $n_k$: intersection number of original contour with unstable thimble associated to $\mathcal{J}_k$.
- $S_\text{eff}[\phi] = S[\tilde{\phi}(\phi)] - \log J(\phi)$: effective action on the thimble including Jacobian.

## Open questions

- Whether the sign problem can be solved (not just alleviated) for QCD at large baryon density — likely NP-hard in general.
- Systematic combination of multiple methods (e.g. ML contour deformation + complex Langevin).
- Extension of ML-based contour deformation methods from scalar theories to full lattice QCD with gauge fields.
- Real-time dynamics (Minkowski path integrals) where sign problem is especially severe.
- Whether tensor network methods can reach large volumes needed for QCD thermodynamics.

## Links

- [[sign-problem]]
- [[overlap-problem]]
- [[silver-blaze-problem]]
- [[lefschetz-thimbles]]
- [[holomorphic-gradient-flow]]
- [[complex-langevin-dynamics]]
- [[contour-deformation]]
- [[reweighting]]
- [[average-phase]]
- [[partition-function-lattice]]
- [[gert-aarts]]
- [[denes-sexty]]
