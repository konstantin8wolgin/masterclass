---
title: "Average phase"
tags: [sign-problem, average-phase, diagnostic]
sources: ["[[lawrence-2023-deep-learning-fermion-sign]]", "[[aarts-2026-lattice-sign-problem-review]]"]
---

## Definition

The average phase $\langle\sigma\rangle$ (or $\langle e^{i\varphi}\rangle_{pq}$) is the expectation value of the complex phase of the Boltzmann factor, taken in the phase-quenched ensemble. It is the primary diagnostic for sign problem severity. When $\langle\sigma\rangle \approx 1$ the sign problem is mild; when $\langle\sigma\rangle \ll 1$ the sign problem is severe.

## Intuition

The full partition function relates to the phase-quenched one by $Z = Z_{pq} \langle e^{i\varphi}\rangle_{pq}$. If $\langle e^{i\varphi}\rangle_{pq}$ is small, computing $Z$ requires accurately measuring a quantity much smaller than its fluctuations — exponentially many samples are needed. Subtractions and contour deformations aim to increase $\langle\sigma\rangle$ by modifying the Boltzmann factor (or integration contour) to reduce phase fluctuations.

## Key results

- $\langle\sigma\rangle \sim e^{-(F-F_{pq})/T}$: the average phase decays exponentially with volume at fixed temperature and density.
- A subtraction $g = \nabla \cdot v$ changes the phase-quenched partition function $Z_{Q,\text{sub}}$ and hence the average phase, without changing physical expectation values ([[lawrence-2023-deep-learning-fermion-sign]]).
- The heavy-dense subtraction $g_\text{hd}$ analytically removes the dominant source of phase fluctuations at large $\mu$, yielding a substantial improvement in $\langle\sigma\rangle$ at $\mu/m_F \gtrsim 0.5$.
- The goal of subtraction methods is to make $\langle\sigma\rangle_\text{sub}$ as close to 1 as possible; a perfect subtraction would give $\langle\sigma\rangle_\text{sub} = 1$.

## Connections

- [[sign-problem]] — the average phase is the standard measure of sign problem severity.
- [[subtraction-method]] — directly modifies the average phase.
- [[contour-deformation]] — also modifies the average phase by changing the integration contour.

## Open questions

- Whether the average phase alone is the right diagnostic, or whether the SNR of specific observables is more informative.
