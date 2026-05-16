---
title: "Complex Langevin dynamics"
tags: [complex-langevin, sign-problem, stochastic, holomorphic]
sources: ["[[aarts-2026-lattice-sign-problem-review]]", "[[ihssen-2026-solving-sign-problems-piks]]"]
---

## Definition

Complex Langevin (CL) dynamics extends the real Langevin stochastic process to a complex field $z = x + iy \in \mathbb{C}$ when the action $S(x)$ is complex. The drift term $-\partial S/\partial x$ becomes complex and drives the field off the real axis; the resulting stationary distribution is claimed to reproduce the target $e^{-S}$ expectation values under certain correctness criteria.

## Intuition

The Langevin equation $dx/dt = -\partial S/\partial x + \eta$ generates a Markov process on the real line whose stationary distribution is $p \propto e^{-S}$ when $S$ is real. When $S$ is complex, the drift $-\partial S/\partial z$ is complex, so the process explores $\mathbb{C}$ rather than $\mathbb{R}$. The idea is that the correct expectation values are recovered even though the trajectory is in $\mathbb{C}$, because holomorphic observables satisfy the same equations. However, convergence is not guaranteed, and correctness must be verified by separate criteria (e.g. fall-off of the probability distribution).

## Key results

- CL does not always converge to the correct distribution: runaway trajectories (excursions to large $|z|$) and convergence to wrong attractors are known failure modes.
- Correctness criteria: the probability distribution $P(x,y)$ must decay sufficiently fast for large $|z|$ [Aarts et al.].
- CL has been successfully applied to QCD at high temperature and imaginary $\mu$, and to the Bose gas.
- CL fails for the real-time harmonic oscillator in certain parameter regimes — the PIK-architecture works in these cases ([[ihssen-2026-solving-sign-problems-piks]]).
- CL is related to holomorphic gradient flow: both involve the complexified drift $\partial S/\partial z$; the CL drift is the reverse-sign version of the flow.

## Connections

- [[sign-problem]] — CL is the main alternative to contour deformation methods.
- [[lefschetz-thimbles]] — the holomorphic gradient flow defining thimbles uses the same complex drift as CL.
- [[physics-informed-kernels]] — PIKfolds work in regimes where CL fails; the two approaches are complementary.
- [[aarts-2026-lattice-sign-problem-review]] — comprehensive review including CL.

## Open questions

- Whether adaptive step-size methods or gauge cooling can completely eliminate convergence failures.
- Whether CL can be combined with PIK-architecture for theories where neither alone is sufficient.
- Correctness proof for CL applied to full lattice QCD at non-zero baryon density.
