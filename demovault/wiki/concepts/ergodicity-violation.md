---
title: "Ergodicity violation"
tags: [ergodicity, hmc, monte-carlo, potential-barriers]
sources: ["[[temmen-2025-fully-ergodic-radial-updates]]", "[[aarts-2026-lattice-sign-problem-review]]"]
---

## Definition

An ergodicity violation occurs in an MCMC simulation when the Markov chain cannot reach all configurations of non-zero probability weight, because the trajectory is blocked by infinite (or very large) potential barriers. The simulation then fails to sample the correct distribution and observable measurements are biased.

## Intuition

HMC evolves configurations along molecular dynamics trajectories guided by the force $-\partial S/\partial\phi$. If $S$ diverges to $+\infty$ on a codimension-1 manifold in configuration space (e.g. where the fermion determinant vanishes), the force becomes infinite and the trajectory is repelled, never crossing the barrier. Regions of high probability on the other side of the barrier are never visited. The simulation is trapped in one sector of configuration space.

## Key results

- In the Hubbard model with exponential discretisation, the fermion matrix $M[\phi|\kappa]$ has manifolds of zeros that divide configuration space into disconnected sectors. Standard HMC cannot cross between sectors.
- Topological freezing in lattice QCD with chiral fermions is the same phenomenon: the topological charge cannot change because the barrier between topological sectors is infinite for chiral fermions.
- Radial updates ([[radial-updates]]) resolve this by proposing multiplicative jumps $\phi' = e^\gamma\phi$ that are agnostic to the diverging force — the proposal does not use the gradient.
- The acceptance probability $\alpha_R = \min(1, e^{-\Delta S + d\gamma})$ includes the Jacobian term $e^{d\gamma}$ (here $d$ is the field space dimension) and still satisfies detailed balance.
- Coarsening the MD integration (large step size) can also facilitate barrier crossings but at the cost of low acceptance rate and growing autocorrelation times — not a scalable solution.

## Connections

- [[hybrid-monte-carlo]] — the standard algorithm susceptible to ergodicity violations.
- [[radial-updates]] — the proposed fix via multiplicative MH proposals.
- [[sign-problem]] — related: topological freezing is a form of ergodicity violation in lattice QCD, also connected to the sign problem at non-zero density.

## Open questions

- General conditions under which a given theory will exhibit ergodicity violations with HMC.
- Whether parallel tempering or flow-based methods can substitute for radial updates in compact field spaces.
