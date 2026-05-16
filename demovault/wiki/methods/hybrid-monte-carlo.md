---
title: "Hybrid Monte Carlo (HMC)"
tags: [hmc, monte-carlo, molecular-dynamics, lattice-qft]
sources: ["[[temmen-2025-fully-ergodic-radial-updates]]", "[[aarts-2026-lattice-sign-problem-review]]"]
---

## Definition

Hybrid Monte Carlo (HMC) is an MCMC algorithm that generates field configurations $\phi$ distributed according to $p[\phi] \propto e^{-S[\phi]}$ by evolving an augmented system $(\phi, \pi)$ using molecular dynamics and then applying a Metropolis acceptance step to correct for discretisation errors.

## Intuition

The idea is to treat the action $S[\phi]$ as a potential energy and introduce fictitious conjugate momenta $\pi_i$ to get a total artificial Hamiltonian $H[\phi,\pi] = \frac{1}{2}\sum_i \pi_i^2 + S[\phi]$. Hamilton's equations evolve $(\phi,\pi)$ along trajectories that conserve $H$, exploring configuration space efficiently. Because the Hamiltonian is conserved exactly (in continuous time), any endpoint of a trajectory has the same weight $e^{-H}$ as the starting point, giving 100% acceptance. In practice, discretisation errors break exact conservation, and the Metropolis step corrects for this.

## Algorithm

1. Refresh momenta: draw $\pi_i \sim \mathcal{N}(0,1)$ independently.
2. Evolve $(\phi, \pi)$ for time $T$ using the leapfrog (or other symplectic) integrator with step size $\epsilon$ and $N_\text{MD} = T/\epsilon$ steps:
   - $\pi_{i+1/2} = \pi_i - \frac{\epsilon}{2}\frac{\partial S}{\partial\phi_i}$
   - $\phi_{i+1} = \phi_i + \epsilon\, \pi_{i+1/2}$
   - $\pi_{i+1} = \pi_{i+1/2} - \frac{\epsilon}{2}\frac{\partial S}{\partial\phi_{i+1}}$
3. Accept $(\tilde{\phi}, \tilde{\pi})$ with probability $\alpha_\text{HMC} = \min(1, e^{-\Delta H})$ where $\Delta H = H[\tilde{\phi},\tilde{\pi}] - H[\phi,\pi]$.
4. If rejected, keep $\phi' = \phi$.

Key parameter: trajectory length $T = \frac{\pi}{2}\sqrt{U\beta/N_t}$ (for Hubbard model) eliminates harmonic autocorrelations.

## Key results

- HMC satisfies detailed balance and is ergodic in principle; in practice it is susceptible to ergodicity violations from infinite potential barriers ([[ergodicity-violation]]).
- The leapfrog integrator is area-preserving and reversible, ensuring the Metropolis correction is valid.
- For the Hubbard model, manifolds of vanishing fermion determinant create infinite barriers that trap standard HMC — radial updates ([[radial-updates]]) are added to restore ergodicity.
- Autocorrelation times grow with system size; the leapfrog trajectory length is tuned to keep acceptance $\gtrsim 70\%$.

## Connections

- [[radial-updates]] — complementary update that restores ergodicity for HMC.
- [[ergodicity-violation]] — the failure mode of HMC in the presence of potential barriers.
- [[sign-problem]] — HMC cannot be applied directly when $e^{-S}$ is complex; requires phase quenching + reweighting.

## Open questions

- Optimal choice of integrator and trajectory length for different theories.
- Whether higher-order symplectic integrators significantly reduce autocorrelation for lattice QCD.
