---
title: "Radial updates"
tags: [radial-updates, metropolis-hastings, ergodicity, hmc]
sources: ["[[temmen-2025-fully-ergodic-radial-updates]]"]
---

## Definition

A radial update is a multiplicative Metropolis-Hastings (MH) proposal that rescales the entire field configuration $\phi = (\phi_1, \ldots, \phi_d)$ by a log-normally distributed factor $e^\gamma$, where $\gamma \sim \mathcal{N}(0, \sigma_R^2)$. The resulting proposal $\phi' = e^\gamma\phi$ scales the radial coordinate $R = \sqrt{\sum_i \phi_i^2}$ by $e^\gamma$ while preserving all angular information.

## Intuition

HMC crosses barriers by following molecular dynamics trajectories guided by the force $-\partial S/\partial\phi$. Where the force diverges (infinite potential barrier), the trajectory is repelled and cannot cross. A multiplicative proposal is agnostic to the force: it simply jumps the configuration to a new radius in field space. Because the proposal is symmetric on a log scale ($\gamma$ and $-\gamma$ are equally probable), the acceptance probability has a clean form with only the action change and the Jacobian.

## Algorithm

Per HMC step, perform one (or more) radial updates:

1. Draw $\gamma \sim \mathcal{N}(0, \sigma_R^2)$.
2. Propose $\phi' = e^\gamma \phi$ (rescale all field components by the same factor $e^\gamma$).
3. Compute the change in action $\Delta S = S[\phi'] - S[\phi]$.
4. Accept with probability $\alpha_R = \min(1, e^{-\Delta S + d\gamma})$ where $d$ is the dimensionality of the field (Jacobian of the multiplicative map).
5. If rejected, keep $\phi' = \phi$.

**Key parameter:** $\sigma_R$ — the proposal standard deviation. Larger $\sigma_R$ allows larger jumps but lower acceptance. Optimal $\sigma_R \propto 1/\sqrt{d}$.

## Key results

- The combined HMC + radial update algorithm satisfies detailed balance: $p[\phi]\Omega(\phi\to\phi') = p[\phi']\Omega(\phi'\to\phi)$.
- For the toy model $p[x] \propto \prod_i \cos^2(x_i)e^{-\beta x_i^2}$, standalone HMC is trapped in the central sector; radial updates with $\sigma_R = 1.75$ resolve the ergodicity violation completely.
- For the two-site Hubbard model: HMC is confined to the middle diagonal band; radial updates with $\sigma_R = 1.8$ explore the full configuration space.
- Autocorrelation times are reduced simultaneously with the restoration of ergodicity.
- Radial updates are specific to non-compact fields; compact fields (e.g. $U(1)$ link variables) require a different approach.

## Connections

- [[hybrid-monte-carlo]] — radial updates are interleaved with HMC steps.
- [[ergodicity-violation]] — the problem that radial updates solve.
- [[radial-update-acceptance]] — the acceptance probability formula.
- [[hubbard-action]] — primary application in this context.

## Open questions

- Optimal ratio of radial updates to HMC steps as a function of system size.
- Extension to compact field spaces.
- Theoretical bound on the improvement in autocorrelation time.
