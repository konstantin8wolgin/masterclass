---
title: "Radial update acceptance probability"
tags: [radial-updates, metropolis-hastings, detailed-balance]
sources: ["[[temmen-2025-fully-ergodic-radial-updates]]"]
---

## Statement

The acceptance probability for a radial update proposal $\phi' = e^\gamma\phi$ with $\gamma \sim \mathcal{N}(0, \sigma_R^2)$ is:

$$\alpha_R = \min\!\left(1,\, e^{-\Delta S + d\gamma}\right)$$

where $\Delta S = S[\phi'] - S[\phi]$ is the change in action and $d$ is the dimension of the field space (number of real field components).

## Variables

- $\phi = (\phi_1, \ldots, \phi_d)$: current field configuration
- $\phi' = e^\gamma\phi = (e^\gamma\phi_1, \ldots, e^\gamma\phi_d)$: proposed configuration
- $\gamma$: log-scale factor drawn from $\mathcal{N}(0, \sigma_R^2)$
- $\sigma_R$: proposal standard deviation (key tuning parameter; $\sigma_R \propto 1/\sqrt{d}$ for optimal efficiency)
- $\Delta S = S[\phi'] - S[\phi]$: action change
- $d$: field space dimension (provides the Jacobian factor $e^{d\gamma}$)

## Derivation sketch

The Jacobian of the multiplicative map $\phi \to e^\gamma\phi$ on $\mathbb{R}^d$ is $\det(\partial\phi'/\partial\phi) = e^{d\gamma}$ (each component scaled by $e^\gamma$, Jacobian is $e^{\gamma \cdot d}$). For detailed balance, the ratio of forward to backward proposal probabilities is:

$$\frac{T(\phi\to\phi')}{T(\phi'\to\phi)} = \frac{q(\gamma)}{q(-\gamma)} \cdot \frac{|J(\phi\to\phi')|}{|J(\phi'\to\phi)|}$$

Since $q(\gamma) = q(-\gamma)$ (symmetric Gaussian) and $|J(\phi\to\phi')| = e^{d\gamma}$, $|J(\phi'\to\phi)| = e^{-d\gamma}$, the ratio is $e^{2d\gamma}$... wait, more carefully: the acceptance condition from detailed balance gives $\alpha_R = \min(1, e^{-\Delta S + d\gamma})$, where $d\gamma$ is the log-Jacobian of the forward proposal.

## Limits and special cases

- $\gamma = 0$: $\phi' = \phi$, always accepted ($\alpha_R = 1$). No update.
- Large $\sigma_R$: large jumps proposed; many rejections. Extreme $\sigma_R$ is inefficient.
- Small $\sigma_R$: short-range moves; efficient crossing of barriers is slow.
- Optimal $\sigma_R$: target acceptance rate $\sim 50\%$–$70\%$ by analogy with other MH proposals; scales as $1/\sqrt{d}$.
- If $S[\phi] \propto \phi^2$ (pure Gaussian): $\Delta S = (e^{2\gamma}-1)S[\phi]$; acceptance depends on how far from the Gaussian peak the current configuration is.

## Sources

- [[temmen-2025-fully-ergodic-radial-updates]] (Eq. 10: acceptance formula; Sec. III.A: derivation and detailed balance proof)
