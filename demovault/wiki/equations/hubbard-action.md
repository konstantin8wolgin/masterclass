---
title: "Hubbard model lattice action"
tags: [hubbard-model, lattice-action, fermion-determinant, exponential-discretisation]
sources: ["[[temmen-2025-fully-ergodic-radial-updates]]"]
---

## Statement

After applying the Hubbard-Stratonovich transformation and integrating out fermions, the effective Hubbard action for the auxiliary field $\phi_{tx}$ (one per site $x$ and time slice $t$) is:

$$S[\phi] = \frac{1}{2U\Delta_t}\sum_{t,x}\phi_{tx}^2 - \log\det M[\phi|\kappa]\det M[-\phi|-\kappa]$$

The fermion matrix $M[\phi|\kappa]$ in the exponential discretisation is:

$$M[\phi|\kappa]_{tx,t'y} = \delta_{t,t'}\delta_{x,y} - (e^{\kappa h})_{xy}\,e^{i\phi_{tx}}\mathcal{B}_t\,\delta_{t',t+1}$$

where $h$ is the hopping matrix, $\mathcal{B}_t = +1$ for $0 < t < N_t$ and $\mathcal{B}_0 = -1$ (anti-periodic boundary conditions).

## Variables

- $\phi_{tx}$: real auxiliary field (Hubbard-Stratonovich field) at time slice $t$, spatial site $x$
- $U$: on-site interaction strength (Hubbard $U$)
- $\kappa$: hopping parameter
- $\Delta_t = \beta/N_t$: temporal lattice spacing; $\beta = 1/T$ inverse temperature; $N_t$ number of time slices
- $h$: hopping matrix (nearest-neighbour on spatial lattice)
- $M[\phi|\kappa]$: fermion matrix for spin-up; $M[-\phi|-\kappa]$ for spin-down
- $\Phi_x = \sum_t \phi_{tx}$: summed auxiliary field per site (used for visualising ergodicity problem)

## Derivation sketch

The Hubbard Hamiltonian $H = H_K + H_U$ with $H_K = -\kappa\sum_{\langle x,y\rangle}(a^\dagger_x a_y - b^\dagger_x b_y)$ and $H_U = \frac{U}{2}\sum_x(a^\dagger_x a_x - b^\dagger_x b_x)^2$ is discretised in Euclidean time as $\beta/N_t$ slices. The Suzuki-Trotter decomposition followed by the continuous Hubbard-Stratonovich transformation $\exp(-\frac{\Delta_t U}{2}n^2) \propto \int d\phi\,\exp(-\frac{1}{2U\Delta_t}\phi^2 - i\phi n)$ decouples the quartic fermion interaction. Integrating out the fermions yields the fermionic determinant.

## Limits and special cases

- Strong coupling $U/\kappa \gg 1$: $\det M \propto \cos^{N_t}(\Phi_x/2)$ (approximately); infinite barriers at $\Phi_x = (2k+1)\pi$ create the ergodicity problem for HMC.
- Two-site model ($N_x = 2$, $N_t = 1$): configuration space is the plane $(\phi_1, \phi_2)$; barriers are diagonal bands at $\phi_1 = \pm\pi/2$, $\phi_2 = \pm\pi/2$.
- $U = 18$, $\kappa = 1$, $\beta = 1$: parameters used in demonstrations.

## Sources

- [[temmen-2025-fully-ergodic-radial-updates]] (Eq. 12–16; Sec. IV.A)
