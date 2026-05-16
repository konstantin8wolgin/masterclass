---
title: "Staggered Thirring model action"
tags: [thirring-model, lattice-action, fermion-sign-problem, staggered-fermions]
sources: ["[[lawrence-2023-deep-learning-fermion-sign]]"]
---

## Statement

The lattice action for the staggered Thirring model in $1+1$ dimensions with $U(1)$-valued link fields $A_\nu(x)$ is:

$$S = \sum_{x,\nu=0,1} \frac{2}{g^2}(1 - \cos A_\nu(x)) - \log\det K[A]$$

The Dirac matrix $K[A]$ is:

$$K[A]_{xy} = m\delta_{xy} + \frac{1}{2}\sum_{\nu=0,1}\left[\eta_\nu e^{iA_\nu(x)+\mu\delta_{\nu,0}}\delta_{x+\hat\nu,y} - \eta_\nu e^{-iA_\nu(y)-\mu\delta_{\nu,0}}\delta_{y+\hat\nu,x}\right]$$

The staggered fermion phases are $\eta_0 = (-1)^{\delta_{0,x_0}} = 1$ and $\eta_1 = (-1)^{x_0}$. Fermions live on a $\beta \times L$ lattice with bare mass $m$ and coupling $g^2$.

## Variables

- $A_\nu(x) \in [0, 2\pi)$: logarithm of $U(1)$ link field at site $x$, direction $\nu$
- $g^2$: bare coupling (strong coupling $g^2 = 1.0$ used in demonstrations)
- $m$: bare fermion mass ($m = 0.05$ in demonstrations)
- $\mu$: chemical potential (source of sign problem when non-zero)
- $\beta, L$: temporal and spatial extent of the lattice
- $\eta_\nu$: staggered fermion phase factor

## Derivation sketch

The Thirring model is a $1+1$D relativistic fermion with quartic interaction $(\bar\psi\gamma^\mu\psi)^2$. Bosonising the four-fermion interaction via a $U(1)$ gauge field $A_\mu$ and discretising with staggered fermions gives the lattice action above. The term $\frac{2}{g^2}(1-\cos A_\nu)$ is the pure gauge (bosonised) kinetic term; $\log\det K[A]$ is the fermion contribution after integrating out the fermions.

## Limits and special cases

- $0+1$ dimensional Thirring model (mechanics): $S_\text{mech} = \sum_t \frac{1}{2g^2}(1-\cos A(t)) - \log\det K[A]$; exactly solvable, partition function involves modified Bessel functions.
- $\mu \neq 0$: $\det K[A] \in \mathbb{C}$ since $K$ is no longer anti-Hermitian; sign problem arises.
- Heavy-dense limit $\beta\mu \gg 1$: $\det K \approx e^{\beta L\mu}(2^{-\beta L}e^{i\sum_x A_0(x)} + \mathcal{O}(e^{-\beta\mu}))$; used for analytic subtractions.
- Bare parameters for demonstrations: $m = 0.05$, $g^2 = 1.0$, lightest fermion mass $am_F = 0.35(2)$, lightest boson mass $am_B = 0.33(1)$.

## Sources

- [[lawrence-2023-deep-learning-fermion-sign]] (Eq. 3–6; Sec. II)
