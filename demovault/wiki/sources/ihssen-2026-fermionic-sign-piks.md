---
title: "Fermionic Sign Problem with Physics-Informed Kernels"
authors: [Friederike Ihssen, Renzo Kapust, Jan M. Pawlowski]
year: 2026
type: paper
tags: [physics-informed-kernels, sign-problem, hubbard-model, fermionic, gaussian-mixture, one-site, cdf-method]
source_file: raw/papers/FermionikPIK.pdf
---

## Summary

Working note extending the complex PIK architecture ([[ihssen-2026-solving-sign-problems-piks]]) to fermionic sign problems, using the one-site Hubbard model as the primary testbed. The paper establishes two complementary approaches: (1) an analytic CDF-based reformulation of the Wegner equation, yielding the PIK velocity as $\dot\phi_t = -(1/p_t)(d/dt)F_t$ where $F_t$ is the CDF of $p_t$; and (2) a PIKfold Mixture approach that applies a separate PIK per Gaussian component. The sign problem severity is analyzed via $Q_{sp} = Z/Z_{pq}$ for three contour choices. A KL-divergence optimization fitting a Gaussian mixture to $p$ is proposed as a practical criterion for the PIKfold Mixture.

---

## 1. One-site Hubbard model

**Hamiltonian**: $H = U\hat{n}_\uparrow\hat{n}_\downarrow - \mu\hat{n}$ where $\hat{n}_\sigma$ is the spin-$\sigma$ number operator.

**Partition sum** (eq. 2):
$$\mathcal{N} = \mathrm{Tr}\!\left[e^{-\beta H}\right] = 1 + 2e^{\beta\mu} + e^{\beta(2\mu - U)}$$

This follows directly from the four Fock-space eigenvalues. The three terms correspond to zero, single (two spin states), and double occupation.

**Path integral** (eq. 3): after the [[hubbard-stratonovich-transformation]]:
$$\mathcal{N} = \sqrt{\frac{\beta}{2\pi U}}\int_{-\infty}^{\infty} d\phi\, e^{-S(\phi)}, \qquad S(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + e^{\beta(\mu + i\phi + U/2)}\right]$$

This is [[one-site-hubbard-fermionic-action]]. The imaginary coupling $i\phi$ makes $S(\phi) \in \mathbb{C}$ for $\phi \in \mathbb{R}$, causing the sign problem.

---

## 2. Partition sum as Gaussian mixture (eq. 4)

The distribution $p(\phi) = e^{-S(\phi)}$ is a mixture of three Gaussians:
$$\mathcal{N} = \sqrt{\frac{\beta}{2\pi U}}\int d\phi\left[e^{-\frac{\beta}{2U}\phi^2} + 2e^{\beta(\mu+U/2+i\phi)}e^{-\frac{\beta}{2U}\phi^2} + e^{2\beta(\mu+U/2+i\phi)}e^{-\frac{\beta}{2U}\phi^2}\right] \tag{4}$$

Each term is a Gaussian shifted by a purely imaginary amount. The key identity (eq. 7) that makes each term integrable:
$$\int_{-\infty}^{+\infty} dx\, e^{-ax^2 + bx + c} = \sqrt{\frac{\pi}{a}}\, e^{\frac{b^2}{4a} + c}$$

The physics of the mixture: part of the interesting structure comes from the distribution being a mixture of distributions rather than a single one. The logarithm in $S(\phi)$ (which will later host a determinant in multi-site models) has poles in the complexified space. See [[gaussian-mixture-decomposition]] and [[gaussian-composition-effective-action]].

---

## 3. Sign problem severity (Sec. 1, eq. 6)

$$Q_{sp} = \frac{\int_{-\infty}^{+\infty} d\phi\, I(\phi)}{\int_{-\infty}^{+\infty} d\phi\, |I(\phi)|}$$

where $I(\phi) = e^{-S(\phi)}$ is the full integrand. $Q_{sp} = 1$: no sign problem; $Q_{sp} \to 0$: intractable.

Three contour choices analyzed (Fig. 1):
- **Original**: $\phi \in \mathbb{R}$. $Q_{sp}$ drops from 1 to 0 as $\mu$ increases past 0.
- **Particle-hole**: $\phi \to -\phi + 2iU$. Improves $Q_{sp}$ for large $\mu$.
- **Half-filling**: $\phi \to \phi + iU$. $Q_{sp} = 1$ at $\mu = U/2$; best contour at half-filling.

The contour choices correspond to different starting points for a PIK flow.

---

## 4. Poles in complexified field space (eqs. 8–9)

The action has logarithmic singularities where $1 + e^{\beta(\mu + i\phi_\star + U/2)} = 0$:
$$\mathrm{Re}(\phi_\star) = \frac{\pi + 2\pi k}{\beta}, \quad k \in \mathbb{Z} \qquad \text{and} \qquad \mathrm{Im}(\phi_\star) = \mu + \frac{U}{2}$$

Poles sit on a horizontal line in $\mathbb{C}$ at height $\mathrm{Im}(\phi_\star) = \mu + U/2$, spaced by $2\pi/\beta$. The sign problem onset (Fig. 2) is linked to the pole position crossing $\mathrm{Im}(\phi) = 0$ (half-filling) and then $\mathrm{Im}(\phi) = 2U$ (where a new connected real contour emerges).

---

## 5. PIKfold parametrization (eqs. 16–17)

A one-Pikfold approach uses the action path:
$$S_t(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + e^{t\beta(\mu + i\phi + U/2)}\right] \tag{16}$$

Under this parametrization the poles move as:
$$\mathrm{Re}(\phi_\star) = -\frac{\pi + 2\pi k}{t\beta}, \qquad \mathrm{Im}(\phi_\star) = \mu + \frac{U}{2} \tag{17}$$

Poles start at $|\mathrm{Re}(\phi_\star)| = \infty$ (safe from the manifold at $t=0$) and move inward as $t \to 1$. This controls [[pikfold-ripping]].

---

## 6. Analytical solution via CDF reformulation (Sec. 2.1, eqs. 10–14)

The continuity equation for the flowing PDF:
$$\frac{dp_t(\phi)}{dt} = -\frac{\partial}{\partial\phi}\!\left(\dot\phi_t(\phi)\,p_t(\phi)\right) \tag{10}$$

with boundary condition $\lim_{\phi\to\pm\infty} p_t \dot\phi_t = 0$ (eq. 11). Integrating from $-\infty$ to $\phi$:

$$\frac{d}{dt}\int_{-\infty}^\phi p_t(\psi)\,d\psi = -\dot\phi_t(\phi)\,p_t(\phi) \tag{12}$$

Since $F_t(\phi) = \int_{-\infty}^\phi p_t(\psi)\,d\psi$ (the CDF, eq. 13), this gives the **closed-form PIK velocity** (eq. 14):

$$\boxed{\dot\phi_t(\phi) = -\frac{1}{p_t(\phi)}\frac{d}{dt}F_t(\phi)}$$

For the Gaussian mixture $p_t = \sum_i w_i p_{t,i}$:
$$F_t(\phi) = \sum_i w_i F_{t,i}(\phi) \tag{15}$$

Each $F_{t,i}$ is the CDF of a Gaussian = an error function, analytically continued to $\mathbb{C}$. The analytic solution is therefore given as sums and ratios of Gaussians and error functions.

---

## 7. One-PIKfold solution: basis expansion (eq. 18–19)

To solve the Wegner equation ([[wegner-equation]]) numerically, expand the kernel:
$$\dot\phi_t(\phi) = \sum_{i=0}^N k_{i,t}\,K_{i,t}(\phi) \tag{18}$$

Basis functions (eq. 19):
$$K_{i,t}(\phi) = \begin{cases} \phi^i/i! & i < N_\mathrm{poly} \\ 1/(\phi - \phi_\star^{(n)}(t)) & i \geq N_\mathrm{poly} \end{cases}$$

The polynomial part captures smooth structure; the rational (pole basis) part captures the singularity structure of the complexified action near $\phi_\star^{(n)}(t)$. Hermite polynomials were also tested but performed worse than the monomial basis.

The PIKfold grows in the imaginary direction for larger chemical potentials (Figs. 3, 5). A single PIKfold can soften but not fully cure the sign problem in the critical $\mu$ regimes (consistent with general results for the Thirring model, [[lawrence-2023-deep-learning-fermion-sign]]).

---

## 8. PIKfold Mixture (Sec. 2.1.3)

Since a single PIKfold cannot cure the sign problem for all $\mu$, write $p$ as a Gaussian mixture:
$$p(\phi) = \sum_k w_k \hat\rho_k(\phi), \qquad \hat\rho_k(\phi) = e^{-\frac{1}{2}(\phi - \psi_k)^T \mathcal{C}_k^{-1}(\phi - \psi_k)} \tag{22}$$

with complex covariance $\mathcal{C}_k$ and mean $\psi_k$. Apply a separate PIK to each component $\hat\rho_k$; each component PIK maps one Gaussian to a sign-problem-free contour via:
$$(\mathcal{C}_{k,xy}, \psi_{k,x} \in \mathbb{C}) \;\to\; (\hat G_{k,xy}, m_{k,x} \in \mathbb{R})$$

The weights $w_k$ can be made real throughout (by shifts of the mean and covariance parameters).

---

## 9. KL divergence optimization (Sec. 2.1.3, eqs. 24–29)

Fit a modelled distribution $\underline{p}$ to the true $p$ via KL divergence. For the one-site Hubbard model, the modelled distribution is parametrized as (eq. 28):
$$\underline{p}(\phi) = \frac{1}{\underline{\mathcal{N}}}\sum_k e^{\omega_k}\, e^{-\frac{1}{2\sigma^2}(\phi - im_k)^2}$$

for real parameters $\sigma, m_k, \omega_k$ (purely imaginary centers $im_k$ reflect the Gaussian mixture structure of $p$). The natural variance choice is $\sigma^2 = U/\beta$.

The complex KL divergence (eq. 26):
$$``D_{KL}(\underline{p}\,\|\,p)" = \mathbb{E}_{k\sim\mathrm{Cat}(w),\,\varphi\sim\underline\rho_k}\!\left[\frac{1}{2}\!\left(\frac{1}{\sigma^2} - \frac{1}{\underline{\sigma}^2}\right)\phi_k^2(\varphi) - \log\frac{\det M(\phi_k(\varphi))}{\det\underline{M}(\phi_k(\varphi))}\right] + \log\mathcal{N} - \log\underline{\mathcal{N}}$$

Since the KL divergence is complex, the practical optimization criterion (eq. 29):
$$``D_{KL}^{|\cdot|}(\underline{p}\,\|\,p)" = \mathbb{E}\!\left[|S(\phi_k(\varphi)) - S(\phi_k(\varphi))|^2\right]$$

has a global zero at the optimal parameters and can be efficiently computed without a sign problem. KL divergence history during optimization is shown in Fig. 6; convergence to machine precision is achieved for the one-site model at $U=1$, $\beta=10$, $\mu=1$.

---

## Notation summary

| Symbol | Meaning |
|--------|---------|
| $\mathcal{N}$ | partition function $= 1+2e^{\beta\mu}+e^{\beta(2\mu-U)}$ |
| $Q_{sp}$ | sign problem severity $= Z/Z_{pq}$ |
| $F_t(\phi)$ | CDF of $p_t$ at flow time $t$ |
| $\phi_\star^{(k)}$ | poles at $\mathrm{Re}=(2k+1)\pi/t\beta$, $\mathrm{Im}=\mu+U/2$ |
| $K_{i,t}(\phi)$ | kernel basis: monomials + pole terms $1/(\phi-\phi_\star^{(n)})$ |
| $\underline{p}(\phi)$ | modelled Gaussian mixture distribution |
| $m_k, \omega_k$ | imaginary centers and log-weights of modelled Gaussians |

---

## Open questions

- Convergence of the CDF-based analytic PIK for theories with multiple poles per period strip.
- Extension to finite lattice (multi-site Hubbard): does the mixture structure persist in a useful form?
- Systematic comparison of the single-PIKfold (basis expansion) vs PIKfold Mixture approaches in terms of sign problem suppression.
- Whether the $|D_{KL}|^2$ criterion selects the globally optimal Gaussian mixture or only a local optimum.

---

## Links

- [[physics-informed-kernels]]
- [[physics-informed-rg-flows]]
- [[sign-problem]]
- [[wegner-equation]]
- [[one-site-hubbard-fermionic-action]]
- [[gaussian-mixture-decomposition]]
- [[gaussian-composition-effective-action]]
- [[hubbard-stratonovich-transformation]]
- [[pikfold-ripping]]
- [[hubbard-action]]
- [[average-phase]]
- [[ihssen-2026-solving-sign-problems-piks]]
- [[fukushima-2026-one-site-hubbard-piks]]
