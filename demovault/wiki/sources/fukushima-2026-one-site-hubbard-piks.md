---
title: "One-site Hubbard model with PIKs"
authors: [Kenji Fukushima, Friederike Ihssen, Renzo Kapust, Jan M. Pawlowski]
year: 2026
type: paper
tags: [physics-informed-kernels, sign-problem, hubbard-model, fermionic, qcd-proxy, lefschetz-thimbles, gaussian-mixture]
source_file: raw/papers/FermionPIK.pdf
---

## Summary

Full paper by Fukushima, Ihssen, Kapust, Pawlowski (in preparation) treating the one-site Hubbard model as a QCD proxy via the Dirac determinant structure $\prod_\alpha \det(1+Le^{\beta\mu})\det(1+L^\dagger e^{-\beta\mu})$. The central result is that the Boltzmann weight decomposes into three shifted Gaussians (the Gaussian composition of the effective action), each of which admits an exact analytic PIK solution with pure imaginary velocity $\dot\phi = iU$. The paper introduces PIKfold ripping, PIK equivalence classes, and the two-step PIK-pair resolution. Appendix A establishes the explicit connection between complex PIKs and Lefschetz thimbles, including the key identity $S_t(\phi_t(\phi_0)) - \mathrm{Tr}\log[\partial\phi_t/\partial\phi_0] = S_0(\phi_0)$.

---

## 1. QCD connection

The Dirac determinant in QCD at finite density reads:
$$\prod_{\text{phase space}} \det(1 + Le^{\beta\mu})\det(1 + L^\dagger e^{-\beta\mu})$$

The one-site Hubbard model realizes this structure with $L = e^{i\phi}$ (a $U(1)$ element), making it a tractable QCD proxy. Studying the sign problem here is directly relevant to the finite-density QCD sign problem. See [[sign-problem]].

---

## 2. Partition function and path integral

The one-site Hubbard Hamiltonian is $\hat{H} = U\hat{n}_\uparrow\hat{n}_\downarrow - \mu\hat{n}$. Its Fock-space spectrum $\{0, -\mu, -\mu, U-2\mu\}$ gives the exact partition function:
$$Z = \mathrm{Tr}[e^{-\beta H}] = 1 + 2e^{\beta\mu} + e^{\beta(2\mu - U)}. \tag{2}$$

After the [[hubbard-stratonovich-transformation]], the path integral representation is:
$$Z = \mathcal{N}\int d\phi\, \left(1 + e^{\beta(i\phi + \mu + U/2)}\right)^2 e^{-\beta\phi^2/2U} \tag{3}$$

with normalization $\mathcal{N} = \sqrt{\beta/(2\pi U)}$. The full action [[one-site-hubbard-fermionic-action]] (eq. 7):
$$S(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + e^{\beta(i\phi + \mu + U/2)}\right]$$

---

## 3. Gaussian decomposition of the partition function (eq. 6)

Expanding the squared bracket in the path integral gives:
$$Z = Z_0 + Z_1 + Z_2$$

with:
- $Z_0 = \mathcal{N}\int d\phi\, e^{-\beta\phi^2/2U}$ — free Gaussian, weight 1 (empty state)
- $Z_1 = 2\mathcal{N}\int d\phi\, e^{-\beta\phi^2/2U + i\beta\phi + \beta(\mu+U/2)}$ — shifted Gaussian, weight $2e^{\beta(\mu+U/2)}$ (single occupation, two spin states)
- $Z_2 = \mathcal{N}\int d\phi\, e^{-\beta\phi^2/2U + 2i\beta\phi + 2\beta(\mu+U/2)}$ — double-shifted Gaussian, weight $e^{2\beta(\mu+U/2)}$ (double occupation)

Each $Z_j$ is a Gaussian shifted by $j\cdot iU$ in $\phi$-space and integrable in closed form. See [[gaussian-mixture-decomposition]] and [[gaussian-composition-effective-action]].

---

## 4. Singular points / poles (eq. 8)

The action has logarithmic poles where $1 + e^{\beta(\mu + i\phi + U/2)} = 0$:
$$\phi_\star^{(n)} = \frac{(2n+1)\pi}{\beta} + i\!\left(\mu + \frac{U}{2}\right), \qquad n \in \mathbb{Z}$$

These are arranged on a horizontal line in $\mathbb{C}$ at $\mathrm{Im}(\phi_\star) = \mu + U/2$, periodically spaced by $2\pi/\beta$ along the real axis. Their position depends on $\mu$: at half-filling ($\mu = -U/2$) the poles sit on the real axis.

---

## 5. PIKfold ripping and path choices (Sec. IV)

When a pole $\phi_\star(t)$ crosses the sampling manifold $\mathcal{M}_t$ during RG-time evolution, the PIKfold fails — it rips. See [[pikfold-ripping]].

Two distinct choices of RG-time path $S_t$ give different pole dynamics:

**Path 1** (scale $\beta$ by $t$):
$$S_t(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + e^{t\beta(i\phi+\mu+U/2)}\right] \tag{9a}$$
Poles: $\mathrm{Re}(\phi_\star) = (2n+1)\pi/(t\beta)$, $\mathrm{Im}(\phi_\star) = \mu + U/2$. Poles start at $|\mathrm{Re}(\phi_\star)| = \infty$ (safely away from any manifold near the real axis) and move inward as $t$ increases. Ripping occurs if the moving pole reaches the manifold.

**Path 2** (scale $U/2$ by $t$, eq. 10a):
$$S_t(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + e^{\beta(i\phi+\mu+tU/2)}\right] \tag{10a}$$
Poles: $\mathrm{Re}(\phi_\star) = (2n+1)\pi/\beta$, $\mathrm{Im}(\phi_\star) = \mu + tU/2$. Poles move in the imaginary direction as $t$ increases.

For $\mu \lesssim \mu + U/2$, both paths keep the manifold away from poles; for $\mu \approx 0.1$ (near half-filling), Path 1 causes ripping.

---

## 6. Effective action decomposition for PIK (Sec. V, eq. 15)

The full $t$-dependent weight splits as [[gaussian-composition-effective-action]]:
$$e^{-S_t(\phi)} = e^{-\frac{\beta}{2U}\phi^2} + 2e^{t\beta(\mu + U/2 + i\phi) - \frac{\beta}{2U}\phi^2} + e^{2t\beta(\mu + U/2 + i\phi) - \frac{\beta}{2U}\phi^2}$$

At $t=0$: pure free Gaussian. At $t=1$: exact target weight $e^{-S(\phi)}$. Each term is a Gaussian shifted by $j\cdot itU$ in $\phi$-space.

---

## 7. Analytic PIK solution (eqs. 16–18)

For each shifted Gaussian component $j$, the Wegner equation ([[wegner-equation]]) admits the exact solution:
$$\dot\phi_t = iU \qquad \Longrightarrow \qquad \phi_t(\phi_0) = \phi_0 + iUt$$

The field moves with constant imaginary velocity $iU$ — a **pure imaginary shift**. This is an exact analytic kernel, no numerical approximation required. Compare to the free-theory complex shift in [[ihssen-2026-solving-sign-problems-piks]] (eq. 62: $\dot\phi = -i/m^2$).

Derivation: for the $j$-th component action $S_t^{(j)} = \frac{\beta}{2U}\phi^2 - jt\beta(\mu+U/2+i\phi) + \text{const}$, the Wegner equation with ansatz $\dot\phi = c$ (constant) forces $c = iU$.

---

## 8. Observable estimator via weighted sum (eq. 20)

For samples $\{\phi_i\}$ drawn from the free Gaussian $p_0 \propto e^{-\beta\phi^2/2U}$:
$$\langle\mathcal{O}\rangle = \frac{\sum_i \left[\mathcal{O}(\phi_i) + 2\mathcal{N}_1\,\mathcal{O}(\phi_i + iU) + \mathcal{N}_2\,\mathcal{O}(\phi_i + 2iU)\right]}{1 + 2\mathcal{N}_1 + \mathcal{N}_2}$$

with occupation weights $\mathcal{N}_1 = e^{\beta(\mu+U/2)}$, $\mathcal{N}_2 = e^{2\beta(\mu+U/2)}$. The three field arguments are the images under the three analytic PIK flows at transport time $t=1$.

---

## 9. Two-step PIK-pair resolution (Sec. VI)

When the direct PIK path causes ripping, the two-step procedure:

**Step 1**: Find $p_t$ that prevents ripping, using a log-sum-exp parametrization:
$$\mathcal{T}_t = \log\sum_i c_{i,t}\,\rho_{i,t}(\phi), \qquad c_{0,t} = 1, \quad \rho_{0,t} = \bar\rho_0 \tag{25a}$$
with $t$-dependent centers $\phi^{(i)}_{m,t}$ and covariances $\mathcal{C}^{(i)}_{nm,t}$ (Gaussian components, eq. 25b).

**Step 2**: Efficiently sample from $p_t$ using the analytic Gaussian structure.

---

## 10. PIK equivalence classes (eq. 22)

Two distributions $p \simeq q$ iff there exists a PIKfold $\mathcal{P} = \{\mathcal{M}_t\}$ connecting them:
$$p \simeq q \;:\; \exists\,(p_t, \mathcal{P}) \;\text{ with }\; p_0 = p,\; p_1 = q$$

This is an equivalence relation on the space of distributions. Distributions in the same class can be sampled using the same PIK architecture; the fermionic and bosonic terms of the Hubbard model lie in different equivalence classes.

---

## 11. PIKs vs. Lefschetz thimbles (Appendix A)

Key differences between LT and complex PIK:
- LT requires summing over many thimbles with inter-thimble phases; complex PIK lives on a **single** manifold.
- PIK is automatically weight-preserving ($\theta = 0$ everywhere on $\mathcal{M}_t$); LT leaves a residual phase from the Jacobian of the real parametrisation.
- Constructing LT via PIK: define auxiliary action $\tilde{S}_t$ targeting the thimble $\mathcal{J}_\sigma$; kernel equation (A15).

Key identity (eqs. A13/B13) — action plus log-Jacobian is conserved along the flow:
$$S_t\!\left(\phi_t(\phi_0)\right) - \mathrm{Tr}\log\!\left[\frac{\partial\phi_t}{\partial\phi_0}\right] = S_0(\phi_0)$$

At $t=1$, this gives the thimble sampling weight: $e^{-S_0(\phi_0)}$ exactly when the Jacobian is accounted for.

See [[lefschetz-thimbles]], [[contour-deformation]].

---

## 12. General complex sampling framework (Appendix B)

Proof of the weight-preservation identity for complex distributions in the sampling-pair formalism. The PIK architecture maps $(p_0(\phi), \mathcal{M}_0) \to (p_1(\phi), \mathcal{M}_1)$ with $p_1 = p$ and preserves all expectation values. Generalizes [[ihssen-2026-solving-sign-problems-piks]] Section 2 to the fermionic context.

---

## Notation summary

| Symbol | Meaning |
|--------|---------|
| $U$ | on-site Coulomb interaction |
| $\mu$ | chemical potential |
| $\beta$ | inverse temperature |
| $\phi$ | Hubbard-Stratonovich auxiliary field |
| $Z_j$ | $j$-th Gaussian component of partition function ($j=0,1,2$) |
| $\phi_\star^{(n)}$ | poles of $e^{-S}$ in complexified field space |
| $\mathcal{N}_j = e^{j\beta(\mu+U/2)}$ | occupation weight for $j$-fold occupied state |
| $\mathcal{T}_t$ | log-sum-exp PIRG action for two-step resolution |

---

## Open questions

- Extension to multi-site Hubbard ($N_x > 1$): does $\dot\phi = iU$ persist or acquire corrections from the fermion determinant?
- Optimal path choice (Path 1 vs Path 2) for minimizing ripping risk across the $(\beta, \mu, U)$ parameter space.
- LT construction via PIK in higher-dimensional theories: tractability and computational cost.
- Whether PIK equivalence classes provide a useful topological classification of sign problems.

---

## Links

- [[physics-informed-kernels]]
- [[physics-informed-rg-flows]]
- [[sign-problem]]
- [[wegner-equation]]
- [[lefschetz-thimbles]]
- [[contour-deformation]]
- [[one-site-hubbard-fermionic-action]]
- [[gaussian-mixture-decomposition]]
- [[gaussian-composition-effective-action]]
- [[hubbard-stratonovich-transformation]]
- [[pikfold-ripping]]
- [[hubbard-action]]
- [[ihssen-2026-solving-sign-problems-piks]]
- [[ihssen-2026-fermionic-sign-piks]]
