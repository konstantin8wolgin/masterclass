---
title: "Addressing the sign-problem in Euclidean path integrals with radial basis function neural networks"
authors: [Gábor Balassa]
year: 2026
type: paper
tags: [sign-problem, rbf-networks, path-integral, finite-density, silver-blaze, complex-scalar-field]
source_file: raw/papers/2510.01695v3.pdf
---

## Summary

Expands the non-linear interaction term $e^{-S_{1,x}}$ of the Euclidean path integral as a sum of Gaussian (radial basis function) kernels, then performs a momentum-space diagonalisation that factorises the partition function. This converts an $\mathcal{O}(K^{N^4})$ problem into $\mathcal{O}(KN^4)$. Applied to complex scalar $\phi^4$ at finite chemical potential in $3+1$D; finds $\mu_c = 1.17 \pm 0.018$.

**When to use**: when the action splits into a tractable quadratic part (diagonal in momentum space) and a nonlinear self-interaction part that can be well-approximated by Gaussians.

---

## 1. Continuum and lattice action

Continuum complex scalar $\phi^4$ in $3+1$D:
$$S = \int d^4x\left[|\partial_\mu\phi|^2 + (m^2-\mu^2)|\phi|^2 + \mu(\phi^*\partial_4\phi - \phi\partial_4\phi^*) + \lambda(|\phi|^2)^2\right]. \tag{1}$$

Write $\phi = \frac{1}{\sqrt{2}}(\phi_1 + i\phi_2)$ as two real fields. On an $N^4$ lattice with unit spacing, the discretised action splits as:
$$S = \sum_x\left[S_{0,x}(\mu) + S_{1,x}(m,\lambda)\right], \tag{2}$$
where:

**$\mu$-dependent / kinetic part** (source of the sign problem):
$$S_{0,x}(\mu) = 4\sum_{a=1,2}\phi_{a,x}^2 - \sum_{i=1}^3\sum_{a=1,2}\phi_{a,x}\phi_{a,x+\hat{i}} - \cosh(\mu)\sum_{a=1,2}\phi_{a,x}\phi_{a,x+\hat{4}} - i\sinh(\mu)\left[\phi_{1,x}\phi_{2,x+\hat{4}} - \phi_{2,x}\phi_{1,x+\hat{4}}\right]. \tag{3}$$

**Self-interaction part** (to be approximated by RBF):
$$S_{1,x}(m,\lambda) = \frac{m^2}{2}\sum_{a=1,2}\phi_{a,x}^2 + \frac{\lambda}{4}\!\left(\sum_{a=1,2}\phi_{a,x}^2\right)^{\!2}. \tag{4}$$

The imaginary term $i\sinh(\mu)[\ldots]$ in $S_{0,x}$ makes the full action complex and causes the sign problem.

---

## 2. RBF expansion

Approximate the nonlinear interaction factor at each lattice site as:
$$e^{-S_{1,x}(m,\lambda)} \approx \sum_{k=1}^K a_k\,e^{-A\sum_{a=1,2}(\phi_{a,x} - c_{a,k})^2}. \tag{5}$$

Parameters:
- $K$: number of Gaussian kernels
- $A > 0$: common width parameter (kept constant across kernels)
- $c_k = (c_{1,k}, c_{2,k})$: 2D centres in the $(\phi_1, \phi_2)$ plane
- $a_k$: weight coefficients (to be determined by least squares)

**Fitting the weights**: generate $N_T$ training points $(\phi_{1,x}, \phi_{2,x})$ uniformly in the range where $e^{-S_{1,x}}$ is non-negligible. Build the coefficient matrix $H(c_k, A) \in \mathbb{R}^{N_T \times K}$ and solve:
$$\hat{a} = \arg\min_a\|H(c_k, A)a - y\|^2, \tag{6}$$
where $y_i = e^{-S_{1,x}(m,\lambda)}$ evaluated at training points. Use SVD with regularisation $r_c = 0.001$ (truncate very small singular values).

**Practical parameters** (paper benchmark, $m=1$, $\lambda=2$): $K=10^2$ centres on $[-1,1]\times[-1,1]$ grid, $\Delta c = 2/(10-1)$, $A=3$. Training range: $\phi_{1,x},\phi_{2,x} \in [-2.5, 2.5]$, $N_T = 1000$.

**Model selection**: scan $K \in [5^2, 30^2]$ and $A \in [1, 15]$ (for silver blaze: $A \in [0.5, 3.1]$, $C \in [0.5, 4.1]$), choose configuration minimising RMSE on a held-out test set. Optimal for silver blaze calculation: $K=10^2$, $A \approx 1.7$, $C \approx 2.9$ (centres in $[-C,C]\times[-C,C]$).

---

## 3. Partition function after RBF

After substituting the RBF expansion, the partition function becomes:
$$Z = \int\prod_x d\phi_{1,x}\,d\phi_{2,x}\,e^{-S_{0,x}(\mu)}\cdot F[\phi_{a,x}], \tag{7}$$
where
$$F[\phi_{a,x}] = \sum_{k=1}^K a_k\,e^{-A\sum_a\phi_{a,x}^2 + 2A\sum_a c_{a,k}\phi_{a,x} - A\sum_a c_{a,k}^2}. \tag{8}$$

The sum $S_{0,x}$ can be written in matrix form as $S_{0,x} \to \frac12\phi^T M(\mu)\phi$ for a general $\mu$-dependent quadratic matrix $M(\mu)$. Then:
$$Z \approx \sum_{k\in\{1,\ldots,K\}^{N^4}}\hat{a}_k\int\mathcal{D}\phi\,e^{-\frac12\phi^T W\phi + 2A(\tilde{c}_k^T\phi) - A(1^T\tilde{c}_k^2)}, \tag{9}$$
where $W = [M(\mu) - 2\bar{A}]$, $\bar{A}$ is a diagonal matrix with $A$ on the diagonal, and $\hat{a}_k = \prod_x a_{k_x}$ is the product of all weight factors for a specific kernel combination.

---

## 4. Momentum-space diagonalisation

Apply a unitary similarity transformation $U$ that diagonalises $W = M(\mu) - 2\bar{A}$ in the 2D $(\tilde\Phi_{1,p}, \tilde\Phi_{2,p})$ momentum space, generating a block-diagonal matrix $\tilde{W} = U^TWU$.

Transform fields: $\tilde\Phi = U^T\phi$, where $\tilde\Phi_p$ are the momentum-space fields at momenta $p_i = 2\pi n_i/N$, $-N/2 < n_i \leq N/2$.

After transformation the partition function factorises over momenta:
$$Z \to \prod_p\sum_{k=1}^K Q_k[\tilde\Phi_{a,p}, \mu, A, c_{a,k}, a_k], \tag{11}$$
where each factor $Q_k$ is a Gaussian integral over the $2\times2$ block at momentum $p$.

**Closed-form result** (eq 12):
$$Z \approx \sum_{k\in\{1,\ldots,K\}^{N^4}}\hat{a}_k\prod_p\frac{2\pi}{\sqrt{\mathcal{A}_p^2+\mathcal{B}_p^2}}\exp\!\left[-A(c_{1,k_p}^2+c_{2,k_p}^2) + \frac{2A^2}{\mathcal{A}_p^2+\mathcal{B}_p^2}\left\{\mathcal{A}_p[(U^T\hat{c}_{1,k})_p]^2 + [(U^T\hat{c}_{2,k})_p]^2\right\}\right], \tag{12}$$
where:
$$\mathcal{A}_p = -2A + 4\sum_{i=1}^3\sin^2\!\!\left(\frac{p_i}{2}\right) + 2\left[1-\cosh(\mu)\cos(p_4)\right], \tag{13a}$$
$$\mathcal{B}_p = 2\sinh(\mu)\sin(p_4). \tag{13b}$$

**Complexity reduction**: from $\mathcal{O}(K^{N^4})$ (naive product over all kernel combinations) to $\mathcal{O}(KN^4)$ after choosing appropriate symmetric centre distributions.

**RBF approximated partition function** (simplified form for symmetric centres around 0):
$$Z^\text{rbf} = \prod_p\sum_{k=1}^K a_k\int\mathcal{D}\tilde\Phi_p\exp\!\left[-\frac12\sum_{a,b}\tilde\Phi_{a,p}\tilde{W}_{ab,p}\tilde\Phi_{b,p} + 2A\sum_a c_{a,k_p}\tilde\Phi_{a,p} - A\sum_a c_{a,k_p}^2\right]. \tag{15}$$

---

## 5. Observables and effective potential

**Number density** via numerical derivative:
$$\langle n\rangle = \frac{1}{V}\frac{\partial\ln Z^\text{rbf}}{\partial\mu} \approx \frac{\ln Z^\text{rbf}(\mu+\Delta\mu) - \ln Z^\text{rbf}(\mu)}{\Delta\mu}, \quad \Delta\mu = 0.001. \tag{19}$$

**Effective potential** $V_\text{eff}(\phi_1)$ from Legendre transform. Add source $J_a\phi_{a,x}$ to $S_{1,x}$ (eq 16), set $J_2=0$ by O(2) symmetry. Classical field:
$$\langle\phi_1\rangle_J = \frac{\partial\ln Z^\text{rbf}[J_1]}{\partial J_1}\bigg|_{J_1=J}. \tag{18}$$
Effective potential:
$$V_\text{eff}(\phi_1) \propto \int d\langle\phi\rangle\, J(\langle\phi_1\rangle). \tag{17}$$

**Critical chemical potential**: fit $J(\langle\phi_1\rangle) = a\langle\phi_1\rangle^3 + b\langle\phi_1\rangle$; $\mu_c$ is determined by the zero-crossing of the coefficient $b$ as a function of $\mu$. Result: $\mu_c = 1.17 \pm 0.018$.

---

## 6. Algorithm (Julia implementation outline)

```julia
# Step 1: Train RBF network at each lattice site
function fit_rbf(m, λ, K, A, c_grid, N_T)
    # Sample N_T training points in field range
    ϕ_train = rand(Uniform(-2.5, 2.5), N_T, 2)
    y = [exp(-S1(m, λ, ϕ)) for ϕ in eachrow(ϕ_train)]
    # Build coefficient matrix H[i,k] = exp(-A * sum((ϕ_i - c_k)^2))
    H = [exp(-A * sum((ϕ_train[i,:] - c_grid[k,:]).^2)) for i=1:N_T, k=1:K]
    # Solve least squares with SVD regularisation
    a = H \ y   # or: svd-based with truncation
    return a
end

# Step 2: Build W matrix and diagonalise in momentum space
# W = M(μ) - 2A*I  (tridiagonal structure in momentum space)
# Compute 𝒜_p and ℬ_p for each momentum p

# Step 3: Compute ln Z^rbf via eq. 12 (product over momenta, sum over k)
# Use log-sum-exp trick for numerical stability

# Step 4: Compute ⟨n⟩ = ∂ ln Z / ∂μ by finite difference
```

**Key Julia notes**:
- Momentum sum: nested loop over $p \in [-N/2, N/2]^4$, total $N^4$ terms
- Use `FFTW.jl` for the $U$ similarity transformation
- `log-sum-exp` trick essential to avoid overflow in product form
- For model selection: compute RMSE on test set for each $(K, A, C)$

---

## Notation summary

| Symbol | Meaning |
|--------|---------|
| $S_{0,x}(\mu)$ | $\mu$-dependent part (kinetic + chemical potential, source of sign problem) |
| $S_{1,x}(m,\lambda)$ | self-interaction part, approximated by RBF |
| $a_k$ | RBF weight coefficients (from least squares) |
| $A$ | Gaussian width parameter (constant across kernels) |
| $c_{a,k}$ | centre of the $k$-th kernel in direction $a$ |
| $W = M(\mu)-2\bar{A}$ | matrix diagonalised by $U$ to factorise $Z$ |
| $\mathcal{A}_p, \mathcal{B}_p$ | real/imaginary parts of inverse propagator at momentum $p$ (eq 13) |
| $Z^\text{rbf}$ | RBF-approximated partition function |
| $\mu_c = 1.17 \pm 0.018$ | critical chemical potential for Bose condensation |

---

## Open questions

- Extension to fermionic theories (Grassmann fields complicate $S_{0,x}$)
- Scaling of required $K$ with volume and coupling $\lambda$
- Error control for observables beyond $\ln Z$
- Generalisation to gauge theories

---

## Links

- [[sign-problem]]
- [[silver-blaze-problem]]
- [[rbf-path-integral]]
- [[partition-function-lattice]]
- [[complex-langevin-dynamics]]
