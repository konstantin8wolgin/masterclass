---
title: "Solving sign problems with physics-informed kernels"
authors: [Friederike Ihssen, Renzo Kapust, Jan M. Pawlowski]
year: 2026
type: paper
tags: [physics-informed-kernels, sign-problem, complex-distribution, pikfold, zero-dimensional, real-time]
source_file: raw/papers/2603.03159v1.pdf
---

## Summary

Extends the PIK architecture (see [[ihssen-2025-generative-sampling-piks]]) to systems with complex probability distributions, addressing the sign problem. The core property is **weight-preservation**: the PIK flow maps a complex distribution $p(\varphi) \in \mathbb{C}$ to a real, positive distribution on a deformed manifold $\mathcal{M}_1$ (a PIKfold), on which standard Monte Carlo is efficient. Demonstrated on zero-dimensional $\phi^4$ theory with complex action and on the real-time harmonic oscillator, including cases where Complex Langevin and Lefschetz thimbles fail.

---

## 1. Complex sampling problem

Consider a complex scalar field on a lattice $\mathcal{D}$ with complexified field space:
$$\phi : \mathcal{D} \to \mathcal{C} = \mathbb{C}^\mathcal{D}. \tag{6}$$
The distribution is:
$$p(\varphi) = e^{-S(\varphi)} \in \mathbb{C}, \qquad S(\varphi) = \hat{S}(\varphi) + \log\mathcal{N}, \quad \hat{S}(\varphi=0) = 0. \tag{1a}$$
The partition function and measure:
$$Z = \int_\mathcal{M} d\varphi\, p(\varphi) = 1, \qquad \int_\mathcal{M} d\varphi = \prod_i \int_\mathbb{R} d\varphi_i. \tag{1b/c}$$

Since $p \in \mathbb{C}$, importance sampling is impossible directly. A **sampling pair** $(p_0(\phi), \mathcal{M}_0)$ needs to be found where sampling is efficient and sign-free.

---

## 2. PIK map and weight-preservation

The PIK architecture maps one sampling pair to another:
$$(p_0(\phi), \mathcal{M}_0) \longrightarrow (p_1(\phi), \mathcal{M}_1), \qquad p_1(\phi) = p(\phi). \tag{2}$$

The key property: the map is **weight-preserving**,
$$\int_{\mathcal{M}_1} d\mu_1\, p_1(\phi)\,\mathcal{O}(\phi) \stackrel{!}{=} \int_\mathcal{M} d\mu\, p(\phi)\,\mathcal{O}(\varphi(\phi)), \tag{3}$$
for general operators $\mathcal{O}$. Here $d\mu_t$ is the volume element of $\mathcal{M}_t$.

This means: if the sampling pair at $t=0$ is sign-problem free, **all** pairs $(p_t, \mathcal{M}_t)$ are sign-problem free.

---

## 3. Flow equation and complex Wegner equation

The flowing distribution $p_t(\phi) = e^{-S_t(\phi)}$ satisfies:
$$\frac{dp_t(\phi)}{dt} = -\frac{\partial}{\partial\phi_i}\!\left[\dot\phi_{t,i}(\phi)\, p_t(\phi)\right], \quad t \in [0,1], \tag{4}$$
where index $i \in \mathcal{D}$ is summed over data points. The $t$-dependent distribution is $p_t = e^{-S_t}$ with $S_t = \hat{S}_t + \log\mathcal{N}_t$.

The corresponding **Wegner equation** (same form as for real PIKs):
$$\frac{dS_t(\phi)}{dt} + \dot\phi_{t,i}(\phi)\,\frac{\partial S_t(\phi)}{\partial\phi_i} = \frac{\partial\dot\phi_{t,i}(\phi)}{\partial\phi_i}. \tag{7}$$

For complex $dS_t/dt$ and/or $\partial S_t/\partial\phi$, the PIK $\dot\phi(\phi)$ is complex and the differential manifold $\partial_t\mathcal{M}_t \subset \mathbb{C}^\mathcal{D}$ is complex-valued.

---

## 4. PIKfold condition

The **PIKfold** $\mathcal{M}_1$ is the final sampling manifold on which $p_1(\phi)\,d\mu_1$ is real and positive. For $\mathcal{M}_0 = \mathbb{R}$ (real initial sampling), the PIKfold is determined by the **weight-preserving ODE** (1D case):
$$\frac{d\phi_I}{d\phi_R} = -\frac{\text{Im}[p_t(\phi)]}{\text{Re}[p_t(\phi)]} = \tan\!\left(\text{Im}[S_t(\phi_R + i\phi_I)]\right). \tag{16}$$

Here $\phi = \phi_R + i\phi_I$ are the real and imaginary coordinates. Equation (16) defines $\mathcal{M}_t$ as a monotonic function $\phi_I(\phi_R)$ covering $\mathbb{R}$.

**Existence criterion**: a real PIKfold exists if and only if equation (16) has a well-defined solution at $t=1$. Non-existence at $t=1$ means no real PIKfold for the given initial pair — one must choose a different $(p_0, \mathcal{M}_0)$ or work with generalised PIKfolds.

**Boundary condition**: $\text{Re}[S_t(\phi)] > 0$ as $|\phi|\to\infty$ on $\mathcal{M}_t$, ensuring convergence of the contour integral.

---

## 5. Sampling with the complex PIK

The RG-time evolution of samples is:
$$\phi_t(\phi_0) = \phi_0 + \int_0^t ds\,\dot\phi_s(\phi_s(\phi_0)). \tag{8b / 34}$$

The observable estimator at any $t$ using $N_s$ samples from $p_0$:
$$\langle\mathcal{O}\rangle_{p_t} \approx \frac{1}{N_s}\sum_{n=1}^{N_s}\mathcal{O}(\phi_t(\phi_0^{(n)})). \tag{36}$$

At $t=1$, this gives the observable with respect to the target distribution $p(\varphi)$.

---

## 6. Zero-dimensional $\phi^4$ with complex coupling

**Target action** (main text, single-site):
$$\hat{S}_t(\phi) = \frac{1}{2}\phi^2 + \frac{1}{4}\phi^4 + (1+it)\phi, \quad \phi \in \mathbb{C}. \tag{17}$$
At $t=0$: real $\hat{S}_0 = \frac12\phi^2 + \frac14\phi^4 + \phi$ (real sampling possible). At $t=1$: complex action $\hat{S}_1 = \frac12\phi^2+\frac14\phi^4+(1+i)\phi$.

The initial sampling manifold is $\mathcal{M}_0 = \mathbb{R}$; a real PIKfold $\mathcal{M}_1$ exists for this action (verified numerically). The PIKfold is depicted in paper Fig. 1 (black line).

**Supplement S3.2 — complex $\mathbb{Z}_2$-symmetric $\phi^4$**: interpolating path for testing complex mass:
$$\hat{S}_t(\phi) = m_t^2\frac{\phi^2}{2} + \lambda\frac{\phi^4}{4}, \quad m_t^2 = m_0^2 + \Delta m^2 t, \tag{65}$$
with $m_0^2 = \lambda = 1$ and $\Delta m^2 \in \mathbb{C}$.

**Free theory with complex current** (analytically solvable):
$$\hat{S}_t(\phi) = \frac{1}{2}m^2\phi^2 + it\phi. \tag{58}$$
The Wegner equation gives an analytic PIK:
$$\dot\phi_t(\phi) = -\frac{i}{m^2}, \qquad \phi_t(\phi_0) = \phi_0 - \frac{i}{m^2}t. \tag{62}$$
This is a pure complex shift — the field moves parallel to the imaginary axis.

**Normalisation**: $\frac{d\log\mathcal{N}_t}{dt} = -\frac{t}{m^2}$, so $\log\mathcal{N}_t = -\frac{t^2}{2m^2} + \log\mathcal{N}_0$, with $\mathcal{N}_0 = \sqrt{2\pi/m^2}$.

---

## 7. Real-time harmonic oscillator

Goal: compute real-time (Minkowski) observables from Euclidean samples by rotating couplings.

**Action path** (from Euclidean $t=0$ to Minkowski $t=1$):
$$\hat{S}_t(\phi) = \sum_{i\in\mathcal{D}}\frac{\kappa_t}{2}(\phi_{i+1}-\phi_i)^2 - \left(m_t^2+\epsilon\right)\frac{\phi_i^2}{2}, \tag{19}$$
with time-dependent couplings:
$$\kappa_t = e^{-i\frac{\pi}{2}t}, \qquad m_t^2 = m^2\, e^{+i\frac{\pi}{2}t}. \tag{20}$$

At $t=0$: $\kappa_0 = 1$, $m_0^2 = m^2$ → standard Euclidean harmonic oscillator.
At $t=1$: $\kappa_1 = -i$, $m_1^2 = im^2$ → Minkowski signature ($e^{iS_M}$).

Small $\epsilon > 0$ is a mass regularisation ensuring finiteness of the path integral.

**Mass matrix notation**: write $\hat{S}_t(\phi) = \frac{1}{2}\phi_i M_{t,ij}\phi_j$ where $M_t$ is the $N_\text{max}\times N_\text{max}$ tridiagonal matrix with
- diagonal: $-m_t^2 - \epsilon + 2\kappa_t$ (or $+m_t^2+\epsilon$ depending on sign convention — check paper eq 77)
- off-diagonal: $-\kappa_t$

**Analytic PIK for free theory**:
$$\dot\phi_t(\phi) = -\frac{1}{2}M_t^{-1}\frac{dM_t}{dt}\phi. \tag{21 / 81}$$

This is a **linear transformation of $\phi$**: the field is just multiplied by the matrix $-\frac{1}{2}M_t^{-1}\dot{M}_t$ at each time step. No nonlinear kernel needed.

**Numerical parameters**: $N_\text{max} = 30$ and $N_\text{max} = 75$ lattice sites, $m^2 = 0.3$, $\epsilon = \frac{2m}{N_\text{max}}$.

**Sampling procedure**:
1. Sample $10^5$ Euclidean configurations $\phi_0 \sim p_E(\varphi) \propto e^{-\hat{S}_E(\varphi)}$ at $t=0$
2. Integrate ODE $\phi_t = \phi_0 + \int_0^t ds\,\dot\phi_s(\phi_s)$ using Euler or RK4 with $N_t \geq 2^9$ steps (convergence criterion)
3. Evaluate $\langle\phi_0\phi_{x_0}\rangle$ at $t=1$ using the transported samples

**Time-stepping convergence**: the time-discretisation error for the correlation $\langle\varphi_0\varphi_1\rangle$ agrees with the exact solution for $N_t \geq 2^9$ equidistant steps; this does not scale with $N_\text{max}$ (spatial lattice size).

---

## 8. General dimensions — free scalar theory (Supplement S4.1)

For a free scalar on a $d$-dimensional lattice with $|\mathcal{D}|$ sites:
$$\hat{S}_t(\phi) = \sum_{i\in\mathcal{D}}\left[T_i(\phi) + V_i(\phi)\right] = \frac{1}{2}\sum_\mu\kappa_{t,\hat\mu}(\phi_{i+\hat\mu}-\phi_i)^2 + \left(m_t^2+\epsilon\right)\frac{\phi_i^2}{2} = \frac{1}{2}\phi_i M_{t,ij}\phi_j. \tag{77}$$

Couplings for Euclidean→Minkowski rotation:
$$\kappa_{t,\hat\mu} = e^{-i\kappa_{M,\hat\mu}\frac{\pi}{2}t}, \quad m_t^2 = m^2 e^{+i\frac{\pi}{2}t}, \tag{78}$$
where $\kappa_{M,0} = 1$ (time direction) and $\kappa_{M,i} = -1$ for $i=1,\ldots,d-1$ (spatial).

**Normalisation**: $\mathcal{N}_t = \sqrt{(2\pi)^{|\mathcal{D}|}/\det M_t}$.

**Wegner equation for this theory** (eq 80):
$$\frac{\partial\dot\phi_{t,i}}{\partial\phi_i} - \dot\phi_{t,i}(\phi)\,M_t^{ij}\phi_j = \frac{1}{2}\phi_i\frac{dM_t^{ij}}{dt}\phi_j + \frac{d\log\mathcal{N}_t}{dt}.$$

**Solution** (using $\frac{d\log\mathcal{N}_t}{dt} = -\frac{1}{2}\text{Tr}(M_t^{-1}\dot{M}_t)$ via Jacobi's formula):
$$\dot\phi_t(\phi) = -\frac{1}{2}M_t^{-1}\frac{dM_t}{dt}\phi. \tag{81}$$

This requires solving $M_t\dot\phi_t = -\frac{1}{2}\dot{M}_t\phi$ — a linear system per time step, exactly the same structure as the real PIK but now complex.

---

## 9. Real-time proxy $\phi^4$ (Supplement S3.4)

Zero-dimensional case with complex RG-time path (as test of sign problem capability):
$$\hat{S}_t(\phi) = m_t^2\frac{\phi^2}{2} + \lambda_t\frac{\phi^4}{4}, \quad m_t^2 = m^2 e^{-i\frac{\pi}{2}t}, \quad \lambda_t = \lambda\, e^{-i\frac{\pi}{2}t}. \tag{68}$$

At $t=0$: $\hat{S}_0 = \hat{S}_E$ (real Euclidean). At $t=1$: $\hat{S}_1 = -i\hat{S}_M$ (Minkowski prefactor $e^{iS_M}$).

Free theory version ($\lambda=0$) has analytic solution:
$$\dot\phi_t = \frac{i\pi}{4}\phi, \qquad \phi_t(\phi_0) = \phi_0\, e^{i\frac{\pi}{4}t}. \tag{71}$$
Jacobian: $\partial_{\phi_0}\phi_t = e^{i\frac{\pi}{4}t}$, with $\log\mathcal{N}_t = -\frac{\pi}{4}it + \log\mathcal{N}_0$.

---

## 10. Comparison to other methods

| Method | Key equation | Advantage | Limitation |
|--------|-------------|-----------|-----------|
| Complex PIK | Wegner eq. (7) + weight-preservation | Constructive, analytic, no boundary terms | PIKfold existence not guaranteed |
| Complex Langevin (CL) | $dx/dt = -\partial S/\partial x + \eta$ | Simple to implement | Convergence not guaranteed; fails for real-time oscillator |
| Lefschetz thimbles | $\partial_t\tilde\phi = -\overline{(\partial S/\partial\tilde\phi)}$ | Eliminates intra-thimble sign problem | Residual inter-thimble (global) sign problem; Jacobian |
| Holomorphic flow | Same ODE as thimble | Deforms contour toward thimbles | Does not fully solve global sign problem |

CLE relation (Supplement S2.1): the PIK flow equation (4) with operator kernel $\dot\phi_\text{op}(\phi) = \dot\phi(\phi) + K(\phi)\frac{\partial}{\partial\phi}$ reduces to CLE when $\dot\phi = \partial S/\partial\phi$ and $K_{ij}(\phi) = \delta_{ij}$.

---

## Algorithm summary (Julia implementation)

```
# Real PIK for complex action — 0D benchmark
# Inputs: action path S_t(ϕ, t), N_t steps, N_s samples

# Step 1: Compute PIKfold M_1 (for each t, integrate ODE for ϕ_I)
#   dϕ_I/dϕ_R = -Im[p_t(ϕ)] / Re[p_t(ϕ)] = tan(Im[S_t(ϕ_R + i*ϕ_I)])
#   Solve this ODE numerically (e.g. RK4 in ϕ_R) for t=1

# Step 2: Sample ϕ_0 from p_0 on M_0=ℝ (real Gaussian or standard Monte Carlo)

# Step 3: For each time step t_i:
#   - Solve Wegner equation (7) numerically for dot_ϕ_{t_i}
#     (linear system: expand in basis, collect coefficients from training configs)
#   - Integrate field: ϕ_{t_{i+1}} = ϕ_{t_i} + Δt * dot_ϕ_{t_i}(ϕ_{t_i})
#     (use RK4 for better accuracy)

# Step 4: Evaluate observables on transported samples at t=1

# For free theory (harmonic oscillator): dot_ϕ_t(ϕ) = -0.5 * M_t \ (dM_t/dt * ϕ)
#   M_t is tridiagonal → use LinearAlgebra.lu or sparse solvers
```

**Key Julia operations**:
- Complex-valued arrays: `Array{ComplexF64}` for $\phi_t$, $S_t$, $M_t$
- ODE for PIKfold (16): `DifferentialEquations.jl` in 1D, or manual RK4
- Matrix solve for free theory: `M_t \ (dM_dt * ϕ)` where `M_t` is sparse tridiagonal
- For $N_t \geq 512$ steps: use adaptive integrator or fixed RK4

---

## Notation summary

| Symbol | Meaning |
|--------|---------|
| $p(\varphi) \in \mathbb{C}$ | complex probability distribution (sign problem source) |
| $\mathcal{M}_t$ | sampling manifold at RG-time $t$; starts real, deforms into $\mathbb{C}^\mathcal{D}$ |
| PIKfold | the manifold $\mathcal{M}_1$ on which $p_1(\phi)\,d\mu_1 \in \mathbb{R}_{>0}$ |
| $\mathcal{P} = \mathcal{P}_1$, $\mathcal{P}_t = \{\mathcal{M}_s | s\in[0,t]\}$ | set of (sub-)PIKfolds |
| $p_E(\varphi) \sim e^{-\hat{S}_E}$ | Euclidean distribution at $t=0$ |
| $p_M(\varphi) \sim e^{i\hat{S}_M}$ | Minkowski distribution at $t=1$ |
| $\kappa_t = e^{-i\pi t/2}$ | time-dependent coupling rotating Euclidean→Minkowski |
| $M_t$ | mass matrix; $\hat{S}_t = \frac12\phi M_t\phi$ for free theory |

---

## Open questions

- Extension to finite-density fermionic theories (complex fermion determinant)
- Higher-dimensional interacting models (infinitely many cumulants)
- Existence criteria for PIKfolds in general gauge theories
- Whether complex PIK replaces CL in all cases where CL fails
- Full numerical details of complex Wegner equation solver (companion paper in preparation)

---

## Links

- [[physics-informed-kernels]]
- [[physics-informed-rg-flows]]
- [[sign-problem]]
- [[wegner-equation]]
- [[complex-langevin-dynamics]]
- [[lefschetz-thimbles]]
- [[ihssen-2025-generative-sampling-piks]]
