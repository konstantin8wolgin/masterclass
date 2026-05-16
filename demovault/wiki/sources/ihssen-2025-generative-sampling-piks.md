---
title: "Generative sampling with physics-informed kernels"
authors: [Friederike Ihssen, Renzo Kapust, Jan M. Pawlowski]
year: 2025
type: paper
tags: [physics-informed-kernels, generative-network, renormalisation-group, lattice-qft, ood-problem, wegner-equation]
source_file: raw/papers/2510.26678v1.pdf
---

## Summary

Introduces the physics-informed kernel (PIK) architecture for generative Monte Carlo sampling. The architecture is built on PIRG flows: given an analytically prescribed path $S_t$ for the action ($t \in [0,1]$), the kernel $\dot\phi_t$ satisfying the **Wegner equation** is determined by solving an independent linear system at each $t$. This decouples training across layers — no backpropagation through time — and provides systematic OOD error control via the same linear PDE.

---

## 1. Setup and notation

The target distribution lives on a (hyper-)cubic lattice $\mathcal{D} = [0,\ldots,N_\text{max}]^d$ with
$$\phi : \mathcal{D} \to \mathbb{R}^\mathcal{D}, \qquad p(\varphi) = \frac{1}{\mathcal{N}} e^{-\hat{S}(\varphi)}, \qquad S = \hat{S} + \log\mathcal{N}.$$
Lattice volumes are typically $N_\text{max} \approx 10^5$ ($d=1,2,3$) or $N_\text{max} \approx 10^2$ ($d=4$).

The architecture starts from a simple base distribution $p_0(\varphi)$ at $t=0$ (e.g. a free Gaussian) and flows to the target $p_1(\phi) = p(\phi)$ at $t=1$.

---

## 2. PIRG flow equation

The central object is the **continuity / Wegner flow equation** for the flowing probability density $p_t(\phi) = e^{-S_t(\phi)}$:
$$\frac{dp_t(\phi)}{dt} = -\frac{\partial}{\partial\phi_i}\!\left[\dot\phi_{t,i}(\phi)\, p_t(\phi)\right], \qquad t \in [0,1]. \tag{4}$$

This is a total derivative: it preserves the partition sum $Z = \int \mathcal{D}\phi\, p_t(\phi) = \text{const}$ (in the absence of boundary terms). The kernel $\dot\phi_t : I \times \mathcal{C} \to \mathcal{C}$ is the infinitesimal map propagating the field:
$$\frac{d\phi_t}{dt} = \dot\phi_t(\phi_t), \qquad \phi_{t=0} = \varphi. \tag{11}$$

---

## 3. Wegner equation — the core PDE

Substituting $p_t = e^{-S_t}$ into (4) yields the **Wegner equation**:
$$\frac{dS_t(\phi)}{dt} + \dot\phi_t(\phi)\,\frac{\partial S_t(\phi)}{\partial\phi} = \frac{\partial\dot\phi_t(\phi)}{\partial\phi}. \tag{14}$$

**Key insight**: if the path $S_t$ is fixed analytically, this is a **linear first-order PDE for $\dot\phi_t$** alone. The $S_t$-dependent coefficients on the left are analytically known, so the equation decouples across different $t$-values.

The right-hand side is the divergence of $\dot\phi_t$; the left-hand side contains the time derivative of $S_t$ and the directional derivative of $S_t$ along $\dot\phi_t$.

---

## 4. Parametrisation of $S_t$

Split the action as:
$$S_t(\phi) = \hat{S}_t(\phi) + \log\mathcal{N}_t, \tag{16a}$$
where $\hat{S}_t$ is field-dependent and $\log\mathcal{N}_t$ is a field-independent normalisation. Parametrise:
$$\hat{S}_t(\phi) = \sum_i c_{i,t}\,\mathcal{O}_{i,t}(\phi), \tag{16b}$$
with analytically known basis functions $\mathcal{O}_{i,t}$, and scale-dependent coefficients $c_{i,t}$. At the boundaries, $\{c_i, \mathcal{O}_i\}$ must match the base $S_0$ and target $S_1$.

The normalisation flow satisfies:
$$\frac{d\log\mathcal{N}_t}{dt} = -\int\mathcal{D}\phi\, p_t(\phi)\,\frac{d\hat{S}_t(\phi)}{dt}. \tag{17}$$
This integral is expensive; in practice it is eliminated by taking **differences of Wegner equations** (see Section 5 below).

---

## 5. Parametrisation of $\dot\phi_t$ and the linear system

Expand the kernel in a chosen basis $\{K_{j,t}(\phi)\}$:
$$\dot\phi_t(\phi) = \sum_j k_{j,t}\,K_{j,t}(\phi). \tag{18}$$

Substituting into the Wegner equation (14) yields a **linear system**:
$$A_t\, k_t = b_t, \tag{19}$$
where the matrix $A_t$ and vector $b_t$ are determined by the basis $K_{j,t}(\phi)$ and the parametrised $S_t(\phi)$.

Because (19) is linear, finding the kernel reduces to **solving an independent linear system per time step** — local minima are not an issue. This is the core computational task in a PIK implementation.

### Avoiding explicit $d\log\mathcal{N}_t/dt$

To avoid evaluating (17), take the Wegner equation (14) at a target field $\phi$ and subtract the same equation evaluated at a reference configuration $\chi$. The field-independent $d\log\mathcal{N}_t/dt$ cancels, giving:
$$0 = \frac{d\hat{S}_t(\phi)}{dt} - \frac{d\hat{S}_t(\chi)}{dt}
  - \left[\frac{\partial}{\partial\phi} - \frac{\partial\hat{S}_t(\phi)}{\partial\phi}\right]\dot\phi_t(\phi)
  + \left[\frac{\partial}{\partial\phi} - \frac{\partial\hat{S}_t(\phi)}{\partial\phi}\right]\dot\phi_t(\phi)\bigg|_{\phi=\chi}. \tag{20}$$

Choosing $\chi = 0$ gives a vanishing time and field derivative of the action at $\chi$, simplifying this further. This equation contains only analytically computable terms.

---

## 6. OOD correction (residual / collocation)

Define the **linear operator**:
$$\mathcal{L}_\chi\,\dot\phi_t(\phi) := \left[\frac{\partial}{\partial\phi} - \frac{\partial\hat{S}_t(\phi)}{\partial\phi}\right]\dot\phi_t(\phi) - \left[\frac{\partial\hat{S}_t(\phi)}{\partial\phi} - \frac{\partial}{\partial\phi}\right]\dot\phi_t(\phi)\bigg|_{\phi=\chi}. \tag{26}$$

Given an approximate kernel $\dot\phi_t^{(0)}$, compute its **residual**:
$$r_\chi^{(0)}(\phi) = \frac{d\hat{S}_t(\phi)}{dt} - \frac{d\hat{S}_t(\chi)}{dt} - \mathcal{L}_\chi\,\dot\phi_t^{(0)}(\phi). \tag{27}$$

Find a **correction** $\Delta\dot\phi_t^{(0)}$ satisfying:
$$\mathcal{L}_\chi\,\Delta\dot\phi_t^{(0)}(\phi) = r_\chi^{(0)}(\phi). \tag{28}$$

This is again a linear PDE, solved by collocation with compactly supported radial basis functions (Wendland functions, see Section 7). The corrected kernel is:
$$\dot\phi_t^{(1)}(\phi) = \dot\phi_t^{(0)}(\phi) + \Delta\dot\phi_t^{(0)}(\phi). \tag{29}$$

This procedure can be iterated — the linear structure guarantees systematic improvement.

---

## 7. Collocation correction with Wendland RBF

The correction term $\Delta\dot\phi_t^{(0)}$ is expanded as (Appendix B):
$$\Delta\dot\phi_t^{(0)}(\phi) = \sum_{i=1}^M w_i\,W_i(\phi,\phi_i). \tag{B2}$$

The basis functions are **anti-symmetrised Wendland functions** (to preserve $\mathbb{Z}_2$ symmetry):
$$W_i(\phi,\phi_i) = W\!\left(\frac{|\phi-\phi_i|}{\sigma}\right) - W\!\left(\frac{|\phi+\phi_i|}{\sigma}\right), \tag{B4}$$
where $W$ is the 4-times differentiable Wendland function:
$$W(r) = (1-r)_+^{l+2}\left[(l^2+4l+3)r^2 + (3l+6)r + 3\right], \quad l = \lfloor d/2\rfloor + 3. \tag{B3}$$
Here $(1-r)_+ = \max(1-r,0)$ ensures compact support. Shape parameter: $\sigma \approx 0.41$.

Centres $\{\phi_i\}_{i=1}^M$ are chosen in regions where the approximate kernel deviates most from the true one ($M=11$ in the paper). The weights $w_i$ are found by solving the overdetermined linear system at $N=501$ collocation points via QR decomposition.

---

## 8. Zero-dimensional benchmark ($\phi^4$ theory, $d=0$)

**Target action**:
$$\hat{S}(\phi) = \frac{1}{2}m^2\phi^2 + \frac{\lambda}{4!}\phi^4, \qquad \phi \in \mathbb{R}. \tag{21}$$
Parameters used: base $(m_0^2,\lambda_0) = (1, 0)$ (unimodal Gaussian), target $(m_1^2,\lambda_1) = (-2, 1/2)$ (bimodal).

**Linear interpolation path**:
$$\hat{S}_t(\phi) = \frac{1}{2}m^2(t)\phi^2 + \frac{\lambda(t)}{4}\phi^4, \tag{22a}$$
$$m^2(t) = m_0^2 + t(m_1^2 - m_0^2), \quad \lambda(t) = \lambda_0 + t(\lambda_1 - \lambda_0). \tag{22b}$$

**Basis functions for kernel** (exploiting $\mathbb{Z}_2$ symmetry, $\hat{S}$ symmetric under $\phi\to-\phi$):
$$K_{1,t}(\phi) = \phi, \quad K_{j,t}(\phi) = \sin(\omega_{j,t}\phi) \text{ for } j \geq 2. \tag{25}$$
$N_f = 25$ learnable frequencies $\omega_{j,t}$; coefficients $k_{j,t}$ from solving (20).

**ODE integration**: 4th-order Runge-Kutta with step size $\Delta t = 1/25$, giving $N_t = 51$ time slices (kernels).

**Training**: Adam optimiser, loss = mean squared error of (20) at $\chi = 0$. Learning rate:
- First two kernels: $\eta = 10^{-2}$, 600 epochs, batch size 512
- Subsequent kernels: $\eta = 10^{-3}$, 600 epochs (weights initialised from previous kernel)
- Loss bounded by $10^{-4}$ at saturation.

**Sampling**: draw $N_s = 10^6$ samples from $p_0$ (Gaussian), push through RK4 flow
$$\phi_t = \phi_0 + \int_0^t ds\,\dot\phi_s(\phi_s).$$
Observable at any $t$: $\langle\mathcal{O}\rangle_{p_t} \approx \frac{1}{N_s}\sum_n \mathcal{O}(\phi_t(\phi_0^{(n)}))$.

Speedup vs. HMC (rough estimate): ~20× more independent samples per wall-clock time.

---

## 9. Higher-dimensional $\phi^4$ theory ($d>0$)

Action on $d$-dimensional lattice $\mathcal{D}$ with $V$ sites:
$$\hat{S}_t(\phi) = \sum_{\hat{n}\in\mathcal{D}}\left[\frac{1}{2}(m^2(t)+2d)\phi_{\hat{n}}^2 - \sum_{\mu=1}^d\phi_{\hat{n}}\phi_{\hat{n}+\hat\mu} + \frac{\lambda(t)}{4}\phi_{\hat{n}}^4\right]. \tag{31}$$

**Helmholtz decomposition of the kernel** (reducing $V$ functions to 1):

For higher-dimensional systems, write the kernel as a gradient field plus a divergence-free part:
$$\dot\phi_t(\phi) = -\frac{\partial}{\partial\phi}\Omega_t(\phi) + H_t(\phi), \qquad \frac{\partial}{\partial\phi}[p_t(\phi)H_t(\phi)] = 0. \tag{32}$$

Since $H_t$ does not contribute to the change in $S_t$, one solves the **adapted Wegner equation** for $\Omega_t$ alone:
$$\left[\frac{\partial^2}{\partial\phi^2} - \frac{\partial S_t(\phi)}{\partial\phi}\frac{\partial}{\partial\phi}\right]\Omega_t(\phi) = -\frac{dS_t(\phi)}{dt}. \tag{34}$$

This is a **linear PDE of elliptic type** (canonical Laplacian minus a drift term). Boundary condition:
$$\lim_{R\to\infty}\int_{\partial\mathcal{V}}\mathcal{D}\phi_s\,p_t(\phi_s)\dot\phi_t(\phi_s) = 0, \tag{36}$$
ensuring normalisation. The kernel is then recovered via $\dot\phi_t = -\partial_\phi\Omega_t$.

**Scaling advantage**: for local actions (nearest-neighbour couplings), the operator in (34) maps tensor products of nearest-neighbour terms to themselves, so the number of required basis functions scales with the *spacetime dimension* $d$ rather than the volume $V$.

---

## 10. PIKs as finite-time Langevin corrections (Appendix C)

Connection to Langevin: the Langevin SDE
$$\frac{d\phi_t}{dt} = f_t(\phi) + \sqrt{2}\,\eta_t$$
drives $p_t$ according to
$$\frac{dS_t}{dt} = \frac{\partial}{\partial\phi}f_t - f_t\frac{\partial S_t}{\partial\phi} + \frac{\partial^2 S_t}{\partial\phi^2} - \left(\frac{\partial S_t}{\partial\phi}\right)^2. \tag{C2}$$

Setting $f_t = -\frac{\partial S_t}{\partial\phi} + \dot\phi_t$ and substituting into (C2) retrieves the Wegner equation (14). PIKs thus correspond to **finite-time corrections to the Langevin drift** that force thermalisation to complete in exactly $T=1$ rather than $t\to\infty$.

---

## Algorithm summary (Julia implementation outline)

```
# Inputs: S0 (base action), S1 (target action), N_t (time steps), N_s (samples)
# 1. Define interpolating path: S_t(ϕ, t) via linear or custom schedule
# 2. For each time step t_i = i/N_t  (i = 0, ..., N_t):
#    a. Build matrix A_{t_i} and vector b_{t_i} from basis K_{j,t}(ϕ)
#       evaluated at training configs (or collocation points)
#    b. Solve linear system: k_{t_i} = A_{t_i} \ b_{t_i}
#    c. Store kernel: ϕ̇_{t_i}(ϕ) = Σ_j k_{j,t_i} K_{j,t_i}(ϕ)
#    d. [Optional] Compute residual r; solve correction ΔΦ; update kernel
# 3. Draw N_s samples ϕ_0 from p_0 (Gaussian)
# 4. Integrate field ODE: ϕ_{t+Δt} = RK4(ϕ_t, ϕ̇_t, Δt)
# 5. Evaluate observables: ⟨O⟩ ≈ mean(O.(ϕ_1))
```

**Key Julia operations**:
- `A \ b` (backslash) for the per-layer linear system — sparse if basis is local
- `DifferentialEquations.jl` (RK4) for ODE integration of samples
- Basis $K_{j,t}$: `[ϕ, sin.(ω[j]*ϕ) for j in 2:N_f]` with learnable `ω`
- For $d>0$: solve elliptic PDE (34) — can use `LinearAlgebra.cholesky` if operator is positive definite

---

## Notation summary

| Symbol | Meaning |
|--------|---------|
| $\phi_t(\phi_0)$ | field configuration at RG-time $t$, starting from $\phi_0$ |
| $\dot\phi_t(\phi)$ | physics-informed kernel — infinitesimal map at time $t$ |
| $S_t(\phi)$ | negative log-likelihood at time $t$; $S_0$ = base, $S_1$ = target |
| $\hat{S}_t, \log\mathcal{N}_t$ | field-dependent and normalisation parts of $S_t$ |
| $\mathcal{O}_{i,t}(\phi)$ | basis functions for $\hat{S}_t$ |
| $K_{j,t}(\phi)$ | basis functions for $\dot\phi_t$ |
| $k_{j,t}$ | kernel expansion coefficients (solution of $A_t k_t = b_t$) |
| $\theta_t$ | all optimisable parameters of $\dot\phi_t$ at time $t$ (independent of $\theta_{t'}$) |
| $\chi$ | reference configuration (typically $\chi=0$) for difference Wegner equation |
| $\Omega_t(\phi)$ | scalar potential for Helmholtz-decomposed kernel in $d>0$ |

---

## Open questions

- Efficient basis choice for interacting theories beyond scalar $\phi^4$
- Scaling to very large lattice volumes (basis size vs. $V$)
- Global optimisation of the pair $(F_t, \dot\phi_t)$ — currently local per layer
- Extension to fermionic theories and gauge fields
- Parameter-conditional flows for scanning over $(m^2, \lambda)$ simultaneously

---

## Links

- [[physics-informed-kernels]]
- [[physics-informed-rg-flows]]
- [[wegner-equation]]
- [[sign-problem]]
- [[ihssen-2026-solving-sign-problems-piks]]
