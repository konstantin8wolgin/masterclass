---
title: "Deep Learning of Fermion Sign Fluctuations"
authors: [Scott Lawrence, Yukari Yamauchi]
year: 2023
type: paper
tags: [sign-problem, subtraction, neural-network, thirring-model, lattice-qft]
source_file: raw/papers/2212.14606v2.pdf
---

## Summary

Alleviates the fermion sign problem by subtracting a divergence $g = \nabla\cdot v$ from the Boltzmann factor $e^{-S}$. Any such subtraction integrates to zero (divergence theorem), leaving $Z$ and all observables unchanged but improving the average phase. When analytic ansätze are insufficient, a neural network parametrises the vector field $v$. Demonstrated on the staggered Thirring model in $1+1$D. Key result: every infinitesimal contour deformation is a subtraction; the converse is false.

**When to use**: when the sign problem is mild-to-moderate and one wants a simple, provably-correct modification of the Boltzmann factor without deforming the integration manifold. Neural network approach for non-perturbative $\mu$.

---

## 1. Subtraction method

**Core idea**: find $g(\phi)$ with $\int\mathcal{D}\phi\, g(\phi) = 0$, then define subtracted Boltzmann factor:
$$e^{-S_\text{sub}} \equiv e^{-S} - g(\phi). \tag{12 / general}$$

Guaranteed construction: set $g = \nabla\cdot v$ for any vector field $v$; then $\int g = 0$ by the divergence theorem (for compact fields or fast-decaying $v$).

**Scaled vector field**: set $v = e^{-S}u$ where $u$ is order 1 (not exponentially small). The subtracted Boltzmann factor:
$$e^{-S_u} \equiv e^{-S} - \nabla\cdot(e^{-S}u) = e^{-S}(1 + u\cdot\nabla S - \nabla\cdot u). \tag{13}$$

**Observable estimator** (avoids sign-to-noise in naive ratio):
$$\langle\mathcal{O}\rangle = \left\langle\mathcal{O} - \frac{e^{-S_u}}{e^{-S_v}}\,u\cdot\nabla\mathcal{O}\right\rangle_{S_u}. \tag{15}$$

Empirically avoids severe sign-to-noise (no rigorous proof given).

---

## 2. Connection to contour deformations

For a contour $\tilde\phi_i(\phi) = \phi_i - u_i(\phi)$ parametrised by the real plane, the subtraction induced by an infinitesimal deformation is:
$$g_\text{sub}(\phi) = e^{-S(\tilde\phi(\phi))}\det\frac{\partial\tilde\phi}{\partial\phi} - e^{-S(\phi)} = e^{-S(\phi)}\left[1 + u\cdot\nabla S - \nabla\cdot u\right] - e^{-S(\phi)} + O(u^2). \tag{16/18}$$

**This matches (13) to leading order in $u$**: every infinitesimal contour deformation is a subtraction.

**Converse fails**: not every subtraction is a contour deformation. Example: $S = \cos\theta + \epsilon$ — the sign problem is entirely removed by subtracting $\cos\theta$, but no contour deformation achieves this.

---

## 3. Thirring model

**Lattice action** (staggered, $1+1$D, $\beta\times L$ lattice):
$$S = \sum_{x,\nu=0,1}\frac{2}{g^2}(1-\cos A_\nu(x)) - \log\det K[A], \tag{3}$$
where $A_\nu(x) \in [0, 2\pi)$ are logarithms of $U(1)$-valued link fields.

**Dirac matrix $K[A]$**:
$$K[A]_{xy} = m\delta_{xy} + \frac{1}{2}\sum_{\nu=0,1}\left[\eta_\nu e^{iA_\nu(x)+\mu\delta_{\nu,0}}\delta_{x+\nu,y} - \eta_\nu e^{-iA_\nu(x')-\mu\delta_{\nu,0}}\delta_{x-\nu,y}\right], \tag{4}$$
with staggered phases $\eta_0 = (-1)^{\delta_{x_0,0}}$, $\eta_1 = (-1)^{x_0}$.

Full Dirac matrix entry (eq 6):
$$K[A]_{tt'} = \frac{1}{2}\left[e^{\mu+iA(t)}\delta_{(t+1)t'} - e^{-\mu-iA(t')}\delta_{(t'+1)t} - e^{\mu+iA(t)}\delta_{tN}\delta_{t'1} + e^{-\mu-iA(t')}\delta_{t1}\delta_{t'N}\right] + m\delta_{tt'}.$$

Bare parameters: $m=0.05$, $g^2=1.0$ (strong coupling). Sign problem from $\mu$-dependent imaginary terms.

Simple contour (shift): $\tilde{A}_0 = A_0 + i\mu$, $\tilde{A}_1 = A_1$. Removes $\mu$ from Dirac matrix completely; used as baseline.

---

## 4. Analytic subtractions

**Heavy-dense (large $\mu$) subtraction** — leading term of det $K$ via polymer expansion:
$$g_\text{hd}(A) = e^{\frac{2}{g^2}\sum_{x,\nu}\cos A_\nu(x)}\cdot 2^{-\beta L}e^{\beta L\mu+t\sum_x A_{0}(x)} - e^{\beta L\mu}(4\pi)^{-\beta L}\left[I_0(2/g^2)I_1(2/g^2)\right]^{\beta L}. \tag{22}$$

**Gradient descent subtraction**: functional derivative of $Z_{Q,u}$ w.r.t. $u$ gives:
$$\frac{\delta}{\delta u(\phi)}Z_{Q,u} = -i\nabla\,\text{Im}\,S. \tag{24}$$
One gradient step: $u \propto \frac{\delta}{\delta u}Z_{Q,u} = -i\nabla\,\text{Im}\,S$.

**Heat equation subtraction**: evolve $e^{-S_t(\phi)}$ under the heat equation $\partial_t e^{-S_t} = -\sum_i\partial_{\phi_i}^2 e^{-S_t}$ for one step $\Delta t$:
$$e^{-S} \to e^{-S} - \Delta t\sum_i\partial_{\phi_i}^2 e^{-S}. \tag{29}$$

---

## 5. Neural network parametrisation

Network represents the scaled vector field $u_i(\phi)$:
- **Inputs**: $\sin A_\nu(x)$, $\cos A_\nu(x)$ — total $4V$ inputs ($V$ = number of lattice sites)
- **Architecture**: densely connected multi-layer perceptron, each hidden layer width $4V$, CELU activation
- **Outputs**: $4V$ real numbers = real and imaginary parts of $2V$ components of complex $u$
- **Cost function** (maximise average phase):
  $$C[u(\lambda)] = -\langle\sigma\rangle_u, \qquad \frac{\partial C}{\partial\lambda} = -\langle\text{Re}\,S\rangle_Q. \tag{30/31}$$
- **Gradient** (quenched expectation, efficiently computed): $\frac{\partial}{\partial\lambda}\log Z_Q^{(\lambda)} = -\langle\text{Re}\,S\rangle_Q$.
- **Training**: ADAM with scheduled learning rate:
  $$\eta(n) = \begin{cases}10^{-4} + \frac{10^{-2}-10^{-4}}{200}n & 1\leq n\leq 200 \\ 10^{-2}\cdot 0.1^{(n-200)/(C\times 10^2)} & n > 200\end{cases}. \tag{32}$$
- At each step: 100 samples used for gradient estimate; samples re-used until average reweighting $R = \frac{1}{100}\sum_k|\log\frac{e^{-S(\phi_k)} - g(\phi_k)}{e^{-S(\phi_k)} - \tilde{g}(\phi_k)}| > 0.3$.

---

## Algorithm outline (Julia)

```julia
# Subtraction method for Thirring model
# 1. Define vector field u: ℝ^{4V} → ℝ^{4V} (neural network or analytic)
# 2. Subtracted Boltzmann factor:
#    e_sub(ϕ, u) = exp(-S(ϕ)) * (1 + dot(u(ϕ), ∇S(ϕ)) - divergence(u, ϕ))
# 3. Sample from e_sub using MCMC
# 4. Compute observables via eq. (15):
#    ⟨O⟩ = mean(O(ϕ) - exp(-S_u(ϕ))/exp(-S_v(ϕ)) * dot(u(ϕ), ∇O(ϕ)))
# 5. Train u by gradient descent on C[u] = -⟨σ⟩_u
#    Gradient: ∂C/∂λ = -⟨Re S⟩_Q (quenched, cheap to compute)
```

---

## Notation summary

| Symbol | Meaning |
|--------|---------|
| $g(\phi)$ | subtraction function, $\int g = 0$ |
| $v = e^{-S}u$ | scaled vector field; $u$ is $\mathcal{O}(1)$ |
| $e^{-S_u} = e^{-S}(1 + u\cdot\nabla S - \nabla\cdot u)$ | subtracted Boltzmann factor |
| $\langle\mathcal{O}\rangle_Q$ | quenched expectation value |
| $Z_{Q,u}$ | quenched partition function after subtraction |
| $g_\text{hd}$ | heavy-dense subtraction (large-$\mu$ analytic) |
| $C[u] = -\langle\sigma\rangle_u$ | cost function (average phase to maximise) |

---

## Open questions

- No rigorous argument for why estimator (15) avoids sign-to-noise
- Scaling to higher-dimensional lattices
- Handling Boltzmann factors with zeros (scaling trick breaks down)
- Optimal combination with contour deformation methods

---

## Links

- [[sign-problem]]
- [[subtraction-method]]
- [[contour-deformation]]
- [[average-phase]]
- [[thirring-model-action]]
- [[holomorphic-gradient-flow]]
