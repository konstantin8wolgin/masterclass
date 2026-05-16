# 0D phi4 Sign Problem — PIK Flow Simulation

## Context

Zero-dimensional phi4 theory with a complex action (real-time path integral). The action is
rotated from Euclidean (t=0) to Minkowski (t=1) along a one-parameter path. The sign problem is
addressed via the PIK (physics-informed kernel) flow from Ihssen, Kapust, Pawlowski 2026,
specifically Supplement S3.4. Observables are computed using the quenched reweighting formula.
Implementation in Julia.

Primary reference: `wiki/sources/ihssen-2026-solving-sign-problems-piks.md` (Supplement S3.4).
Supporting: `wiki/equations/wegner-equation.md`, `wiki/concepts/sign-problem.md`.

---

## Physics

### Action path

    S_t(phi) = m_t^2/2 * phi^2 + lambda_t/4 * phi^4

    m_t^2  = m^2  * exp(-i*pi*t/2)
    lambda_t = lambda * exp(-i*pi*t/2)

Parameters: m^2 = lambda = 1, t in [0, 1].

At t=0: real Euclidean action, no sign problem.
At t=1: S_1 = -i * S_E (Minkowski prefactor e^{i*S_M}), severe sign problem.

### Key identity

    dS_t/dt = (-i*pi/2) * S_t(phi)

Both couplings carry exp(-i*pi*t/2), so the t-derivative factors out -i*pi/2 exactly.

### Wegner equation

    d(phi_dot_t)/dphi - dS_t/dphi * phi_dot_t = dS_t/dt = (-i*pi/2) * S_t(phi)

This is a linear first-order ODE in phi_dot_t as a function of phi (with t as a parameter).

Closed-form solution via integrating factor exp(-S_hat_t), boundary condition phi_dot_t -> 0 at +-inf:

    phi_dot_t(phi) = (i*pi/2) * exp(S_hat_t(phi)) * integral_{phi}^{+inf} [S_hat_t(phi') - <S_hat_t>] * exp(-S_hat_t(phi')) dphi'

where S_hat_t is the field-dependent part of the action (no normalization constant), and:

    <S_hat_t> = integral S_hat_t(phi) exp(-S_hat_t(phi)) dphi / integral exp(-S_hat_t(phi)) dphi

The centering by <S_hat_t> comes from the normalization flow: d log Z_t/dt = (i*pi/2) <S_hat_t>.
Without it the formula is wrong — verified by free-theory check below.
The tail integral converges because Re[S_hat_t(phi)] -> +inf as |phi| -> inf for t in [0, 1).

### Free-theory validation (lambda=0)

Analytic result (Ihssen 2026 Eq. 71):

    phi_dot_t = (i*pi/4) * phi
    phi_t(phi_0) = phi_0 * exp(i*pi*t/4)

Used in script 01 to validate the numerical integral formula.

### Quenched distribution

    Q_t(phi) proportional to exp(-Re[S_t(phi)])

### Observable formula (reweighting)

    <O>_t = <O * exp(-i * Im[S_t])>_Q / <exp(-i * Im[S_t])>_Q

where both expectation values are computed by direct numerical quadrature.
Observables: O = phi^2 (two-point function) and O = |phi|.

---

## Architecture

Single Julia project. Five scripts with JLD2 checkpoints.

```
demovault/projects/phi4/
├── CLAUDE.md         (not a project file — handled separately)
├── Project.toml
├── data/
│   ├── action_check.jld2
│   ├── wegner.jld2
│   ├── flow.jld2
│   └── observables.jld2
├── 01_action_check.jl
├── 02_wegner.jl
├── 03_flow.jl
├── 04_observables.jl
└── 05_summary.jl
```

---

## Script responsibilities

### 01_action_check.jl

- Define S(phi, t), dS_dt(phi, t) = (-i*pi/2)*S(phi,t), dS_dphi(phi, t) as pure functions
- Verify key identity at a grid of (phi, t) pairs: dS_dt ≈ (-i*pi/2) * S
- For lambda=0: compute phi_dot_t from the integral formula at several phi values,
  compare to analytic i*pi/4 * phi; print pass/fail for each
- Save parameters and reference values to data/action_check.jld2

### 02_wegner.jl

- Load data/action_check.jld2
- phi_grid in [-5, 5], N_phi = 500; t_grid: 50 values in [0, 1]
- Compute phi_dot_t(phi) via QuadGK for all (t, phi) pairs
  Upper limit of tail integral: cutoff where exp(-Re[S_t(L)]) < 1e-12
- Save (phi_grid, t_grid, phi_dot) to data/wegner.jld2
- Plot Re[phi_dot_t] and Im[phi_dot_t] vs phi for t = 0, 0.25, 0.5, 0.75, 1

### 03_flow.jl

- Load data/wegner.jld2; build interpolant for phi_dot_t(phi) in phi and t
- 20 initial values phi_0 uniformly in [-3, 3]
- Integrate dphi/dt = phi_dot_t(phi_t) with fixed-step RK4, N_steps = 512
- phi_t is complex for t > 0
- Save complex trajectories (shape: N_init x N_t_steps) to data/flow.jld2
- Plot trajectories in the complex phi-plane (Re vs Im), colored by t value
- Mark the t=1 endpoints (PIKfold) prominently

### 04_observables.jl

- Independent — does not need flow or Wegner data
- t_grid: 100 values in [0, 1]
- For each t, compute via QuadGK (rtol = 1e-8):
    num_phi2   = integral phi^2   * exp(-i*Im[S]) * exp(-Re[S]) dphi
    num_absphi = integral |phi|   * exp(-i*Im[S]) * exp(-Re[S]) dphi
    denom      = integral           exp(-i*Im[S]) * exp(-Re[S]) dphi
    norm_Q     = integral           exp(-Re[S])                 dphi
    <phi^2>    = num_phi2   / denom
    <|phi|>    = num_absphi / denom
    avg_sign   = denom / norm_Q
- Save (t_grid, phi2, absphi, avg_sign) to data/observables.jld2
- Plot Re[<phi^2>], Re[<|phi|>] vs t
- Plot |avg_sign| vs t (sign problem severity)

### 05_summary.jl

- Load data/observables.jld2, data/flow.jld2, data/wegner.jld2
- Print formatted table at t = 0, 0.25, 0.5, 0.75, 1.0:
  columns: t | Re<phi^2> | Im<phi^2> | Re<|phi|> | |avg_sign|
- Assemble 2x2 combined figure:
    top-left:     Re/Im phi_dot_t vs phi (selected t values, from wegner.jld2)
    top-right:    PIKfold trajectories in complex plane (from flow.jld2)
    bottom-left:  Re[<phi^2>] and Re[<|phi|>] vs t
    bottom-right: |avg_sign| vs t

---

## Numerical parameters

- phi_grid: [-5, 5], N_phi = 500
- t_grid (Wegner): 50 points in [0, 1]
- t_grid (observables): 100 points in [0, 1]
- Flow: N_init = 20, N_steps = 512
- QuadGK rtol: 1e-8
- Tail cutoff: L where exp(-Re[S_t(L)]) < 1e-12

---

## Dependencies

- QuadGK.jl — adaptive 1D quadrature
- JLD2.jl — checkpoint serialization (preserves ComplexF64 arrays)
- Plots.jl — output figures
- Interpolations.jl — bilinear interpolant for phi_dot_t(phi, t) used in script 03
- (No DifferentialEquations.jl needed — manual RK4 in script 03)

---

## Validation checks

1. Script 01: free-theory phi_dot_t from integral vs analytic i*pi/4 * phi
2. Script 04: at t=0, avg_sign = 1 (no sign problem)
3. Script 04: at t=0, <phi^2> and <|phi|> match known Euclidean integrals
   (for m^2=lambda=1: compute reference with QuadGK using real S_0)
4. Script 03: trajectories are Z2-symmetric (phi_0 and -phi_0 give mirror trajectories)
   because S_t is even in phi
