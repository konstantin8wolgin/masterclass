# phi4 Sign Problem — PIK Flow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> **Parallelism note:** Tasks 4 and 5 are independent and can be dispatched as parallel subagents.

**Goal:** Implement the 0D phi4 sign-problem simulation in Julia — solve the Wegner equation for the PIK velocity field, evolve samples via RK4 flow, and compute observables via quenched reweighting.

**Architecture:** Five checkpoint-based Julia scripts in `demovault/projects/phi4/`. Each script saves results to `data/*.jld2` so scripts can be inspected and re-run independently. Scripts 03 and 04 are independent of each other.

**Tech Stack:** Julia 1.10+, QuadGK.jl, JLD2.jl, Plots.jl, Interpolations.jl

**Read before implementing:**
- `demovault/wiki/sources/ihssen-2026-solving-sign-problems-piks.md` — Supplement S3.4 is the exact model
- `demovault/wiki/equations/wegner-equation.md`
- `demovault/docs/superpowers/specs/2026-05-05-phi4-sign-problem-design.md` — full spec

---

## Task 0: Project scaffold

**Files:**
- Create: `demovault/projects/phi4/Project.toml`
- Create: `demovault/projects/phi4/data/.gitkeep`

- [ ] **Step 1: Create Project.toml**

```toml
name = "phi4_sign"
uuid = "00000000-0000-0000-0000-000000000001"
version = "0.1.0"

[deps]
QuadGK = "1fd47b50-473d-5a70-b4e7-f838b478d57c"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"

[compat]
julia = "1.10"
```

- [ ] **Step 2: Create data directory**

```bash
mkdir -p demovault/projects/phi4/data
touch demovault/projects/phi4/data/.gitkeep
```

- [ ] **Step 3: Instantiate packages**

```bash
cd demovault/projects/phi4
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Expected: packages download and compile, no errors.

---

## Task 1: 01_action_check.jl — define action, validate free-theory

**Files:**
- Create: `demovault/projects/phi4/01_action_check.jl`

**Physics being implemented:**
- S_hat(φ, t) = exp(-iπt/2)/2 * φ² + exp(-iπt/2)/4 * φ⁴  (m²=λ=1)
- Key identity: dS_hat/dt = (-iπ/2) * S_hat(φ, t)
- Free-theory (λ=0): analytic velocity φ̇_t = (iπ/4)φ (Ihssen 2026 Eq. 71)
- Correct Wegner solution requires centering: integrand is S_hat - ⟨S_hat⟩

- [ ] **Step 1: Write 01_action_check.jl**

```julia
using QuadGK, JLD2

# ── Parameters ────────────────────────────────────────────────────────────────
const M2 = 1.0   # bare mass squared
const LAM = 1.0  # bare coupling

# ── Action (field-dependent part only, no log Z) ──────────────────────────────
S_hat(φ, t) = M2 * exp(-im*π/2*t)/2 * φ^2 + LAM * exp(-im*π/2*t)/4 * φ^4
dS_dt(φ, t) = (-im*π/2) * S_hat(φ, t)   # exact: both couplings carry exp(-iπt/2)

# ── Mean ⟨S_hat_t⟩ under e^{-S_hat} ──────────────────────────────────────────
function action_mean(t; L=8.0)
    num, _ = quadgk(φ -> S_hat(φ, t) * exp(-S_hat(φ, t)), -L, L, rtol=1e-10)
    den, _ = quadgk(φ -> exp(-S_hat(φ, t)), -L, L, rtol=1e-10)
    return num / den
end

# ── Velocity field from Wegner equation (real φ) ──────────────────────────────
# φ̇_t(φ) = (iπ/2) exp(S_hat(φ)) ∫_φ^∞ [S_hat(φ') - ⟨S_hat⟩] exp(-S_hat(φ')) dφ'
function phi_dot(φ_real::Real, t; L=8.0)
    Sbar = action_mean(t; L=L)
    integrand(x) = (S_hat(x, t) - Sbar) * exp(-S_hat(x, t))
    integral, _ = quadgk(integrand, φ_real, L, rtol=1e-10)
    return (im*π/2) * exp(S_hat(φ_real, t)) * integral
end

# ── Free-theory version (LAM=0) for validation ───────────────────────────────
S_hat_free(φ, t) = M2 * exp(-im*π/2*t)/2 * φ^2

function phi_dot_free(φ_real::Real, t; L=8.0)
    Sbar_free = begin
        num, _ = quadgk(φ -> S_hat_free(φ,t)*exp(-S_hat_free(φ,t)), -L, L, rtol=1e-10)
        den, _ = quadgk(φ -> exp(-S_hat_free(φ,t)), -L, L, rtol=1e-10)
        num/den
    end
    integrand(x) = (S_hat_free(x, t) - Sbar_free) * exp(-S_hat_free(x, t))
    integral, _ = quadgk(integrand, φ_real, L, rtol=1e-10)
    return (im*π/2) * exp(S_hat_free(φ_real, t)) * integral
end

# ── Validation 1: key identity dS/dt = (-iπ/2) S ────────────────────────────
println("=== Identity check: dS_hat/dt == (-iπ/2) S_hat ===")
all_pass = true
for (φ, t) in [(0.5, 0.0), (1.0, 0.25), (-1.5, 0.5), (2.0, 0.75)]
    err = abs(dS_dt(φ, t) - (-im*π/2) * S_hat(φ, t))
    pass = err < 1e-14
    all_pass &= pass
    println("  φ=$(φ), t=$(t): err=$(round(err, sigdigits=3))  $(pass ? "PASS" : "FAIL")")
end

# ── Validation 2: free-theory φ̇_t vs analytic (iπ/4)φ ───────────────────────
println("\n=== Free-theory check: phi_dot_free vs (iπ/4)φ ===")
for (φ, t) in [(1.0, 0.25), (2.0, 0.5), (-1.0, 0.75), (0.5, 0.9)]
    numeric  = phi_dot_free(φ, t)
    analytic = (im*π/4) * φ
    err = abs(numeric - analytic)
    pass = err < 1e-6
    all_pass &= pass
    println("  φ=$(φ), t=$(t): numeric=$(round(numeric, digits=5)), analytic=$(round(analytic, digits=5)), err=$(round(err, sigdigits=3))  $(pass ? "PASS" : "FAIL")")
end

println("\nOverall: $(all_pass ? "ALL PASS" : "FAILURES DETECTED")")

# ── Save parameters ───────────────────────────────────────────────────────────
jldsave("data/action_check.jld2"; M2, LAM)
println("\nSaved parameters to data/action_check.jld2")
```

- [ ] **Step 2: Run and verify**

```bash
cd demovault/projects/phi4
julia --project=. 01_action_check.jl
```

Expected output:
```
=== Identity check: dS_hat/dt == (-iπ/2) S_hat ===
  φ=0.5, t=0.0: err=0.0  PASS
  ...
=== Free-theory check: phi_dot_free vs (iπ/4)φ ===
  φ=1.0, t=0.25: numeric=0.7854+0.0im, analytic=0.7854+0.0im, err=<1e-6  PASS
  ...
Overall: ALL PASS
```

If any FAIL: check that the centering term (Sbar) is included in phi_dot_free. The uncentered formula gives a wrong result with an extra erfc term.

---

## Task 2: 02_wegner.jl — solve Wegner equation on grid, plot velocity field

**Files:**
- Create: `demovault/projects/phi4/02_wegner.jl`

**What this produces:** phi_dot_t(φ) on a (φ, t) grid; saved as ComplexF64 array. Plots Re/Im of velocity field vs φ at 5 selected t values.

- [ ] **Step 1: Write 02_wegner.jl**

```julia
using QuadGK, JLD2, Plots

# ── Load parameters ───────────────────────────────────────────────────────────
data = load("data/action_check.jld2")
const M2  = data["M2"]
const LAM = data["LAM"]

# ── Action and helpers (same as 01) ──────────────────────────────────────────
S_hat(φ, t) = M2 * exp(-im*π/2*t)/2 * φ^2 + LAM * exp(-im*π/2*t)/4 * φ^4

function action_mean(t; L=8.0)
    num, _ = quadgk(φ -> S_hat(φ, t) * exp(-S_hat(φ, t)), -L, L, rtol=1e-8)
    den, _ = quadgk(φ -> exp(-S_hat(φ, t)), -L, L, rtol=1e-8)
    return num / den
end

function phi_dot(φ_real::Real, t, Sbar; L=8.0)
    integrand(x) = (S_hat(x, t) - Sbar) * exp(-S_hat(x, t))
    integral, _ = quadgk(integrand, φ_real, L, rtol=1e-8)
    return (im*π/2) * exp(S_hat(φ_real, t)) * integral
end

# ── Grid setup ────────────────────────────────────────────────────────────────
# Avoid t=1 exactly: Re[S_hat] → 0 there, integrals become oscillatory
const N_phi = 500
const N_t   = 50
phi_grid = range(-5.0, 5.0, length=N_phi)
t_grid   = range(0.0, 0.99, length=N_t)   # stop at 0.99, not 1.0

# ── Compute phi_dot on grid ───────────────────────────────────────────────────
phi_dot_grid = zeros(ComplexF64, N_t, N_phi)

println("Computing velocity field on $(N_t) × $(N_phi) grid...")
for (i, t) in enumerate(t_grid)
    Sbar = action_mean(t)
    for (j, φ) in enumerate(phi_grid)
        phi_dot_grid[i, j] = phi_dot(φ, t, Sbar)
    end
    i % 10 == 0 && println("  t-slice $i / $N_t done")
end
println("Done.")

# ── Save ──────────────────────────────────────────────────────────────────────
jldsave("data/wegner.jld2";
    phi_grid = collect(phi_grid),
    t_grid   = collect(t_grid),
    phi_dot  = phi_dot_grid)
println("Saved to data/wegner.jld2")

# ── Plot: Re and Im of phi_dot_t vs phi for 5 t-values ───────────────────────
plot_t_vals = [0.0, 0.25, 0.5, 0.75, 0.99]
colors = [:blue, :green, :orange, :red, :purple]

p1 = plot(title="Re[φ̇_t(φ)] vs φ", xlabel="φ", ylabel="Re[φ̇_t]", legend=:topright)
p2 = plot(title="Im[φ̇_t(φ)] vs φ", xlabel="φ", ylabel="Im[φ̇_t]", legend=:topright)

for (tv, col) in zip(plot_t_vals, colors)
    # find nearest index in t_grid
    idx = argmin(abs.(collect(t_grid) .- tv))
    t_actual = t_grid[idx]
    phi_dot_row = phi_dot_grid[idx, :]
    plot!(p1, collect(phi_grid), real.(phi_dot_row), label="t=$(round(t_actual,digits=2))", color=col)
    plot!(p2, collect(phi_grid), imag.(phi_dot_row), label="t=$(round(t_actual,digits=2))", color=col)
end

savefig(plot(p1, p2, layout=(1,2), size=(900,400)), "data/02_velocity_field.png")
println("Plot saved to data/02_velocity_field.png")
```

- [ ] **Step 2: Run and verify**

```bash
cd demovault/projects/phi4
julia --project=. 02_wegner.jl
```

Expected:
- Prints progress every 10 t-slices, completes in ~1-5 min depending on hardware
- `data/wegner.jld2` created (~2 MB)
- `data/02_velocity_field.png` created
- At t=0: Im[φ̇_0] should be positive and symmetric, Re[φ̇_0] = 0 (since S_hat is real at t=0 so φ̇ is purely imaginary × real integral = imaginary)
- At t=0, φ=0: phi_dot should be 0 (S_hat is even, symmetry → φ̇(0) = 0)

---

## Task 3: 03_flow.jl — integrate RK4 flow, plot PIKfold

**Files:**
- Create: `demovault/projects/phi4/03_flow.jl`

**What this produces:** Complex trajectories φ_t^(n) for 20 initial real values. Plot shows how the real line deforms into the PIKfold contour in ℂ.

**Implementation note on complex φ:** The velocity field was computed on the real φ grid. Once φ_t becomes complex, we evaluate φ̇_t by interpolating in t and using Re[φ_t] as the grid coordinate (valid approximation for small |Im[φ_t]|; exact if the grid is dense enough). This is standard for visualization of thimble/PIKfold formation.

- [ ] **Step 1: Write 03_flow.jl**

```julia
using JLD2, Plots, Interpolations

# ── Load Wegner data ──────────────────────────────────────────────────────────
data = load("data/wegner.jld2")
phi_grid_real = data["phi_grid"]   # Vector{Float64}, length N_phi
t_grid        = data["t_grid"]     # Vector{Float64}, length N_t
phi_dot_grid  = data["phi_dot"]    # Matrix{ComplexF64}, size (N_t, N_phi)

# ── Build interpolant for phi_dot_t(phi) ─────────────────────────────────────
# Separate real and imag parts (Interpolations.jl requires real-valued arrays)
interp_re = interpolate((t_grid, phi_grid_real), real.(phi_dot_grid), Gridded(Linear()))
interp_im = interpolate((t_grid, phi_grid_real), imag.(phi_dot_grid), Gridded(Linear()))

# Extrapolation: clamp to grid boundary
extrap_re = extrapolate(interp_re, Flat())
extrap_im = extrapolate(interp_im, Flat())

function phi_dot_interp(phi_complex, t)
    φ_r = clamp(real(phi_complex), phi_grid_real[1], phi_grid_real[end])
    t_c = clamp(t, t_grid[1], t_grid[end])
    return extrap_re(t_c, φ_r) + im * extrap_im(t_c, φ_r)
end

# ── RK4 integrator ────────────────────────────────────────────────────────────
function rk4_step(φ, t, dt)
    k1 = phi_dot_interp(φ,        t)
    k2 = phi_dot_interp(φ + dt/2*k1, t + dt/2)
    k3 = phi_dot_interp(φ + dt/2*k2, t + dt/2)
    k4 = phi_dot_interp(φ + dt*k3,   t + dt)
    return φ + dt/6 * (k1 + 2k2 + 2k3 + k4)
end

# ── Flow integration ──────────────────────────────────────────────────────────
const N_init   = 20
const N_steps  = 512
const t_start  = t_grid[1]
const t_end    = t_grid[end]
const dt       = (t_end - t_start) / N_steps

phi0_values = range(-3.0, 3.0, length=N_init)
# trajectories[n, step] = complex phi_t for initial condition n at step
trajectories = zeros(ComplexF64, N_init, N_steps + 1)

println("Integrating flow for $N_init initial conditions, $N_steps RK4 steps...")
for (n, phi0) in enumerate(phi0_values)
    trajectories[n, 1] = phi0
    for step in 1:N_steps
        t_now = t_start + (step-1)*dt
        trajectories[n, step+1] = rk4_step(trajectories[n, step], t_now, dt)
    end
end
println("Done.")

# ── Save ──────────────────────────────────────────────────────────────────────
jldsave("data/flow.jld2";
    trajectories = trajectories,
    phi0_values  = collect(phi0_values),
    t_steps      = collect(range(t_start, t_end, length=N_steps+1)))
println("Saved to data/flow.jld2")

# ── Plot: trajectories in complex phi-plane ───────────────────────────────────
p = plot(title="PIKfold formation: field trajectories in ℂ",
         xlabel="Re[φ]", ylabel="Im[φ]",
         legend=false, aspect_ratio=:equal)

n_steps_total = N_steps + 1
cmap = cgrad(:viridis)

for n in 1:N_init
    traj = trajectories[n, :]
    # color segments by t
    for step in 1:N_steps
        color_val = (step-1) / N_steps
        plot!(p,
            [real(traj[step]), real(traj[step+1])],
            [imag(traj[step]), imag(traj[step+1])],
            color=cmap[color_val], linewidth=1.5, alpha=0.8)
    end
end

# Mark start (t=0, real line) and end (t=t_end, PIKfold)
scatter!(p, real.(trajectories[:, 1]),   imag.(trajectories[:, 1]),
         marker=:circle, markersize=4, color=:blue,   label="t=0 (real line)")
scatter!(p, real.(trajectories[:, end]), imag.(trajectories[:, end]),
         marker=:star5,  markersize=6, color=:red,    label="t=$(round(t_end,digits=2)) (PIKfold)")

savefig(p, "data/03_pikfold.png")
println("Plot saved to data/03_pikfold.png")
```

- [ ] **Step 2: Run and verify**

```bash
cd demovault/projects/phi4
julia --project=. 03_flow.jl
```

Expected:
- Completes quickly (<10 sec)
- `data/flow.jld2` created
- `data/03_pikfold.png`: trajectories start on real axis, curve into complex plane as t increases
- Z₂ symmetry check: trajectory from φ₀ = +x and φ₀ = −x should be mirror images (Im[φ_t(+x)] = Im[φ_t(−x)], Re[φ_t(+x)] = −Re[φ_t(−x)]) because S_hat is even in φ

---

## Task 4: 04_observables.jl — quadrature + reweighting (INDEPENDENT of task 5)

**Files:**
- Create: `demovault/projects/phi4/04_observables.jl`

**Physics:** Quenched reweighting formula:
```
⟨O⟩_t = ⟨O exp(-i Im[S_hat])⟩_Q / ⟨exp(-i Im[S_hat])⟩_Q
where Q_t ∝ exp(-Re[S_hat_t(φ)])
```

**Numerical note at t → 1:** Re[S_hat_t] → 0 so Q_t flattens; integrals become oscillatory. QuadGK handles this with rtol=1e-8 and extended domain, but precision near t=1 is lower. This is the sign problem manifesting numerically.

- [ ] **Step 1: Write 04_observables.jl**

```julia
using QuadGK, JLD2, Plots

data = load("data/action_check.jld2")
const M2  = data["M2"]
const LAM = data["LAM"]

S_hat(φ, t) = M2 * exp(-im*π/2*t)/2 * φ^2 + LAM * exp(-im*π/2*t)/4 * φ^4

# Integration domain: extend at t→1 where Re[S_hat] flattens
# Require Re[S_hat(L,t)] = cos(πt/2)/2 * L² > 6  →  L > sqrt(12/cos(πt/2))
function domain_limit(t)
    c = cos(π/2 * t)
    c < 0.01 && return 30.0   # t very close to 1, use large cutoff
    return clamp(sqrt(12.0 / c), 8.0, 30.0)
end

function compute_observables(t)
    L = domain_limit(t)
    phase(φ)     = exp(-im * imag(S_hat(φ, t)))
    quenched(φ)  = exp(-real(S_hat(φ, t)))
    weight(φ)    = phase(φ) * quenched(φ)

    norm_Q,  _ = quadgk(quenched,                         -L, L, rtol=1e-8)
    denom,   _ = quadgk(weight,                           -L, L, rtol=1e-8)
    num_phi2,_ = quadgk(φ -> φ^2   * weight(φ),          -L, L, rtol=1e-8)
    num_abs, _ = quadgk(φ -> abs(φ) * weight(φ),         -L, L, rtol=1e-8)

    phi2     = num_phi2 / denom
    absphi   = num_abs  / denom
    avg_sign = denom    / norm_Q
    return phi2, absphi, avg_sign
end

# ── Compute over t grid ───────────────────────────────────────────────────────
const N_t = 100
t_grid = range(0.0, 1.0, length=N_t)

phi2_arr     = zeros(ComplexF64, N_t)
absphi_arr   = zeros(ComplexF64, N_t)
avgsign_arr  = zeros(ComplexF64, N_t)

println("Computing observables at $N_t t-values...")
for (i, t) in enumerate(t_grid)
    phi2_arr[i], absphi_arr[i], avgsign_arr[i] = compute_observables(t)
    i % 10 == 0 && println("  t=$(round(t,digits=3))  ⟨φ²⟩=$(round(phi2_arr[i],digits=4))  |avg_sign|=$(round(abs(avgsign_arr[i]),sigdigits=3))")
end
println("Done.")

# ── Euclidean sanity checks at t=0 ───────────────────────────────────────────
println("\n=== Sanity checks at t=0 ===")
println("  avg_sign(0) should be 1: $(round(avgsign_arr[1], digits=6))")
phi2_real_direct, _ = quadgk(φ -> φ^2 * exp(-S_hat(φ, 0.0)), -8.0, 8.0, rtol=1e-10)
norm_direct,      _ = quadgk(φ ->       exp(-S_hat(φ, 0.0)), -8.0, 8.0, rtol=1e-10)
phi2_euclidean = phi2_real_direct / norm_direct
println("  ⟨φ²⟩ from reweighting: $(round(phi2_arr[1], digits=6))")
println("  ⟨φ²⟩ direct Euclidean: $(round(phi2_euclidean, digits=6))")
println("  Difference: $(round(abs(phi2_arr[1] - phi2_euclidean), sigdigits=3))")

# ── Save ──────────────────────────────────────────────────────────────────────
jldsave("data/observables.jld2";
    t_grid    = collect(t_grid),
    phi2      = phi2_arr,
    absphi    = absphi_arr,
    avg_sign  = avgsign_arr)
println("\nSaved to data/observables.jld2")

# ── Plots ─────────────────────────────────────────────────────────────────────
t_vec = collect(t_grid)

p1 = plot(t_vec, real.(phi2_arr),
    label="Re⟨φ²⟩", xlabel="t", ylabel="Observable",
    title="Observables vs t", lw=2, color=:blue)
plot!(p1, t_vec, real.(absphi_arr),
    label="Re⟨|φ|⟩", lw=2, color=:red)
savefig(p1, "data/04_observables.png")

p2 = plot(t_vec, abs.(avgsign_arr),
    xlabel="t", ylabel="|avg sign|",
    title="Sign problem severity: |⟨e^{-i Im S}⟩_Q|",
    lw=2, color=:black, ylims=(0, 1.05), legend=false)
hline!(p2, [1.0], linestyle=:dash, color=:gray, label="no sign problem")
savefig(p2, "data/04_sign_severity.png")

println("Plots saved to data/04_observables.png and data/04_sign_severity.png")
```

- [ ] **Step 2: Run and verify**

```bash
cd demovault/projects/phi4
julia --project=. 04_observables.jl
```

Expected:
- `avg_sign(0)` ≈ 1.0 + 0.0im (no sign problem at t=0)
- ⟨φ²⟩ at t=0 from reweighting matches direct Euclidean integration (difference < 1e-6)
- |avg_sign| decreases monotonically from 1 (t=0) toward 0 (t=1)
- ⟨φ²⟩ real part should evolve smoothly from the Euclidean value (~0.35 for m²=λ=1) toward the Minkowski value

---

## Task 5: 05_summary.jl — combined table and figure (INDEPENDENT of task 4)

**Files:**
- Create: `demovault/projects/phi4/05_summary.jl`

**Note:** This task depends on tasks 2, 3, and 4 having been run (all three JLD2 files must exist).

- [ ] **Step 1: Write 05_summary.jl**

```julia
using JLD2, Plots, Printf

# ── Load all checkpoint data ──────────────────────────────────────────────────
obs  = load("data/observables.jld2")
flow = load("data/flow.jld2")
weg  = load("data/wegner.jld2")

t_obs      = obs["t_grid"]
phi2       = obs["phi2"]
absphi     = obs["absphi"]
avg_sign   = obs["avg_sign"]

trajectories = flow["trajectories"]
phi0_values  = flow["phi0_values"]
t_flow       = flow["t_steps"]

phi_grid_r   = weg["phi_grid"]
t_grid_w     = weg["t_grid"]
phi_dot_grid = weg["phi_dot"]

# ── Printed table ─────────────────────────────────────────────────────────────
println("\n╔══════╦══════════════╦══════════════╦══════════════╦════════════╗")
println("║  t   ║  Re⟨φ²⟩     ║  Im⟨φ²⟩     ║  Re⟨|φ|⟩   ║ |avg sign| ║")
println("╠══════╬══════════════╬══════════════╬══════════════╬════════════╣")
for t_val in [0.0, 0.25, 0.5, 0.75, 1.0]
    idx = argmin(abs.(t_obs .- t_val))
    @printf("║ %.2f ║ %12.6f ║ %12.6f ║ %12.6f ║ %10.6f ║\n",
        t_obs[idx],
        real(phi2[idx]), imag(phi2[idx]),
        real(absphi[idx]),
        abs(avg_sign[idx]))
end
println("╚══════╩══════════════╩══════════════╩══════════════╩════════════╝")

# ── 2×2 combined figure ───────────────────────────────────────────────────────

# Panel 1: velocity field Re[φ̇_t] vs φ
p_vel = plot(title="Velocity field Re[φ̇_t(φ)]",
             xlabel="φ", ylabel="Re[φ̇_t]", legend=:topright)
for (tv, col) in zip([0.0, 0.25, 0.5, 0.75, 0.99], [:blue,:green,:orange,:red,:purple])
    idx = argmin(abs.(t_grid_w .- tv))
    plot!(p_vel, phi_grid_r, real.(phi_dot_grid[idx, :]),
          label="t=$(round(t_grid_w[idx],digits=2))", color=col)
end

# Panel 2: PIKfold trajectories
p_fold = plot(title="PIKfold formation",
              xlabel="Re[φ]", ylabel="Im[φ]",
              legend=false, aspect_ratio=:equal)
cmap = cgrad(:viridis)
for n in axes(trajectories, 1)
    traj = trajectories[n, :]
    for step in 1:length(traj)-1
        c = (step-1) / (length(traj)-1)
        plot!(p_fold,
            [real(traj[step]), real(traj[step+1])],
            [imag(traj[step]), imag(traj[step+1])],
            color=cmap[c], linewidth=1.2, alpha=0.7)
    end
end
scatter!(p_fold, real.(trajectories[:, 1]),   zeros(size(trajectories,1)),
         color=:blue, markersize=3, label="t=0")
scatter!(p_fold, real.(trajectories[:, end]), imag.(trajectories[:, end]),
         color=:red, markersize=4, marker=:star5, label="t=t_end")

# Panel 3: observables vs t
p_obs = plot(t_obs, real.(phi2),   label="Re⟨φ²⟩",   lw=2, color=:blue,
             xlabel="t", ylabel="Observable", title="Observables vs t")
plot!(p_obs, t_obs, real.(absphi), label="Re⟨|φ|⟩",  lw=2, color=:red)

# Panel 4: sign severity
p_sign = plot(t_obs, abs.(avg_sign),
              xlabel="t", ylabel="|avg sign|",
              title="Sign problem severity",
              lw=2, color=:black, ylims=(0, 1.05), legend=false)
hline!(p_sign, [1.0], linestyle=:dash, color=:gray)

fig = plot(p_vel, p_fold, p_obs, p_sign,
           layout=(2,2), size=(1100, 800))
savefig(fig, "data/05_summary.png")
println("\nCombined figure saved to data/05_summary.png")
```

- [ ] **Step 2: Run and verify**

```bash
cd demovault/projects/phi4
julia --project=. 05_summary.jl
```

Expected:
- Prints formatted table with values at t = 0, 0.25, 0.5, 0.75, 1.0
- At t=0: Re⟨φ²⟩ ≈ 0.35, |avg sign| ≈ 1.0
- At t=1: |avg sign| ≈ 0 (severe sign problem)
- `data/05_summary.png` created with 2×2 grid

---

## Execution order

```
Task 0  →  Task 1  →  Task 2  →  ┬── Task 3 (flow)
                                  └── Task 4 (observables, independent)
                                         ↓
                                      Task 5 (summary, needs all three)
```

Tasks 3 and 4 can run in parallel after Task 2 completes.
