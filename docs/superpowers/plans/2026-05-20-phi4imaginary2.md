# phi4imaginary2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a new parallel Julia project at `demovault/projects/phi4imaginary2` that implements a clean Tikhonov-regularised, real-axis phi4 imaginary-time rotation benchmark without modifying `demovault/projects/phi4imaginary`.

**Architecture:** Build `phi4imaginary2` as a small library plus numbered scripts. The library owns the action, Tikhonov regulator, stable symmetry-preserving real-axis Wegner velocity, RK4 flow, direct quadrature references, and validation helpers. Scripts only orchestrate data generation and plots. The transported-flow estimator is the primary result; direct Tikhonov real-axis quadrature is a reference benchmark.

**Tech Stack:** Julia 1.10, `QuadGK`, `JLD2`, `Plots`, `Interpolations`, stdlib `Test`, `Random`, `Statistics`.

---

## Grill-Me Decisions Already Resolved

1. `phi4imaginary2` is Tikhonov-only. Do not preserve bare/unregularised runtime modes.
2. Velocity is real-axis only, evaluated for real initial coordinates. Do not claim a full complex-plane kernel.
3. Observables use the transported-flow estimator as primary output and direct Tikhonov real-axis quadrature as reference.
4. Existing `demovault/projects/phi4imaginary` must not be modified.

## File Structure

- Create: `demovault/projects/phi4imaginary2/Project.toml`  
  Project dependencies.
- Create: `demovault/projects/phi4imaginary2/src/Phi4Imaginary2.jl`  
  Core implementation: parameters, action, derivatives, stable velocity, flow, quadrature, symmetry checks.
- Create: `demovault/projects/phi4imaginary2/test/runtests.jl`  
  Enforced tests for invariants and smoke-level numerical behavior.
- Create: `demovault/projects/phi4imaginary2/scripts/01_velocity.jl`  
  Compute and save velocity grid plus symmetry diagnostics.
- Create: `demovault/projects/phi4imaginary2/scripts/02_flow.jl`  
  Integrate real initial points using the real-axis velocity map.
- Create: `demovault/projects/phi4imaginary2/scripts/03_observables.jl`  
  Compute transported-flow observables and quadrature references.
- Create: `demovault/projects/phi4imaginary2/scripts/04_summary.jl`  
  Generate summary table and figure.
- Create: `demovault/projects/phi4imaginary2/README.md`  
  State scope, regulator, limitations, run order, and expected validations.
- Create: `demovault/projects/phi4imaginary2/data/.gitkeep`  
  Keep output directory in git.

## Non-Negotiable Corrections

- Use `S_reg(phi,t) = S_bare(phi,t) + eps * phi^2 / 2`.
- Use `dS_reg_dt(phi,t) = dS_bare_dt(phi,t) = (-im*pi/2) * S_bare(phi,t)`.
- Do not multiply by `(-im*pi/2) * S_reg`; the Tikhonov term has no `t` derivative.
- Compute the velocity for `phi >= 0` with the tail integral and define negative values by odd symmetry.
- Include an exact `phi = 0` grid point by using an odd number of phi samples.
- Never label a plot as PIKfold or Lefschetz thimble unless the PIKfold ODE is actually solved. In this project, use labels like "transported real-axis samples".
- Do not use `abs(phi)` as a primary contour observable. Primary observables are holomorphic powers `phi^2` and `phi^4`.
- Every script must use `@__DIR__`-derived paths so it can run from any working directory.
- Every JLD2 output must save metadata: `m2`, `lam`, `eps`, `n_phi`, `n_t`, and the producing script name.

---

### Task 1: Create Project Skeleton

**Files:**
- Create: `demovault/projects/phi4imaginary2/Project.toml`
- Create: `demovault/projects/phi4imaginary2/src/Phi4Imaginary2.jl`
- Create: `demovault/projects/phi4imaginary2/test/runtests.jl`
- Create: `demovault/projects/phi4imaginary2/scripts/`
- Create: `demovault/projects/phi4imaginary2/data/.gitkeep`

- [ ] **Step 1: Create directories**

Run:

```bash
mkdir -p demovault/projects/phi4imaginary2/src demovault/projects/phi4imaginary2/test demovault/projects/phi4imaginary2/scripts demovault/projects/phi4imaginary2/data
touch demovault/projects/phi4imaginary2/data/.gitkeep
```

Expected: directories exist; `demovault/projects/phi4imaginary` is unchanged.

- [ ] **Step 2: Write `Project.toml`**

Create `demovault/projects/phi4imaginary2/Project.toml`:

```toml
name = "Phi4Imaginary2"
uuid = "11111111-2222-3333-4444-555555555555"
version = "0.1.0"

[deps]
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"

[compat]
Interpolations = "0.16"
JLD2 = "0.6"
Plots = "1"
QuadGK = "2"
julia = "1.10"
```

- [ ] **Step 3: Write an initial failing test file**

Create `demovault/projects/phi4imaginary2/test/runtests.jl`:

```julia
using Test
using Phi4Imaginary2

@testset "Phi4Imaginary2 skeleton" begin
    params = Params()
    @test params.m2 == 1.0
    @test params.lam == 1.0
    @test params.eps == 1e-3
end
```

- [ ] **Step 4: Run tests and verify failure**

Run:

```bash
julia --project=demovault/projects/phi4imaginary2 -e 'using Pkg; Pkg.test()'
```

Expected: FAIL because `src/Phi4Imaginary2.jl` does not exist or does not define `Params`.

- [ ] **Step 5: Create minimal module**

Create `demovault/projects/phi4imaginary2/src/Phi4Imaginary2.jl`:

```julia
module Phi4Imaginary2

export Params

Base.@kwdef struct Params
    m2::Float64 = 1.0
    lam::Float64 = 1.0
    eps::Float64 = 1e-3
end

end
```

- [ ] **Step 6: Run tests and verify pass**

Run:

```bash
julia --project=demovault/projects/phi4imaginary2 -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add demovault/projects/phi4imaginary2
git commit -m "feat: scaffold phi4imaginary2 project"
```

---

### Task 2: Implement Action and Tikhonov Derivatives

**Files:**
- Modify: `demovault/projects/phi4imaginary2/src/Phi4Imaginary2.jl`
- Modify: `demovault/projects/phi4imaginary2/test/runtests.jl`

- [ ] **Step 1: Add failing action tests**

Append to `test/runtests.jl`:

```julia
@testset "Tikhonov action" begin
    p = Params()
    for t in (0.0, 0.25, 0.5, 1.0), phi in (-2.0, -0.5, 0.0, 0.5, 2.0)
        @test S_bare(phi, t, p) ≈ S_bare(-phi, t, p)
        @test S_reg(phi, t, p) ≈ S_reg(-phi, t, p)
        @test dSreg_dt(phi, t, p) ≈ (-im*pi/2) * S_bare(phi, t, p)
    end
    @test S_reg(2.0, 0.7, p) - S_bare(2.0, 0.7, p) ≈ p.eps * 2.0^2 / 2
end
```

- [ ] **Step 2: Run tests and verify failure**

Run:

```bash
julia --project=demovault/projects/phi4imaginary2 -e 'using Pkg; Pkg.test()'
```

Expected: FAIL with `S_bare` / `S_reg` / `dSreg_dt` undefined.

- [ ] **Step 3: Implement action functions**

Replace `src/Phi4Imaginary2.jl` with:

```julia
module Phi4Imaginary2

export Params, S_bare, S_reg, dSreg_dt, dSreg_dphi, integration_bound

Base.@kwdef struct Params
    m2::Float64 = 1.0
    lam::Float64 = 1.0
    eps::Float64 = 1e-3
end

phase(t) = exp(-im * pi * t / 2)

function S_bare(phi, t, p::Params=Params())
    return phase(t) * (p.m2 * phi^2 / 2 + p.lam * phi^4 / 4)
end

function S_reg(phi, t, p::Params=Params())
    return S_bare(phi, t, p) + p.eps * phi^2 / 2
end

function dSreg_dt(phi, t, p::Params=Params())
    return (-im * pi / 2) * S_bare(phi, t, p)
end

function dSreg_dphi(phi, t, p::Params=Params())
    return phase(t) * (p.m2 * phi + p.lam * phi^3) + p.eps * phi
end

function integration_bound(t, p::Params=Params())
    c = max(cos(pi * t / 2), 0.0) + p.eps
    return clamp(sqrt(56.0 / c), 8.0, 80.0)
end

end
```

- [ ] **Step 4: Run tests and verify pass**

Run:

```bash
julia --project=demovault/projects/phi4imaginary2 -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add demovault/projects/phi4imaginary2/src/Phi4Imaginary2.jl demovault/projects/phi4imaginary2/test/runtests.jl
git commit -m "feat: add tikhonov action"
```

---

### Task 3: Implement Stable Odd Real-Axis Wegner Velocity

**Files:**
- Modify: `demovault/projects/phi4imaginary2/src/Phi4Imaginary2.jl`
- Modify: `demovault/projects/phi4imaginary2/test/runtests.jl`

- [ ] **Step 1: Add failing velocity tests**

Append to `test/runtests.jl`:

```julia
@testset "real-axis velocity invariants" begin
    p = Params()
    for t in (0.0, 0.25, 0.5, 0.75, 1.0)
        @test abs(phi_dot(0.0, t, p)) < 1e-10
        for phi in (0.25, 0.75, 1.5, 2.5, 3.5)
            vp = phi_dot(phi, t, p)
            vm = phi_dot(-phi, t, p)
            @test isfinite(real(vp))
            @test isfinite(imag(vp))
            @test abs(vp + vm) < 1e-7 * max(1.0, abs(vp))
        end
    end
end
```

- [ ] **Step 2: Run tests and verify failure**

Run:

```bash
julia --project=demovault/projects/phi4imaginary2 -e 'using Pkg; Pkg.test()'
```

Expected: FAIL with `phi_dot` undefined.

- [ ] **Step 3: Implement centered derivative mean and odd velocity**

Replace `src/Phi4Imaginary2.jl` with:

```julia
module Phi4Imaginary2

using QuadGK

export Params, S_bare, S_reg, dSreg_dt, dSreg_dphi, integration_bound
export dSdt_mean, phi_dot_positive, phi_dot, velocity_grid, odd_symmetry_error

Base.@kwdef struct Params
    m2::Float64 = 1.0
    lam::Float64 = 1.0
    eps::Float64 = 1e-3
end

phase(t) = exp(-im * pi * t / 2)

function S_bare(phi, t, p::Params=Params())
    return phase(t) * (p.m2 * phi^2 / 2 + p.lam * phi^4 / 4)
end

function S_reg(phi, t, p::Params=Params())
    return S_bare(phi, t, p) + p.eps * phi^2 / 2
end

function dSreg_dt(phi, t, p::Params=Params())
    return (-im * pi / 2) * S_bare(phi, t, p)
end

function dSreg_dphi(phi, t, p::Params=Params())
    return phase(t) * (p.m2 * phi + p.lam * phi^3) + p.eps * phi
end

function integration_bound(t, p::Params=Params())
    c = max(cos(pi * t / 2), 0.0) + p.eps
    return clamp(sqrt(56.0 / c), 8.0, 80.0)
end

weight(phi, t, p::Params) = exp(-S_reg(phi, t, p))

function dSdt_mean(t, p::Params=Params(); L=integration_bound(t, p), rtol=1e-9)
    num, _ = quadgk(x -> dSreg_dt(x, t, p) * weight(x, t, p), -L, L; rtol=rtol)
    den, _ = quadgk(x -> weight(x, t, p), -L, L; rtol=rtol)
    return num / den
end

function phi_dot_positive(phi_nonnegative::Real, t, p::Params=Params();
                          L=integration_bound(t, p), mean=dSdt_mean(t, p; L=L),
                          rtol=1e-9)
    phi_nonnegative < 0 && error("phi_dot_positive requires phi >= 0")
    phi_nonnegative == 0 && return zero(ComplexF64)
    integrand(x) = (dSreg_dt(x, t, p) - mean) * weight(x, t, p)
    tail, _ = quadgk(integrand, phi_nonnegative, L; rtol=rtol)
    return exp(S_reg(phi_nonnegative, t, p)) * tail
end

function phi_dot(phi::Real, t, p::Params=Params(); kwargs...)
    phi == 0 && return zero(ComplexF64)
    v = phi_dot_positive(abs(phi), t, p; kwargs...)
    return phi > 0 ? v : -v
end

function velocity_grid(phi_grid, t_grid, p::Params=Params())
    out = zeros(ComplexF64, length(t_grid), length(phi_grid))
    for (i, t) in pairs(t_grid)
        L = integration_bound(t, p)
        mean = dSdt_mean(t, p; L=L)
        for (j, phi) in pairs(phi_grid)
            out[i, j] = phi_dot(phi, t, p; L=L, mean=mean)
        end
    end
    return out
end

function odd_symmetry_error(phi_grid, row)
    n = length(phi_grid)
    return maximum(abs(row[j] + row[n + 1 - j]) for j in eachindex(phi_grid))
end

end
```

- [ ] **Step 4: Run tests and verify pass**

Run:

```bash
julia --project=demovault/projects/phi4imaginary2 -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add demovault/projects/phi4imaginary2/src/Phi4Imaginary2.jl demovault/projects/phi4imaginary2/test/runtests.jl
git commit -m "feat: add stable odd velocity"
```

---

### Task 4: Implement Flow and Observable Estimators

**Files:**
- Modify: `demovault/projects/phi4imaginary2/src/Phi4Imaginary2.jl`
- Modify: `demovault/projects/phi4imaginary2/test/runtests.jl`

- [ ] **Step 1: Add failing flow and quadrature tests**

Append to `test/runtests.jl`:

```julia
@testset "flow and observables" begin
    p = Params()
    phi_grid = collect(range(-4.0, 4.0, length=161))
    t_grid = collect(range(0.0, 1.0, length=21))
    vgrid = velocity_grid(phi_grid, t_grid, p)
    traj = integrate_flow([-1.0, 0.0, 1.0], phi_grid, t_grid, vgrid; n_steps=40)
    @test size(traj) == (3, 41)
    @test traj[2, end] == 0
    @test abs(traj[1, end] + traj[3, end]) < 1e-6 * max(1.0, abs(traj[3, end]))
    ref = quadrature_observables(0.25, p)
    @test haskey(ref, :phi2)
    @test haskey(ref, :phi4)
    @test isfinite(real(ref.phi2))
    @test isfinite(imag(ref.phi2))
end
```

- [ ] **Step 2: Run tests and verify failure**

Run:

```bash
julia --project=demovault/projects/phi4imaginary2 -e 'using Pkg; Pkg.test()'
```

Expected: FAIL with `integrate_flow` and `quadrature_observables` undefined.

- [ ] **Step 3: Add flow and quadrature implementation**

Append these exports near the existing export block in `src/Phi4Imaginary2.jl`:

```julia
export interpolate_velocity, rk4_step, integrate_flow, flow_observables, quadrature_observables
```

Append these functions before the final `end`:

```julia
function interpolate_velocity(phi::Complex, t, phi_grid, t_grid, vgrid)
    tr = clamp(real(t), first(t_grid), last(t_grid))
    xr = clamp(real(phi), first(phi_grid), last(phi_grid))
    ti = searchsortedlast(t_grid, tr)
    xi = searchsortedlast(phi_grid, xr)
    ti = clamp(ti, 1, length(t_grid) - 1)
    xi = clamp(xi, 1, length(phi_grid) - 1)
    t0, t1 = t_grid[ti], t_grid[ti + 1]
    x0, x1 = phi_grid[xi], phi_grid[xi + 1]
    at = (tr - t0) / (t1 - t0)
    ax = (xr - x0) / (x1 - x0)
    v00 = vgrid[ti, xi]
    v10 = vgrid[ti + 1, xi]
    v01 = vgrid[ti, xi + 1]
    v11 = vgrid[ti + 1, xi + 1]
    return (1 - at) * ((1 - ax) * v00 + ax * v01) +
           at * ((1 - ax) * v10 + ax * v11)
end

function rk4_step(phi, t, dt, phi_grid, t_grid, vgrid)
    f(z, s) = interpolate_velocity(z, s, phi_grid, t_grid, vgrid)
    k1 = f(phi, t)
    k2 = f(phi + dt * k1 / 2, t + dt / 2)
    k3 = f(phi + dt * k2 / 2, t + dt / 2)
    k4 = f(phi + dt * k3, t + dt)
    return phi + dt * (k1 + 2k2 + 2k3 + k4) / 6
end

function integrate_flow(phi0_values, phi_grid, t_grid, vgrid; n_steps=512)
    t0, t1 = first(t_grid), last(t_grid)
    dt = (t1 - t0) / n_steps
    traj = zeros(ComplexF64, length(phi0_values), n_steps + 1)
    for (i, phi0) in pairs(phi0_values)
        traj[i, 1] = complex(phi0)
        for step in 1:n_steps
            t = t0 + (step - 1) * dt
            traj[i, step + 1] = rk4_step(traj[i, step], t, dt, phi_grid, t_grid, vgrid)
        end
    end
    return traj
end

function flow_observables(final_values)
    phi2 = mean(final_values .^ 2)
    phi4 = mean(final_values .^ 4)
    return (phi2=phi2, phi4=phi4)
end

function quadrature_observables(t, p::Params=Params(); L=integration_bound(t, p), rtol=1e-9)
    w(x) = weight(x, t, p)
    den, _ = quadgk(w, -L, L; rtol=rtol)
    n2, _ = quadgk(x -> x^2 * w(x), -L, L; rtol=rtol)
    n4, _ = quadgk(x -> x^4 * w(x), -L, L; rtol=rtol)
    qden, _ = quadgk(x -> exp(-real(S_reg(x, t, p))), -L, L; rtol=rtol)
    return (phi2=n2/den, phi4=n4/den, avg_phase=den/qden)
end
```

- [ ] **Step 4: Add `using Statistics`**

Add this near the top of `src/Phi4Imaginary2.jl`:

```julia
using Statistics
```

- [ ] **Step 5: Run tests and verify pass**

Run:

```bash
julia --project=demovault/projects/phi4imaginary2 -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add demovault/projects/phi4imaginary2/src/Phi4Imaginary2.jl demovault/projects/phi4imaginary2/test/runtests.jl
git commit -m "feat: add flow and observables"
```

---

### Task 5: Add Numbered Pipeline Scripts

**Files:**
- Create: `demovault/projects/phi4imaginary2/scripts/01_velocity.jl`
- Create: `demovault/projects/phi4imaginary2/scripts/02_flow.jl`
- Create: `demovault/projects/phi4imaginary2/scripts/03_observables.jl`
- Create: `demovault/projects/phi4imaginary2/scripts/04_summary.jl`

- [ ] **Step 1: Write `scripts/01_velocity.jl`**

```julia
using Phi4Imaginary2, JLD2, Plots

ROOT = normpath(joinpath(@__DIR__, ".."))
DATA = joinpath(ROOT, "data")
mkpath(DATA)

p = Params()
phi_grid = collect(range(-5.0, 5.0, length=161))
t_grid = collect(range(0.0, 1.0, length=65))

vgrid = velocity_grid(phi_grid, t_grid, p)
odd_errors = [odd_symmetry_error(phi_grid, vgrid[i, :]) for i in eachindex(t_grid)]

jldsave(joinpath(DATA, "velocity.jld2");
    phi_grid, t_grid, phi_dot=vgrid, odd_errors,
    m2=p.m2, lam=p.lam, eps=p.eps, script="01_velocity.jl")

p1 = plot(title="Re velocity, Tikhonov real-axis", xlabel="phi", ylabel="Re phi_dot")
p2 = plot(title="Im velocity, Tikhonov real-axis", xlabel="phi", ylabel="Im phi_dot")
for tv in (0.0, 0.25, 0.5, 0.75, 1.0)
    i = argmin(abs.(t_grid .- tv))
    plot!(p1, phi_grid, real.(vgrid[i, :]), label="t=$(round(t_grid[i], digits=2))")
    plot!(p2, phi_grid, imag.(vgrid[i, :]), label="t=$(round(t_grid[i], digits=2))")
end
savefig(plot(p1, p2, layout=(1, 2), size=(1000, 420)), joinpath(DATA, "01_velocity.png"))

println("Saved velocity.jld2 and 01_velocity.png")
println("max odd symmetry error = ", maximum(odd_errors))
```

- [ ] **Step 2: Write `scripts/02_flow.jl`**

```julia
using Phi4Imaginary2, JLD2, Plots

ROOT = normpath(joinpath(@__DIR__, ".."))
DATA = joinpath(ROOT, "data")

d = load(joinpath(DATA, "velocity.jld2"))
phi_grid = d["phi_grid"]
t_grid = d["t_grid"]
vgrid = d["phi_dot"]

phi0_values = collect(range(-3.0, 3.0, length=121))
traj = integrate_flow(phi0_values, phi_grid, t_grid, vgrid; n_steps=512)
t_steps = collect(range(first(t_grid), last(t_grid), length=size(traj, 2)))

jldsave(joinpath(DATA, "flow.jld2");
    phi0_values, trajectories=traj, t_steps,
    m2=d["m2"], lam=d["lam"], eps=d["eps"], script="02_flow.jl")

p = plot(title="Transported real-axis samples", xlabel="Re phi", ylabel="Im phi",
         aspect_ratio=:equal, legend=false)
for i in axes(traj, 1)
    plot!(p, real.(traj[i, :]), imag.(traj[i, :]), alpha=0.45, linewidth=1)
end
scatter!(p, real.(traj[:, 1]), imag.(traj[:, 1]), color=:blue, markersize=2)
scatter!(p, real.(traj[:, end]), imag.(traj[:, end]), color=:red, markersize=2)
savefig(p, joinpath(DATA, "02_flow.png"))

println("Saved flow.jld2 and 02_flow.png")
```

- [ ] **Step 3: Write `scripts/03_observables.jl`**

```julia
using Phi4Imaginary2, JLD2, Plots

ROOT = normpath(joinpath(@__DIR__, ".."))
DATA = joinpath(ROOT, "data")

p = Params()
flow = load(joinpath(DATA, "flow.jld2"))
traj = flow["trajectories"]
final_values = traj[:, end]

flow_obs = flow_observables(final_values)
t_ref = collect(range(0.0, 1.0, length=81))
quad = [quadrature_observables(t, p) for t in t_ref]
quad_phi2 = [q.phi2 for q in quad]
quad_phi4 = [q.phi4 for q in quad]
quad_phase = [q.avg_phase for q in quad]

jldsave(joinpath(DATA, "observables.jld2");
    flow_phi2=flow_obs.phi2,
    flow_phi4=flow_obs.phi4,
    t_ref, quad_phi2, quad_phi4, quad_phase,
    m2=p.m2, lam=p.lam, eps=p.eps, script="03_observables.jl")

p1 = plot(t_ref, real.(quad_phi2), label="Re quadrature phi2", xlabel="t", ylabel="phi2",
          title="Quadrature reference vs transported endpoint")
plot!(p1, t_ref, imag.(quad_phi2), label="Im quadrature phi2", linestyle=:dash)
hline!(p1, [real(flow_obs.phi2)], label="Re flow endpoint", color=:red)
hline!(p1, [imag(flow_obs.phi2)], label="Im flow endpoint", color=:orange, linestyle=:dash)

p2 = plot(t_ref, abs.(quad_phase), label="quadrature avg phase", xlabel="t", ylabel="abs phase",
          title="Tikhonov sign severity")
savefig(plot(p1, p2, layout=(1, 2), size=(1100, 420)), joinpath(DATA, "03_observables.png"))

println("Saved observables.jld2 and 03_observables.png")
println("flow endpoint phi2 = ", flow_obs.phi2)
println("quadrature t=1 phi2 = ", quad_phi2[end])
println("absolute discrepancy = ", abs(flow_obs.phi2 - quad_phi2[end]))
```

- [ ] **Step 4: Write `scripts/04_summary.jl`**

```julia
using JLD2, Plots, Printf

ROOT = normpath(joinpath(@__DIR__, ".."))
DATA = joinpath(ROOT, "data")

vel = load(joinpath(DATA, "velocity.jld2"))
flow = load(joinpath(DATA, "flow.jld2"))
obs = load(joinpath(DATA, "observables.jld2"))

println("=== phi4imaginary2 summary ===")
@printf("regulator: Tikhonov eps = %.3e\n", vel["eps"])
@printf("max velocity odd-symmetry error: %.6e\n", maximum(vel["odd_errors"]))
@printf("flow endpoint phi2: %+.8f %+.8fim\n", real(obs["flow_phi2"]), imag(obs["flow_phi2"]))
@printf("quadrature t=1 phi2: %+.8f %+.8fim\n", real(obs["quad_phi2"][end]), imag(obs["quad_phi2"][end]))
@printf("endpoint discrepancy: %.6e\n", abs(obs["flow_phi2"] - obs["quad_phi2"][end]))

phi_grid = vel["phi_grid"]
t_grid = vel["t_grid"]
vgrid = vel["phi_dot"]
flow_traj = flow["trajectories"]
t_ref = obs["t_ref"]
quad_phi2 = obs["quad_phi2"]

p1 = plot(title="Velocity Im phi_dot", xlabel="phi", ylabel="Im phi_dot")
for tv in (0.0, 0.5, 1.0)
    i = argmin(abs.(t_grid .- tv))
    plot!(p1, phi_grid, imag.(vgrid[i, :]), label="t=$(round(t_grid[i], digits=2))")
end

p2 = plot(title="Transported samples", xlabel="Re phi", ylabel="Im phi", aspect_ratio=:equal, legend=false)
for i in axes(flow_traj, 1)
    i % 6 == 1 && plot!(p2, real.(flow_traj[i, :]), imag.(flow_traj[i, :]), alpha=0.45, linewidth=1)
end

p3 = plot(t_ref, real.(quad_phi2), label="Re quadrature phi2", xlabel="t", ylabel="phi2",
          title="Reference and flow endpoint")
plot!(p3, t_ref, imag.(quad_phi2), label="Im quadrature phi2", linestyle=:dash)
hline!(p3, [real(obs["flow_phi2"])], label="Re flow endpoint", color=:red)
hline!(p3, [imag(obs["flow_phi2"])], label="Im flow endpoint", color=:orange, linestyle=:dash)

savefig(plot(p1, p2, p3, layout=(1, 3), size=(1500, 420)), joinpath(DATA, "04_summary.png"))
println("Saved 04_summary.png")
```

- [ ] **Step 5: Run scripts in order**

Run:

```bash
julia --project=demovault/projects/phi4imaginary2 demovault/projects/phi4imaginary2/scripts/01_velocity.jl
julia --project=demovault/projects/phi4imaginary2 demovault/projects/phi4imaginary2/scripts/02_flow.jl
julia --project=demovault/projects/phi4imaginary2 demovault/projects/phi4imaginary2/scripts/03_observables.jl
julia --project=demovault/projects/phi4imaginary2 demovault/projects/phi4imaginary2/scripts/04_summary.jl
```

Expected:

- `01_velocity.jl` prints a small max odd symmetry error.
- `02_flow.jl` writes finite trajectories.
- `03_observables.jl` prints flow endpoint, quadrature reference, and discrepancy.
- `04_summary.jl` prints metadata and discrepancy.

- [ ] **Step 6: Commit**

```bash
git add demovault/projects/phi4imaginary2/scripts demovault/projects/phi4imaginary2/data
git commit -m "feat: add phi4imaginary2 pipeline"
```

---

### Task 6: Add README and Final Validation

**Files:**
- Create: `demovault/projects/phi4imaginary2/README.md`
- Modify: `demovault/projects/phi4imaginary2/test/runtests.jl`

- [ ] **Step 1: Write README**

Create `demovault/projects/phi4imaginary2/README.md`:

```markdown
# phi4imaginary2

Clean Tikhonov-regularised parallel implementation of the zero-dimensional imaginary phi4 benchmark.

This project intentionally does not modify or depend on `../phi4imaginary`.

## Scope

- Tikhonov-only action:
  `S_reg(phi,t) = exp(-i*pi*t/2) * (m2*phi^2/2 + lam*phi^4/4) + eps*phi^2/2`
- Real-axis Wegner velocity for real initial coordinates.
- Velocity is enforced and tested as odd: `phi_dot(-phi,t) = -phi_dot(phi,t)`.
- Transported-flow estimator is the primary output.
- Direct Tikhonov real-axis quadrature is a reference benchmark.

## Limitations

- This is not a full complex-plane velocity field.
- This does not claim to compute a PIKfold or Lefschetz thimble.
- Non-holomorphic observables are excluded from the primary result.

## Run

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=. scripts/01_velocity.jl
julia --project=. scripts/02_flow.jl
julia --project=. scripts/03_observables.jl
julia --project=. scripts/04_summary.jl
```

Run from `demovault/projects/phi4imaginary2` or use absolute script paths.
```

- [ ] **Step 2: Add checkpoint metadata test**

Append to `test/runtests.jl`:

```julia
@testset "integration bound" begin
    p = Params()
    @test integration_bound(0.0, p) >= 8.0
    @test integration_bound(1.0, p) <= 80.0
    @test integration_bound(1.0, p) > integration_bound(0.0, p)
end
```

- [ ] **Step 3: Run full tests**

Run:

```bash
julia --project=demovault/projects/phi4imaginary2 -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 4: Verify existing project untouched**

Run:

```bash
git status --short demovault/projects/phi4imaginary demovault/projects/phi4imaginary2
```

Expected: only `phi4imaginary2` files appear as new/modified from this implementation session; existing `phi4imaginary` changes must not be introduced by this work.

- [ ] **Step 5: Commit**

```bash
git add demovault/projects/phi4imaginary2/README.md demovault/projects/phi4imaginary2/test/runtests.jl
git commit -m "docs: document phi4imaginary2 scope"
```

---

## Final Verification Checklist

- [ ] `julia --project=demovault/projects/phi4imaginary2 -e 'using Pkg; Pkg.test()'` passes.
- [ ] `scripts/01_velocity.jl` reports small odd-symmetry error and finite velocity values.
- [ ] `scripts/02_flow.jl` reports finite trajectories.
- [ ] `scripts/03_observables.jl` writes both flow and quadrature values.
- [ ] `scripts/04_summary.jl` prints the flow-vs-reference discrepancy.
- [ ] No file under `demovault/projects/phi4imaginary` is edited by this implementation.
- [ ] README clearly states Tikhonov-only and real-axis-only limitations.

## Known Design Tradeoff

The real-axis velocity is used to transport real initial coordinates and is interpolated by `real(phi)` during RK4. This is acceptable for `phi4imaginary2` only because the scope explicitly excludes a full complex-plane kernel. The README and plot labels must avoid claiming "PIKfold", "thimble", or "full complex flow".
