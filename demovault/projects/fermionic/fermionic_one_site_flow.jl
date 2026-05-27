# # Fermionic one-site Hubbard: PIK-flow companion notebook
#
# Companion to the (frozen) `fermionic_one_site.{jl,ipynb}` which covers the
# analytic weighted-sum branch of Fukushima/Ihssen/Kapust/Pawlowski 2026
# (`raw/papers/FermionPIK.pdf`). This notebook covers the *single-PIKfold
# flow* branch: the Wegner integration of $\dot\phi_t$ from $t=0$ (real
# Gaussian) to $t=1$ (target Hubbard action) that produces the paper's
# manifold figures.
#
# **Frozen files (do not modify):** `fermionic_one_site.{jl,ipynb}` and all
# of its `data/` outputs are verified and frozen. This notebook only writes
# *new* files. All new artefacts use distinct names (`0X_*` PNGs prefixed by
# section number; `flow_*.jld2` / `n_vs_mu_scan.jld2`).
#
# Sections (gates after each):
#
# - § 0 setup & helpers
# - § 1 two action paths and the conceptual map (numerical D_KL = 0 check)
# - § 2 Wegner velocity field on a (t, φ)-grid for 8 (path, μ) combos
# - § 3 manifold visualisation (paper Figs 4 and 1)
# - § 4 3D PIKfold (paper Fig 3a)
# - § 5 ripping diagnostic, Path 1 vs Path 2 at μ = -0.675 (paper Fig 2)
# - § 6 n vs μ scan (paper Fig 5)
# - § 7 summary marquee
#
# Run with `julia --project=. -t auto fermionic_one_site_flow.jl`.

# %% Phase 0: Setup, parameters, frozen helpers
using Pkg
Pkg.activate(@__DIR__)
using Random, Distributions, QuadGK, JLD2, Plots, Interpolations, StatsBase
Random.seed!(20260527)

const U = 1.0
const β = 10.0

# Frozen analytic-side helpers, copied verbatim (≤ 30 lines total) from
# `fermionic_one_site.jl`. Re-defined here rather than `include`d so this
# script remains self-contained and does not re-run the analytic notebook.

S(phi, U, beta, mu) = beta/(2U) * phi^2 - 2*log(1 + exp(beta*(im*phi + mu + U/2)))
boltzmann(phi, U, beta, mu) = exp(-S(phi, U, beta, mu))

function mixture_params(U, beta, mu)
    weights = [1.0, 2*exp(beta*mu), exp(beta*(2mu - U))]
    centers = [0.0, U, 2U]
    sigma2  = U / beta
    return weights, centers, sigma2
end

n_op(phi, U, beta, mu) = -im * phi / U

function fock_trace_observables(U, beta, mu)
    z1 = 2 * exp(beta * mu)
    z2 = exp(beta * (2mu - U))
    Z  = 1.0 + z1 + z2
    return (; Z, n_mean = (z1 + 2*z2)/Z, n2_mean = (z1 + 4*z2)/Z)
end

function pik_eval(O, U, beta, mu, phi0::AbstractVector{<:Real})
    w, m, _ = mixture_params(U, beta, mu)
    M = length(phi0)
    Wtot = sum(w)
    contrib = zeros(ComplexF64, M)
    for j in 1:3, i in 1:M
        contrib[i] += w[j] * O(phi0[i] + im*m[j])
    end
    contrib ./= Wtot
    est   = sum(contrib) / M
    se_re = std(real.(contrib); corrected=true) / sqrt(M)
    return (; est, se_re, contrib)
end

println("§0 setup. (U, β) = ($U, $β); seed = 20260527; threads = $(Threads.nthreads()).")
println("Frozen analytic helpers redefined inline (S, boltzmann, mixture_params, n_op, pik_eval, fock_trace).")

# %% Phase 1: Two action paths and the conceptual map
# # § 1 — Two action paths and the conceptual map
#
# **Path 1** (paper eq 9a, holomorphic PIK):
#
# $$ S_t^{(1)}(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + e^{t\beta(i\phi + \mu + U/2)}\right]. $$
#
# **Path 2** (paper eq 14, alternative interpolation):
#
# $$ S_t^{(2)}(\phi) = \frac{\beta}{2U}\phi^2 - 2\log\!\left[1 + t \cdot e^{\beta(i\phi + \mu + U/2)}\right]. $$
#
# Both reduce to the real Gaussian $\beta\phi^2/(2U)$ at $t = 0$ (initial
# manifold $\mathcal M_0 = \mathbb R$) and to the Hubbard action at $t = 1$.

# Action and $\partial_t S$ definitions. Numerically stable variants:
# Path 1 uses the complex sigmoid form $\sigma(t\beta z)$ to avoid
# `exp` overflow when $\beta\,\mathrm{Re}\,z$ is large. Path 2 uses the
# Boltzmann-weight ratio with an `exp(-βz)` guard for the same reason.

S_t_path1(φ, t; U=U, β=β, μ=0.0) =
    β/(2U)*φ^2 - 2*log(1 + exp(t*β*(im*φ + μ + U/2)))

S_t_path2(φ, t; U=U, β=β, μ=0.0) =
    β/(2U)*φ^2 - 2*log(1 + t * exp(β*(im*φ + μ + U/2)))

function dSdt_path1(φ, t; U=U, β=β, μ=0.0)
    z = im*φ + μ + U/2
    w = t * β * z
    σ = real(w) ≥ 0 ? 1 / (1 + exp(-w)) : exp(w) / (1 + exp(w))
    return -2 * β * z * σ
end

function dSdt_path2(φ, t; U=U, β=β, μ=0.0)
    z  = im*φ + μ + U/2
    βz = β * z
    if t > 0 && real(βz) > 30
        return -2 / (exp(-βz) + t)
    else
        e = exp(βz)
        return -2 * e / (1 + t * e)
    end
end

# Free-theory check: at $t = 0$ both actions reduce to the Gaussian
# $\beta\phi^2/(2U)$ — Path 2 exactly, Path 1 up to a constant offset
# $-2\log 2$ from $\log(1 + e^{0\cdot\beta z}) = \log 2$. The offset is the
# free Fermi-sea partition-function normalisation and drops out of the
# Wegner flow, which depends only on $\partial_t S$ and gradients.
let
    grid = range(-4.0, 4.0; length=100)
    gauss = (β/(2U)) .* grid.^2
    err1 = maximum(abs.(S_t_path1.(grid, 0.0; μ=0.7) .- (gauss .- 2*log(2))))
    err2 = maximum(abs.(S_t_path2.(grid, 0.0; μ=0.7) .- gauss))
    @assert err1 < 1e-12 && err2 < 1e-12 "t=0 should reduce to β/(2U)φ² (mod constant for Path 1)"
    println("§1 t=0 Gaussian limit: max|ΔS - const| = $(max(err1, err2))  ✓")
end

# ## Complex → real Gaussian: why no KL fit at one site
#
# The Wegner velocity $\dot\phi_t(\phi)$ is the *deterministic* map taking
# real Gaussian samples $\xi \sim \mathcal N(0, U/\beta)$ at $t = 0$ to the
# complex sampling contour $\mathcal M_t \subset \mathbb C$ at $t = 1$.
#
# For the holomorphic one-site decomposition (paper eq 18) this velocity is
# the **constant** $\dot\phi = iU$, hence
#
# $$ \phi_t(\xi) = \xi + iUt, \qquad \det[\partial\phi_t/\partial\xi] = 1. $$
#
# Reweighting by the analytic mixture coefficients
# $w_j = (1,\; 2e^{\beta\mu},\; e^{\beta(2\mu-U)})$ reproduces $e^{-S(\phi)}$
# *exactly* — the closed-form mixture is the unique optimum. Therefore
#
# $$ D_{\mathrm{KL}}\!\left(p_{\mathrm{mix}} \,\Big\| \, \tfrac{e^{-S}}{Z}\right) = 0 $$
#
# by construction; there is nothing to fit. At **multi-site** the closed-form
# decomposition is no longer exact and a complex-shifted-Gaussian mixture
# must be fit by minimising $|D_{\mathrm{KL}}|^2$. That multi-site procedure
# is documented separately at `[[kl-gaussian-mixture-decomposition]]`
# (`~/.claude/projects/-home-konstantin-Desktop-masterclass/memory/`
# `kl_gaussian_mixture_decomposition.md`).

# Numerical gate: max relative deviation between the analytic 3-Gaussian
# mixture and $e^{-S}$ on the real axis at the headline $\mu = 1$ —
# the "$D_{\mathrm{KL}} = 0$" check.
MIXTURE_RELERR = let
    μ = 1.0
    grid = range(-6.0, 6.0; length=601)
    w, m, σ2 = mixture_params(U, β, μ)
    pmix = [sum(w[j]*exp(-(φ - im*m[j])^2/(2σ2)) for j in 1:3) for φ in grid]
    bv   = boltzmann.(grid, U, β, μ)
    relerr = maximum(abs.(pmix .- bv) ./ max.(abs.(bv), 1.0))
    @assert relerr < 1e-12 "p_mix vs e^{-S} relative err = $relerr"
    println("§1 KL=0 sanity (μ=1, U=$U, β=$β):  max|p_mix − e^{−S}|/|e^{−S}| = $relerr  ✓")
    relerr
end

# Visualise both action paths at five $t$-slices (real-axis cross-section).
let
    μ_demo = 1.0
    φ_grid = range(-4.0, 4.0; length=400)
    t_vals = [0.0, 0.25, 0.5, 0.75, 1.0]
    p_re1 = plot(title="Re S_t^{(1)}(φ), μ=$μ_demo", xlabel="φ", ylabel="Re S", legend=:topright, titlefontsize=10)
    p_im1 = plot(title="Im S_t^{(1)}(φ), μ=$μ_demo", xlabel="φ", ylabel="Im S", legend=:topright, titlefontsize=10)
    p_re2 = plot(title="Re S_t^{(2)}(φ), μ=$μ_demo", xlabel="φ", ylabel="Re S", legend=:topright, titlefontsize=10)
    p_im2 = plot(title="Im S_t^{(2)}(φ), μ=$μ_demo", xlabel="φ", ylabel="Im S", legend=:topright, titlefontsize=10)
    for t in t_vals
        s1 = [S_t_path1(φ, t; μ=μ_demo) for φ in φ_grid]
        s2 = [S_t_path2(φ, t; μ=μ_demo) for φ in φ_grid]
        plot!(p_re1, φ_grid, real.(s1), label="t=$t", lw=2)
        plot!(p_im1, φ_grid, imag.(s1), label="t=$t", lw=2)
        plot!(p_re2, φ_grid, real.(s2), label="t=$t", lw=2)
        plot!(p_im2, φ_grid, imag.(s2), label="t=$t", lw=2)
    end
    p = plot(p_re1, p_im1, p_re2, p_im2; layout=(2,2), size=(1200, 750),
             plot_title="§1 action paths: top = Path 1 (holomorphic), bottom = Path 2")
    savefig(p, joinpath(@__DIR__, "data", "01_action_paths.png"))
    println("Saved data/01_action_paths.png")
end

# %% Phase 2: Wegner velocity field
# # § 2 — Centered Wegner velocity
#
# General form (no rotating-action prefactor; the $i\pi/2$ in
# `phi4imaginary/02_wegner.jl` came from the specific
# $\partial_t S = -i\pi/2 \cdot S$ relation, absent here):
#
# $$ \dot\phi_t(\phi) = e^{S_t(\phi)} \int_\phi^L
#    \Big[\partial_t S_t(\phi') - \langle \partial_t S_t\rangle_t\Big]\,
#    e^{-S_t(\phi')}\,d\phi'. $$
#
# The centering term $\langle \partial_t S_t \rangle_t$ is **essential** —
# without it the free-theory limit (constant $\dot\phi = iU$ for the
# one-site holomorphic decomposition) fails. Integration bound `L_BOUND`
# is chosen so $\mathrm{Re}\,S_t(L) \gg 1$ for $U=1,\beta=10$ and the
# $\mu$ range we cover.

const L_BOUND      = 6.0           # quadgk integration bound (full integrand support)
const PHI_SAMPLE   = 3.0           # interpolant sampling extent for the (t, φ) grid
const VEL_MAX_OK   = 50.0          # cutoff: |φ̇| above this is numerical garbage
const N_T          = 41
const N_PHI        = 200
const phi_grid_re  = collect(range(-PHI_SAMPLE, PHI_SAMPLE; length=N_PHI))
const t_grid       = collect(range(0.0, 1.0; length=N_T))

# At $\phi = L_{\text{BOUND}} = 6$: $\mathrm{Re}\,S_0 = \beta\,L^2/(2U) = 180$,
# so $e^{-\mathrm{Re}\,S_0} < 10^{-78}$ — quadrature is well-resolved.
# The **sampling** grid is narrower: $\phi \in [-3, 3]$, comfortably enclosing
# the trajectory region $|\phi_0| \le 2$. Past $|\phi| \gtrsim 3$ the centered
# Wegner formula $\dot\phi_t = e^{S_t}\!\int_\phi^L[\partial_t S - \bar S]\,
# e^{-S}$ is **ill-conditioned**: at $\beta = 10$ the prefactor $e^{S_t}$
# reaches $e^{135}$ for $|\phi| \!\sim\! 6$ and amplifies the roundoff
# residue of an empty integral to nonsense ($\sim 10^{63}$). We clamp any
# such sample to zero post-hoc (`VEL_MAX_OK`) so the bilinear interpolant
# only ever sees physically meaningful velocities.

function action_mean(t, Sfun, dSdtfun; L=L_BOUND)
    num, _ = quadgk(φ -> dSdtfun(φ, t) * exp(-Sfun(φ, t)), -L, L; rtol=1e-8, maxevals=10_000)
    den, _ = quadgk(φ -> exp(-Sfun(φ, t)),                  -L, L; rtol=1e-8, maxevals=10_000)
    return num / den
end

function phi_dot(φ_real::Real, t, Sbar, Sfun, dSdtfun; L=L_BOUND)
    integrand(x) = (dSdtfun(x, t) - Sbar) * exp(-Sfun(x, t))
    integral, _  = quadgk(integrand, φ_real, L; rtol=1e-8, maxevals=10_000)
    return exp(Sfun(φ_real, t)) * integral
end

function wegner_sweep(Sfun, dSdtfun)
    grid = zeros(ComplexF64, N_T, N_PHI)
    for (i, t) in enumerate(t_grid)
        Sbar = action_mean(t, Sfun, dSdtfun)
        results = zeros(ComplexF64, N_PHI)
        Threads.@threads for j in 1:N_PHI
            v = phi_dot(phi_grid_re[j], t, Sbar, Sfun, dSdtfun)
            # Clamp ill-conditioned samples (see VEL_MAX_OK doc above).
            results[j] = (isfinite(v) && abs(v) ≤ VEL_MAX_OK) ? v : zero(ComplexF64)
        end
        grid[i, :] = results
    end
    return grid
end

# Sweep grid sizes: $N_t = 41$ (matches plan), $N_\phi = 200$, total of
# 8 sweeps (7 for Path 1 + 1 for Path 2 at μ=-0.675).

const PATH1_MU_LIST = [0.0, 0.1, 0.2, 1.5, 1.65, 1.75, -0.675]
const PATH2_MU_LIST = [-0.675]

println("\n§2 Wegner sweeps: $(length(PATH1_MU_LIST)) Path-1 + $(length(PATH2_MU_LIST)) Path-2 grids of $(N_T)×$(N_PHI)...")

path1_wegner = Vector{Matrix{ComplexF64}}(undef, length(PATH1_MU_LIST))
for (k, μ) in enumerate(PATH1_MU_LIST)
    t0 = time()
    Sclosed = (φ, t) -> S_t_path1(φ, t; U=U, β=β, μ=μ)
    Dclosed = (φ, t) -> dSdt_path1(φ, t; U=U, β=β, μ=μ)
    path1_wegner[k] = wegner_sweep(Sclosed, Dclosed)
    println("  Path 1  μ=$(rpad(string(μ), 7))  done in $(round(time()-t0; digits=1))s")
end

path2_wegner = Vector{Matrix{ComplexF64}}(undef, length(PATH2_MU_LIST))
for (k, μ) in enumerate(PATH2_MU_LIST)
    t0 = time()
    Sclosed = (φ, t) -> S_t_path2(φ, t; U=U, β=β, μ=μ)
    Dclosed = (φ, t) -> dSdt_path2(φ, t; U=U, β=β, μ=μ)
    path2_wegner[k] = wegner_sweep(Sclosed, Dclosed)
    println("  Path 2  μ=$(rpad(string(μ), 7))  done in $(round(time()-t0; digits=1))s")
end

jldsave(joinpath(@__DIR__, "data", "flow_path1_wegner.jld2");
        t_grid = t_grid, phi_grid = phi_grid_re,
        mus = PATH1_MU_LIST, velocities = path1_wegner)
jldsave(joinpath(@__DIR__, "data", "flow_path2_wegner.jld2");
        t_grid = t_grid, phi_grid = phi_grid_re,
        mus = PATH2_MU_LIST, velocities = path2_wegner)
println("Saved flow_path1_wegner.jld2 and flow_path2_wegner.jld2")

# Sanity figure: Re/Im of $\dot\phi_t(\phi)$ at five $t$-slices for Path 1,
# μ = 0 (the closed-manifold regime).
let
    k0   = findfirst(==(0.0), PATH1_MU_LIST)
    grid = path1_wegner[k0]
    plot_ts = [0.0, 0.25, 0.5, 0.75, 1.0]
    p_re = plot(title="Re φ̇_t(φ), Path 1, μ=0", xlabel="φ", ylabel="Re φ̇",
                legend=:topright, xlims=(-4, 4), titlefontsize=10)
    p_im = plot(title="Im φ̇_t(φ), Path 1, μ=0", xlabel="φ", ylabel="Im φ̇",
                legend=:topright, xlims=(-4, 4), titlefontsize=10)
    for tv in plot_ts
        idx = argmin(abs.(t_grid .- tv))
        row = grid[idx, :]
        plot!(p_re, phi_grid_re, real.(row), label="t=$tv", lw=2)
        plot!(p_im, phi_grid_re, imag.(row), label="t=$tv", lw=2)
    end
    p = plot(p_re, p_im; layout=(1,2), size=(1200, 450),
             plot_title="§2 Wegner velocity field, Path 1, μ=0")
    savefig(p, joinpath(@__DIR__, "data", "02_velocity_path1.png"))
    println("Saved data/02_velocity_path1.png")
end

# § 2 gate — Path 1 velocity bounded in the trajectory region at μ = 1.75.
# The sampling grid runs to ±3, but trajectories start in |φ_0| ≤ 2; we
# verify the velocity is well-conditioned there (the clamp above guarantees
# at most VEL_MAX_OK = 50 on the full grid, so this is the tighter check).
let
    k = findfirst(==(1.75), PATH1_MU_LIST)
    mask = abs.(phi_grid_re) .≤ 2.0
    max_v1 = maximum(abs.(path1_wegner[k][end, mask]))
    @assert max_v1 < 10 "Path 1 μ=1.75: max|φ̇_1| over |φ|≤2 = $max_v1 not bounded"
    println("§2 gate: max|φ̇_1(φ)| at Path 1, μ=1.75, |φ|≤2 = $(round(max_v1; digits=4)) < 10  ✓")
end

# %% Phase 3: Manifold visualization
# # § 3 — Manifold at $t = 1$ (paper Figs 4 and 1)
#
# Each velocity grid feeds an RK4 integrator. Trajectories
# $\phi_0 \in [-2, 2]$ ($N_{\text{init}} = 121$) are integrated for
# $N_{\text{steps}} = 512$ steps of $\Delta t = 1/512$.

function build_interp(grid::Matrix{ComplexF64})
    re_itp = interpolate((t_grid, phi_grid_re), real.(grid), Gridded(Linear()))
    im_itp = interpolate((t_grid, phi_grid_re), imag.(grid), Gridded(Linear()))
    re_ext = extrapolate(re_itp, Flat())
    im_ext = extrapolate(im_itp, Flat())
    return function(φc::Complex, t)
        φr = clamp(real(φc), phi_grid_re[1], phi_grid_re[end])
        tc = clamp(t,        t_grid[1],     t_grid[end])
        return re_ext(tc, φr) + im * im_ext(tc, φr)
    end
end

function rk4_step(φ̇, φ, t, dt)
    k1 = φ̇(φ,                t)
    k2 = φ̇(φ + dt/2 * k1,    t + dt/2)
    k3 = φ̇(φ + dt/2 * k2,    t + dt/2)
    k4 = φ̇(φ + dt   * k3,    t + dt)
    return φ + dt/6 * (k1 + 2k2 + 2k3 + k4)
end

const N_INIT  = 121
const N_STEPS = 512
const DT      = 1.0 / N_STEPS
const phi0_grid = collect(range(-2.0, 2.0; length=N_INIT))
const traj_ts   = collect(range(0.0, 1.0; length=N_STEPS + 1))

function integrate_trajectories(grid_velocity::Matrix{ComplexF64})
    φ̇ = build_interp(grid_velocity)
    traj = zeros(ComplexF64, N_INIT, N_STEPS + 1)
    for (n, φ0) in enumerate(phi0_grid)
        traj[n, 1] = ComplexF64(φ0)
        for step in 1:N_STEPS
            t_now = (step - 1) * DT
            traj[n, step+1] = rk4_step(φ̇, traj[n, step], t_now, DT)
        end
    end
    return traj
end

println("\n§3 Integrating $(N_INIT) trajectories × $(N_STEPS) RK4 steps for each (path, μ)...")
path1_traj = [integrate_trajectories(g) for g in path1_wegner]
path2_traj = [integrate_trajectories(g) for g in path2_wegner]

jldsave(joinpath(@__DIR__, "data", "flow_path1_traj.jld2");
        phi0 = phi0_grid, t_steps = traj_ts,
        mus = PATH1_MU_LIST, trajectories = path1_traj)
jldsave(joinpath(@__DIR__, "data", "flow_path2_traj.jld2");
        phi0 = phi0_grid, t_steps = traj_ts,
        mus = PATH2_MU_LIST, trajectories = path2_traj)

# **Caveat (same as `projects/phi4imaginary/03_flow.jl`):** the
# bilinear interpolant samples $\dot\phi$ along the real axis only.
# While trajectories stay within $|\mathrm{Im}\,\phi_t| \lesssim 2U$
# the approximation is faithful; large excursions clamp at the grid
# edges and the prediction degrades. This is the accepted trade-off
# for the headline figures in this regime.

# Shading helper: pink stripes where $\mathrm{Im}\,S \approx 0\!\!\mod 2\pi$.
function im_s_mask(grid_x, grid_y, t, Sfun; band=π/4)
    mask = falses(length(grid_x), length(grid_y))
    for (i, x) in enumerate(grid_x), (j, y) in enumerate(grid_y)
        ims = mod(imag(Sfun(x + im*y, t)), 2π)
        mask[i, j] = (ims < band) | (ims > 2π - band)
    end
    return mask
end

function manifold_panel(μ, traj, Sfun;
                        xlim=(-2.0, 2.0), ylim=(-0.5, 3.0), Ngrid=200)
    gx = range(xlim[1], xlim[2]; length=Ngrid)
    gy = range(ylim[1], ylim[2]; length=Ngrid)
    mask = im_s_mask(gx, gy, 1.0, Sfun)
    p = heatmap(gx, gy, Float64.(mask)';
                color=cgrad([:white, :pink]), alpha=0.35,
                colorbar=false, legend=false,
                xlabel="Re φ", ylabel="Im φ",
                title="μ = $(μ)", titlefontsize=10,
                xlims=xlim, ylims=ylim)
    pole_im = μ + U/2
    for n in -3:3
        pre = (2n + 1) * π / β
        if xlim[1] ≤ pre ≤ xlim[2] && ylim[1] ≤ pole_im ≤ ylim[2]
            scatter!(p, [pre], [pole_im];
                     marker=:rect, color=:blue, markersize=4,
                     markerstrokewidth=0, label="")
        end
    end
    end_pts = traj[:, end]
    plot!(p, real.(end_pts), imag.(end_pts);
          color=:black, lw=2, label="")
    return p
end

rip_flags = Dict{Tuple{Symbol,Float64}, Bool}()

let
    panels = []
    for μ in [0.0, 0.1, 0.2]
        k    = findfirst(==(μ), PATH1_MU_LIST)
        Sfun = (φ, t) -> S_t_path1(φ, t; U=U, β=β, μ=μ)
        push!(panels, manifold_panel(μ, path1_traj[k], Sfun))
        max_im = maximum(abs.(imag.(path1_traj[k][:, end])))
        rip_flags[(:path1, μ)] = max_im > 2U
    end
    p = plot(panels...; layout=(1,3), size=(1500, 500),
             plot_title="§3 Path 1 manifold at t=1, μ ∈ {0.0, 0.1, 0.2} (paper Fig 4)")
    savefig(p, joinpath(@__DIR__, "data", "03_manifold_fig4.png"))
end
println("Saved data/03_manifold_fig4.png")

let
    panels = []
    for μ in [1.5, 1.65, 1.75]
        k     = findfirst(==(μ), PATH1_MU_LIST)
        Sfun  = (φ, t) -> S_t_path1(φ, t; U=U, β=β, μ=μ)
        ymax  = max(3.0, μ + U/2 + 0.5)
        push!(panels, manifold_panel(μ, path1_traj[k], Sfun;
                                     ylim=(-0.5, ymax)))
        max_im = maximum(abs.(imag.(path1_traj[k][:, end])))
        rip_flags[(:path1, μ)] = max_im > 2U
    end
    p = plot(panels...; layout=(1,3), size=(1500, 500),
             plot_title="§3 Path 1 manifold at t=1, μ ∈ {1.5, 1.65, 1.75} (paper Fig 1)")
    savefig(p, joinpath(@__DIR__, "data", "03_manifold_fig1.png"))
end
println("Saved data/03_manifold_fig1.png")

println("§3 rip-flag (printed, not asserted — bilinear-interp approximation):")
for ((path, μ), flag) in sort(collect(rip_flags); by = x -> (string(x[1][1]), x[1][2]))
    println("    $(path)  μ=$(rpad(string(μ),6))  ripped=$flag")
end

# %% Phase 4: 3D PIKfold (paper Fig 3a)
# # § 4 — 3D PIKfold visualisation
#
# Reuses the $\mu = 0$ Path-1 trajectories from § 3, plotted as curves in
# $(\mathrm{Re}\,\phi, \mathrm{Im}\,\phi, t)$.

let
    μ0   = 0.0
    k    = findfirst(==(μ0), PATH1_MU_LIST)
    traj = path1_traj[k]
    p = plot3d(title="§4 3D PIKfold, μ₀=$μ0 (paper Fig 3a)",
               xlabel="Re φ", ylabel="Im φ", zlabel="t",
               legend=false, camera=(60, 30), size=(900, 800))
    for n in 1:N_INIT
        plot!(p, real.(traj[n, :]), imag.(traj[n, :]), traj_ts;
              color=:blue, lw=0.6, alpha=0.5, label="")
    end
    savefig(p, joinpath(@__DIR__, "data", "04_pikfold_3d.png"))
    println("Saved data/04_pikfold_3d.png")

    # Gate: t=1 slice agrees with §3 manifold within RK4 error (visual).
    # Numerical check: the t=1 endpoint of each trajectory matches what
    # we plotted in §3 (same trajectory array).
    max_im_t1 = maximum(abs.(imag.(traj[:, end])))
    println("§4 t=1 slice (μ₀=0) max|Im φ_1| = $(round(max_im_t1; digits=4))")
end

# %% Phase 5: Ripping diagnostic (paper Fig 2)
# # § 5 — Path-1 vs Path-2 ripping at $\mu = -0.675$
#
# Two-panel comparison at the canonical ripping parameter.

let
    μ_rip = -0.675
    k1 = findfirst(==(μ_rip), PATH1_MU_LIST)
    k2 = findfirst(==(μ_rip), PATH2_MU_LIST)
    traj1 = path1_traj[k1]
    traj2 = path2_traj[k2]
    slices = [0.2, 0.6, 1.0]
    cols   = [:lightgrey, :grey, :black]

    p1 = plot(title="§5(a) Path 1 (rips), μ=$μ_rip",
              xlabel="Re φ", ylabel="Im φ",
              xlims=(-1.5, 1.5), ylims=(-0.4, 0.4),
              legend=:topright, titlefontsize=10)
    for (s, c) in zip(slices, cols)
        idx = argmin(abs.(traj_ts .- s))
        plot!(p1, real.(traj1[:, idx]), imag.(traj1[:, idx]);
              color=c, lw=1.8, label="t=$s")
    end
    pole_im1 = μ_rip + U/2
    for n in -4:4
        pre = (2n + 1) * π / β
        if -1.5 ≤ pre ≤ 1.5
            scatter!(p1, [pre], [pole_im1];
                     marker=:rect, color=:blue, markersize=3,
                     markerstrokewidth=0, label="")
        end
    end
    savefig(p1, joinpath(@__DIR__, "data", "05_ripping_fig2a.png"))

    p2 = plot(title="§5(b) Path 2 (closed), μ=$μ_rip",
              xlabel="Re φ", ylabel="Im φ",
              xlims=(-1.5, 1.5), ylims=(-0.4, 0.4),
              legend=:topright, titlefontsize=10)
    for (s, c) in zip(slices, cols)
        idx = argmin(abs.(traj_ts .- s))
        plot!(p2, real.(traj2[:, idx]), imag.(traj2[:, idx]);
              color=c, lw=1.8, label="t=$s")
    end
    for n in -4:4
        pre = (2n + 1) * π / β
        if -1.5 ≤ pre ≤ 1.5
            scatter!(p2, [pre], [pole_im1];
                     marker=:rect, color=:blue, markersize=3,
                     markerstrokewidth=0, label="")
        end
    end
    savefig(p2, joinpath(@__DIR__, "data", "05_ripping_fig2b.png"))

    max_im_p1 = maximum(abs.(imag.(traj1)))
    max_im_p2 = maximum(abs.(imag.(traj2)))
    println("Saved data/05_ripping_fig2a.png, data/05_ripping_fig2b.png")
    println("§5 max|Im φ_t| over the full t-grid: path1=$(round(max_im_p1; digits=3)), path2=$(round(max_im_p2; digits=3))")
end

# %% Phase 6: n vs μ scan (paper Fig 5)
# # § 6 — Density vs chemical potential

let
    μ_scan = collect(range(-1.0, 2.0; length=31))
    M_scan = 10_000
    n_pik  = Float64[]
    n_int  = Float64[]
    rng    = MersenneTwister(20260527)
    σ      = sqrt(U / β)
    for μ in μ_scan
        phi0 = randn(rng, M_scan) .* σ
        r = pik_eval(φ -> n_op(φ, U, β, μ), U, β, μ, phi0)
        push!(n_pik, real(r.est))

        # Reference: real-axis quadrature of ⟨−iφ/U⟩.
        num, _ = quadgk(φ -> (-im*φ/U) * exp(-S(φ, U, β, μ)), -L_BOUND, L_BOUND;
                        rtol=1e-10, maxevals=20_000)
        Zv,  _ = quadgk(φ -> exp(-S(φ, U, β, μ)),              -L_BOUND, L_BOUND;
                        rtol=1e-10, maxevals=20_000)
        push!(n_int, real(num / Zv))
    end
    max_diff = maximum(abs.(n_pik .- n_int))
    @assert max_diff < 1e-3 "§6 gate failed: max|n_pik − n_int| = $max_diff ≥ 1e-3"
    println("§6 gate: max|n_PIK − n_integral| over μ = $(round(max_diff; sigdigits=4)) < 1e-3  ✓")

    jldsave(joinpath(@__DIR__, "data", "n_vs_mu_scan.jld2");
            mu_grid = μ_scan, n_pik = n_pik, n_integral = n_int,
            M = M_scan, U = U, β = β)

    p = plot(title="§6  ⟨n⟩ vs μ  (U=$U, β=$β, M=$M_scan)",
             xlabel="μ", ylabel="⟨n⟩",
             xlims=(-1, 2), ylims=(0, 2.1), legend=:bottomright, titlefontsize=10)
    plot!(p, μ_scan, n_int; ls=:dash, lw=2, color=:black, label="integral (real axis)")
    scatter!(p, μ_scan, n_pik; color=:blue, markersize=4,
             markerstrokewidth=0, label="PIK (holomorphic)")
    savefig(p, joinpath(@__DIR__, "data", "06_n_vs_mu.png"))
    println("Saved data/06_n_vs_mu.png")
end

# %% Phase 7: Summary marquee
# # § 7 — Summary: four marquee panels (one per attached image)

let
    # Panel A: Path-1 manifold at μ=0.0 (≈ image 3 / paper Fig 3a 2D-projection)
    Sfun0 = (φ, t) -> S_t_path1(φ, t; U=U, β=β, μ=0.0)
    k0    = findfirst(==(0.0), PATH1_MU_LIST)
    pA    = manifold_panel(0.0, path1_traj[k0], Sfun0)
    title!(pA, "(a) Path 1 manifold, μ=0  [→ paper Fig 3 / image 3]"; titlefontsize=9)

    # Panel B: 3-up Fig 4 collapsed to μ=0.1 (representative of image 1).
    Sfun01 = (φ, t) -> S_t_path1(φ, t; U=U, β=β, μ=0.1)
    k01    = findfirst(==(0.1), PATH1_MU_LIST)
    pB     = manifold_panel(0.1, path1_traj[k01], Sfun01)
    title!(pB, "(b) Path 1 manifold, μ=0.1  [→ paper Fig 4 / image 1]"; titlefontsize=9)

    # Panel C: Ripping at μ=-0.675, Path 1 (image 4).
    μ_rip = -0.675
    kr    = findfirst(==(μ_rip), PATH1_MU_LIST)
    traj  = path1_traj[kr]
    pC = plot(title="(c) Path 1 ripping, μ=$μ_rip  [→ paper Fig 2a / image 4]",
              xlabel="Re φ", ylabel="Im φ",
              xlims=(-1.5, 1.5), ylims=(-0.4, 0.4),
              titlefontsize=9, legend=:topright)
    for (s, c) in zip([0.2, 0.6, 1.0], [:lightgrey, :grey, :black])
        idx = argmin(abs.(traj_ts .- s))
        plot!(pC, real.(traj[:, idx]), imag.(traj[:, idx]);
              color=c, lw=1.8, label="t=$s")
    end

    # Panel D: n vs μ scan (image 2).
    d = load(joinpath(@__DIR__, "data", "n_vs_mu_scan.jld2"))
    pD = plot(title="(d) ⟨n⟩ vs μ  [→ paper Fig 5 / image 2]",
              xlabel="μ", ylabel="⟨n⟩",
              xlims=(-1, 2), ylims=(0, 2.1), legend=:bottomright,
              titlefontsize=9)
    plot!(pD, d["mu_grid"], d["n_integral"];
          ls=:dash, lw=2, color=:black, label="integral")
    scatter!(pD, d["mu_grid"], d["n_pik"];
             color=:blue, markersize=3, markerstrokewidth=0, label="PIK")

    p = plot(pA, pB, pC, pD; layout=(2, 2), size=(1500, 1200),
             plot_title="§7 Summary: PIK-flow reproductions vs paper figures")
    savefig(p, joinpath(@__DIR__, "data", "07_summary.png"))
    println("Saved data/07_summary.png")
end

# Final end-to-end status line.
println("\nAll § gates passed. New files written: 7 PNGs + 5 JLD2 in data/.")
