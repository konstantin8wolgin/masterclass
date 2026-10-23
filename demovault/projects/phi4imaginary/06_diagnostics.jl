using QuadGK, JLD2, Plots

include("reg.jl")
using .Regulators

data = load("data/action_check.jld2")
const M2  = data["M2"]
const LAM = data["LAM"]

const REG_METHOD = :tikhonov
const EPS_REG    = 1e-3

S_hat(φ, t)     = Regulators.S_hat_bare(φ, t; m2=M2, lam=LAM)
S_hat_reg(φ, t) = Regulators.S_reg(φ, t; method=REG_METHOD, ε=EPS_REG, m2=M2, lam=LAM)

weg  = load("data/wegner.jld2")
obs  = load("data/observables.jld2")
phi_grid    = weg["phi_grid"]
t_grid_w    = weg["t_grid"]
phi_dot_grid= weg["phi_dot"]

t_obs       = obs["t_grid"]
phi2_arr    = obs["phi2"]
phi4_arr    = obs["phi4"]
Smean_arr   = obs["Smean"]

# ── 1. Action landscape ──────────────────────────────────────────────────────
phi_h = range(-5.0, 5.0, length=160)
t_h   = range(0.0, 1.0, length=120)
S_re = [real(S_hat(φ, t)) for t in t_h, φ in phi_h]
S_im = [imag(S_hat(φ, t)) for t in t_h, φ in phi_h]

h1 = heatmap(phi_h, t_h, S_re, xlabel="φ", ylabel="t",
             title="Re[S_hat(φ,t)]", color=:RdBu, clims=(-12,12))
h2 = heatmap(phi_h, t_h, S_im, xlabel="φ", ylabel="t",
             title="Im[S_hat(φ,t)]", color=:RdBu, clims=(-12,12))
savefig(plot(h1, h2, layout=(1,2), size=(1100,450)),
        "data/06_action_landscape.png")
println("Saved data/06_action_landscape.png")

# ── 2. Weight landscape |exp(-S)| and arg(exp(-S)) ───────────────────────────
W_abs = [abs(exp(-S_hat(φ, t))) for t in t_h, φ in phi_h]
W_arg = [angle(exp(-S_hat(φ, t))) for t in t_h, φ in phi_h]

h3 = heatmap(phi_h, t_h, log10.(W_abs .+ 1e-30),
             xlabel="φ", ylabel="t",
             title="log₁₀ |exp(-S_hat)|", color=:viridis)
h4 = heatmap(phi_h, t_h, W_arg, xlabel="φ", ylabel="t",
             title="arg(exp(-S_hat))", color=:twilight, clims=(-π,π))
savefig(plot(h3, h4, layout=(1,2), size=(1100,450)),
        "data/06_weight_landscape.png")
println("Saved data/06_weight_landscape.png")

# ── 3. Quiver / velocity flow field in (φ, t) ────────────────────────────────
phi_q = range(-4.5, 4.5, length=18)
t_q   = range(0.05, 0.95, length=14)
xs = Float64[]; ys = Float64[]; us = Float64[]; vs = Float64[]
# Use stored phi_dot_grid: scale arrows by 1/max for visibility
maxv = maximum(abs.(phi_dot_grid))
scale = 0.25 / max(maxv, 1e-12)
for t in t_q, φ in phi_q
    i = argmin(abs.(t_grid_w .- t))
    j = argmin(abs.(phi_grid .- φ))
    v = phi_dot_grid[i, j]
    push!(xs, φ); push!(ys, t)
    push!(us, real(v)*scale); push!(vs, imag(v)*scale)
end
pq = quiver(xs, ys, quiver=(us, vs),
            xlabel="Re[φ]", ylabel="t",
            title="Wegner velocity (Re[φ̇], Im[φ̇]) on (φ,t) grid",
            color=:darkblue, linewidth=0.8)
savefig(pq, "data/06_velocity_field_quiver.png")
println("Saved data/06_velocity_field_quiver.png")

# ── 4. Thimble convergence: Im[end-point] vs phi0 ────────────────────────────
flow = load("data/flow.jld2")
traj = flow["trajectories"]
phi0_values = flow["phi0_values"]
t_steps     = flow["t_steps"]

# Pick a handful of times along the flow
t_pick = [0.25, 0.5, 0.75, t_steps[end]]
pTC = plot(title="Thimble convergence: Im[φ(t)] vs initial real φ₀",
           xlabel="φ₀ (initial real)", ylabel="Im[φ(t)]", legend=:topleft)
cs = [:blue, :green, :orange, :red]
for (tv, c) in zip(t_pick, cs)
    k = argmin(abs.(t_steps .- tv))
    plot!(pTC, phi0_values, imag.(traj[:, k]),
          lw=2, color=c, label="t=$(round(t_steps[k],digits=2))")
end
# free-theory analytic at t_end: Im = phi0 * sin(π t/2)
plot!(pTC, phi0_values, phi0_values .* sin(π*t_steps[end]/2),
      lw=1.5, color=:gray, linestyle=:dash, label="free thimble at t_end")
savefig(pTC, "data/06_thimble_convergence.png")
println("Saved data/06_thimble_convergence.png")

# ── 5. Extra observables: ⟨φ⁴⟩, ⟨S_hat⟩, Binder-like cumulant ────────────────
p5a = plot(t_obs, real.(phi4_arr), label="Re⟨φ⁴⟩", lw=2, color=:blue,
           xlabel="t", ylabel="⟨φ⁴⟩", title="Quartic observable")
plot!(p5a, t_obs, imag.(phi4_arr), label="Im⟨φ⁴⟩", lw=2, color=:blue, linestyle=:dash)

# ⟨S_hat⟩(t) — expected ∝ exp(-iπt/2) · ⟨S_hat⟩(0)
S0 = Smean_arr[1]
analytic_S = @. exp(-im*π/2*t_obs) * S0
p5b = plot(t_obs, real.(Smean_arr), label="Re⟨S_hat⟩ num", lw=2, color=:darkgreen,
           xlabel="t", ylabel="⟨S_hat⟩", title="⟨S_hat⟩(t): numeric vs analytic rotation")
plot!(p5b, t_obs, imag.(Smean_arr), label="Im⟨S_hat⟩ num", lw=2, color=:darkgreen, linestyle=:dash)
plot!(p5b, t_obs, real.(analytic_S), label="Re analytic", lw=1.5, color=:black, linestyle=:dot)
plot!(p5b, t_obs, imag.(analytic_S), label="Im analytic", lw=1.5, color=:gray,  linestyle=:dot)

# Binder-like cumulant Re[⟨φ⁴⟩ - 3⟨φ²⟩²]
binder = phi4_arr .- 3 .* phi2_arr.^2
p5c = plot(t_obs, real.(binder), label="Re", lw=2, color=:purple,
           xlabel="t", ylabel="⟨φ⁴⟩ − 3⟨φ²⟩²", title="Connected 4-pt cumulant")
plot!(p5c, t_obs, imag.(binder), label="Im", lw=2, color=:purple, linestyle=:dash)

savefig(plot(p5a, p5b, p5c, layout=(1,3), size=(1500, 450)),
        "data/06_observables_extra.png")
println("Saved data/06_observables_extra.png")

# ── 6. ε-extrapolation: Re/Im ⟨φ²⟩ vs t for ε ∈ {1e-2, 1e-3, 1e-4} ───────────
function phi2_for_eps(t, ε)
    L = Regulators.L_smooth(t; ε=ε)
    Sr(φ) = Regulators.S_reg(φ, t; method=:tikhonov, ε=ε, m2=M2, lam=LAM)
    w(φ)  = exp(-Sr(φ))
    n, _  = quadgk(φ -> φ^2 * w(φ), -L, L, rtol=1e-8)
    d, _  = quadgk(w,                -L, L, rtol=1e-8)
    return n/d
end

t_eps_grid = range(0.0, 1.0, length=60)
eps_list   = [1e-2, 1e-3, 1e-4]
eps_colors = [:red, :blue, :green]

pE_re = plot(title="Re⟨φ²⟩(t) — ε extrapolation", xlabel="t", ylabel="Re⟨φ²⟩", legend=:topright)
pE_im = plot(title="Im⟨φ²⟩(t) — ε extrapolation", xlabel="t", ylabel="Im⟨φ²⟩", legend=:topright)
for (ε, c) in zip(eps_list, eps_colors)
    vals = [phi2_for_eps(t, ε) for t in t_eps_grid]
    plot!(pE_re, t_eps_grid, real.(vals), lw=2, color=c, label="ε=$(ε)")
    plot!(pE_im, t_eps_grid, imag.(vals), lw=2, color=c, label="ε=$(ε)")
end
# free-theory overlay
plot!(pE_re, t_eps_grid, [real(exp(im*π/2*t)/M2) for t in t_eps_grid],
      lw=1.5, color=:black, linestyle=:dot, label="free")
plot!(pE_im, t_eps_grid, [imag(exp(im*π/2*t)/M2) for t in t_eps_grid],
      lw=1.5, color=:black, linestyle=:dot, label="free")
savefig(plot(pE_re, pE_im, layout=(1,2), size=(1100, 450)),
        "data/06_epsilon_extrapolation.png")
println("Saved data/06_epsilon_extrapolation.png")

println("\nAll diagnostics figures produced.")
