using QuadGK, JLD2, Plots

include("reg.jl")
using .Regulators

# ── Load parameters ───────────────────────────────────────────────────────────
data = load(joinpath(@__DIR__, "data", "action_check.jld2"))
const M2  = data["M2"]
const LAM = data["LAM"]

const REG_METHOD = :tikhonov
const EPS_REG    = 1e-3

# ── Action and helpers (regularised) ─────────────────────────────────────────
S_hat(φ, t)     = Regulators.S_hat_bare(φ, t; m2=M2, lam=LAM)
S_hat_reg(φ, t) = Regulators.S_reg(φ, t; method=REG_METHOD, ε=EPS_REG, m2=M2, lam=LAM)

integration_bound(t) = Regulators.L_smooth(t; ε=EPS_REG)

function action_mean(t; L=integration_bound(t))
    num, _ = quadgk(φ -> S_hat_reg(φ, t) * exp(-S_hat_reg(φ, t)), -L, L;
                    rtol=1e-6, maxevals=5_000)
    den, _ = quadgk(φ -> exp(-S_hat_reg(φ, t)), -L, L;
                    rtol=1e-6, maxevals=5_000)
    return num / den
end

function phi_dot(φ_real::Real, t, Sbar; L=integration_bound(t))
    integrand(x) = (S_hat_reg(x, t) - Sbar) * exp(-S_hat_reg(x, t))
    integral, _ = quadgk(integrand, φ_real, L;
                         rtol=1e-6, maxevals=5_000)
    return (im*π/2) * exp(S_hat_reg(φ_real, t)) * integral
end

# ── Grid setup ────────────────────────────────────────────────────────────────
const N_phi = 120
const N_t   = 25
phi_grid = range(-5.0, 5.0, length=N_phi)
t_grid   = range(0.0, 1.0, length=N_t)

# ── Compute phi_dot on grid (multi-threaded over phi) ─────────────────────────
phi_dot_grid = zeros(ComplexF64, N_t, N_phi)

println("Computing velocity field on $(N_t) × $(N_phi) grid ($(Threads.nthreads()) thread(s))...")
println("Regulator: $(REG_METHOD), ε = $(EPS_REG)")
for (i, t) in enumerate(t_grid)
    Sbar = action_mean(t)
    phi_vec = collect(phi_grid)
    results = zeros(ComplexF64, N_phi)
    Threads.@threads for j in 1:N_phi
        results[j] = phi_dot(phi_vec[j], t, Sbar)
    end
    phi_dot_grid[i, :] = results
    println("  t-slice $i / $N_t  (t=$(round(t, digits=3))) done")
end
println("Done.")

# ── Save ──────────────────────────────────────────────────────────────────────
jldsave(joinpath(@__DIR__, "data", "wegner.jld2");
    phi_grid = collect(phi_grid),
    t_grid   = collect(t_grid),
    phi_dot  = phi_dot_grid,
    reg_method = String(REG_METHOD),
    eps_reg    = EPS_REG)
println("Saved to data/wegner.jld2")

# ── Plot: Re and Im of phi_dot_t vs phi for 5 t-values (linear) ──────────────
plot_t_vals = [0.0, 0.25, 0.5, 0.75, 1.0]
colors = [:blue, :green, :orange, :red, :purple]

phi_vec_plot = collect(phi_grid)
t_vec_plot   = collect(t_grid)

p1 = plot(title="Re[φ̇_t(φ)] vs φ", xlabel="φ", ylabel="Re[φ̇_t]",
          legend=:topright, ylims=(-3, 3), xlims=(-4, 4))
p2 = plot(title="Im[φ̇_t(φ)] vs φ", xlabel="φ", ylabel="Im[φ̇_t]",
          legend=:topright, ylims=(-3, 3), xlims=(-4, 4))

for (tv, col) in zip(plot_t_vals, colors)
    idx = argmin(abs.(t_vec_plot .- tv))
    t_actual = t_grid[idx]
    phi_dot_row = phi_dot_grid[idx, :]
    plot!(p1, phi_vec_plot, real.(phi_dot_row), label="t=$(round(t_actual,digits=2))", color=col)
    plot!(p2, phi_vec_plot, imag.(phi_dot_row), label="t=$(round(t_actual,digits=2))", color=col)
end

savefig(plot(p1, p2, layout=(1,2), size=(900,400)), joinpath(@__DIR__, "data", "02_velocity_field.png"))
println("Plot saved to data/02_velocity_field.png")

# ── Log-scale velocity-field plots ───────────────────────────────────────────
# Avoid log(0) by clamping tiny values
clamp_pos(x) = max(abs(x), 1e-30)

pL1 = plot(title="|Re[φ̇_t(φ)]| vs φ (semilog-y)", xlabel="φ",
           ylabel="|Re[φ̇_t]|", yscale=:log10, legend=:topright)
pL2 = plot(title="|Im[φ̇_t(φ)]| vs φ (semilog-y)", xlabel="φ",
           ylabel="|Im[φ̇_t]|", yscale=:log10, legend=:topright)
pLL1 = plot(title="|Re[φ̇_t]| vs |φ| (log-log, φ>0)", xlabel="|φ|",
            ylabel="|Re[φ̇_t]|", xscale=:log10, yscale=:log10, legend=:topright)
pLL2 = plot(title="|Im[φ̇_t]| vs |φ| (log-log, φ>0)", xlabel="|φ|",
            ylabel="|Im[φ̇_t]|", xscale=:log10, yscale=:log10, legend=:topright)

pos_mask = phi_vec_plot .> 0
phi_pos  = phi_vec_plot[pos_mask]

for (tv, col) in zip(plot_t_vals, colors)
    idx = argmin(abs.(t_vec_plot .- tv))
    t_actual = t_grid[idx]
    row = phi_dot_grid[idx, :]
    plot!(pL1, phi_vec_plot, clamp_pos.(real.(row)), label="t=$(round(t_actual,digits=2))", color=col)
    plot!(pL2, phi_vec_plot, clamp_pos.(imag.(row)), label="t=$(round(t_actual,digits=2))", color=col)
    plot!(pLL1, phi_pos, clamp_pos.(real.(row[pos_mask])), label="t=$(round(t_actual,digits=2))", color=col)
    plot!(pLL2, phi_pos, clamp_pos.(imag.(row[pos_mask])), label="t=$(round(t_actual,digits=2))", color=col)
end

savefig(plot(pL1, pL2, pLL1, pLL2, layout=(2,2), size=(1100,800)),
        joinpath(@__DIR__, "data", "02_velocity_field_log.png"))
println("Plot saved to data/02_velocity_field_log.png")
