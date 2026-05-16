using QuadGK, JLD2, Plots

# ── Load parameters ───────────────────────────────────────────────────────────
data = load(joinpath(@__DIR__, "data", "action_check.jld2"))
const M2_0  = data["M2_0"]
const M2_1  = data["M2_1"]
const LAM_1 = data["LAM_1"]

# ── Action path (same as 01) ──────────────────────────────────────────────────
m2(t)       = M2_0 + t * (M2_1 - M2_0)
lam(t)      = t * LAM_1
S_hat(φ, t) = m2(t)/2 * φ^2 + lam(t)/4 * φ^4
dS_dt(φ, t) = -3/2 * φ^2 + 1/8 * φ^4

function integration_bound(t)
    lam_t = lam(t)
    lam_t > 0.01 && return clamp(sqrt(sqrt(96.0 / lam_t)), 8.0, 30.0)
    m2_t = m2(t)
    m2_t > 0.01 && return clamp(sqrt(56.0 / m2_t), 8.0, 50.0)
    return 50.0
end

# ── Mean ⟨dŜ/dt⟩ under exp(-Ŝ_t) ────────────────────────────────────────────
function action_deriv_mean(t; L=integration_bound(t))
    num, _ = quadgk(φ -> dS_dt(φ,t) * exp(-S_hat(φ,t)), -L, L; rtol=1e-8, maxevals=5_000)
    den, _ = quadgk(φ -> exp(-S_hat(φ,t)),               -L, L; rtol=1e-8, maxevals=5_000)
    return num / den
end

# ── Velocity field ─────────────────────────────────────────────────────────────
# For φ < 0, use ∫_φ^∞ = −∫_{-∞}^φ to avoid catastrophic cancellation:
# integrating from φ to L spans the full distribution and nearly cancels.
function phi_dot(φ::Real, t, Sbar; L=integration_bound(t))
    integrand(x) = (dS_dt(x,t) - Sbar) * exp(-S_hat(x,t))
    if φ >= 0
        # ∫_{-L}^φ ≈ −∫_φ^L; right-tail integral avoids catastrophic cancellation
        integral, _ = quadgk(integrand, φ, L; rtol=1e-8, maxevals=5_000)
        return -exp(S_hat(φ,t)) * integral
    else
        # Left-tail integral ∫_{-L}^φ is naturally small
        integral, _ = quadgk(integrand, -L, φ; rtol=1e-8, maxevals=5_000)
        return exp(S_hat(φ,t)) * integral
    end
end

# ── Grid ──────────────────────────────────────────────────────────────────────
const N_phi = 120
const N_t   = 50
phi_grid = range(-6.0, 6.0, length=N_phi)
t_grid   = range(0.0,  1.0, length=N_t)

# ── Compute velocity field on grid (threaded over φ) ─────────────────────────
phi_dot_grid = zeros(Float64, N_t, N_phi)

println("Computing velocity field on $(N_t)×$(N_phi) grid ($(Threads.nthreads()) thread(s))...")
for (i, t) in enumerate(t_grid)
    Sbar = action_deriv_mean(t)
    phi_vec = collect(phi_grid)
    results = zeros(Float64, N_phi)
    Threads.@threads for j in 1:N_phi
        results[j] = phi_dot(phi_vec[j], t, Sbar)
    end
    phi_dot_grid[i, :] = results
    i % 10 == 0 && println("  t-slice $i / $N_t  (t=$(round(t, digits=3))) done")
end
println("Done.")

# ── Save ──────────────────────────────────────────────────────────────────────
jldsave(joinpath(@__DIR__, "data", "wegner.jld2");
    phi_grid = collect(phi_grid),
    t_grid   = collect(t_grid),
    phi_dot  = phi_dot_grid)
println("Saved to data/wegner.jld2")

# ── Plot: velocity field at 5 t-values ────────────────────────────────────────
plot_ts = [0.0, 0.25, 0.5, 0.75, 1.0]
colors  = [:blue, :green, :orange, :red, :purple]

phi_vec_plot = collect(phi_grid)
t_vec_plot   = collect(t_grid)

p = plot(title="Wegner velocity field φ̇_t(φ)", xlabel="φ", ylabel="φ̇_t(φ)",
         legend=:topleft, size=(700, 450))

for (tv, col) in zip(plot_ts, colors)
    idx = argmin(abs.(t_vec_plot .- tv))
    t_actual = t_grid[idx]
    plot!(p, phi_vec_plot, phi_dot_grid[idx, :],
          label="t=$(round(t_actual, digits=2))", color=col, linewidth=1.5)
end

hline!(p, [0.0], color=:black, linestyle=:dash, label="", linewidth=0.5)
savefig(p, joinpath(@__DIR__, "data", "02_velocity_field.png"))
println("Plot saved to data/02_velocity_field.png")
