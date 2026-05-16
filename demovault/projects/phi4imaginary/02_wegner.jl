using QuadGK, JLD2, Plots

# ── Load parameters ───────────────────────────────────────────────────────────
data = load(joinpath(@__DIR__, "data", "action_check.jld2"))
const M2  = data["M2"]
const LAM = data["LAM"]

# ── Action and helpers (same as 01) ──────────────────────────────────────────
S_hat(φ, t) = M2 * exp(-im*π/2*t)/2 * φ^2 + LAM * exp(-im*π/2*t)/4 * φ^4

function integration_bound(t)
    c = cos(π * t / 2)
    c > 0.01 && return clamp(sqrt(56.0 / c), 8.0, 50.0)
    return 50.0
end

function action_mean(t; L=integration_bound(t))
    num, _ = quadgk(φ -> S_hat(φ, t) * exp(-S_hat(φ, t)), -L, L;
                    rtol=1e-6, maxevals=5_000)
    den, _ = quadgk(φ -> exp(-S_hat(φ, t)), -L, L;
                    rtol=1e-6, maxevals=5_000)
    return num / den
end

function phi_dot(φ_real::Real, t, Sbar; L=integration_bound(t))
    integrand(x) = (S_hat(x, t) - Sbar) * exp(-S_hat(x, t))
    integral, _ = quadgk(integrand, φ_real, L;
                         rtol=1e-6, maxevals=5_000)
    return (im*π/2) * exp(S_hat(φ_real, t)) * integral
end

# ── Grid setup ────────────────────────────────────────────────────────────────
# Reduced grid for tractable runtime; sufficient for RK4 interpolation in 03
const N_phi = 120
const N_t   = 25
phi_grid = range(-5.0, 5.0, length=N_phi)
t_grid   = range(0.0, 0.99, length=N_t)

# ── Compute phi_dot on grid (multi-threaded over phi) ─────────────────────────
phi_dot_grid = zeros(ComplexF64, N_t, N_phi)

println("Computing velocity field on $(N_t) × $(N_phi) grid ($(Threads.nthreads()) thread(s))...")
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
    phi_dot  = phi_dot_grid)
println("Saved to data/wegner.jld2")

# ── Plot: Re and Im of phi_dot_t vs phi for 5 t-values ───────────────────────
plot_t_vals = [0.0, 0.25, 0.5, 0.75, 0.99]
colors = [:blue, :green, :orange, :red, :purple]

phi_vec_plot = collect(phi_grid)
t_vec_plot   = collect(t_grid)

p1 = plot(title="Re[φ̇_t(φ)] vs φ", xlabel="φ", ylabel="Re[φ̇_t]", legend=:topright)
p2 = plot(title="Im[φ̇_t(φ)] vs φ", xlabel="φ", ylabel="Im[φ̇_t]", legend=:topright)

for (tv, col) in zip(plot_t_vals, colors)
    idx = argmin(abs.(t_vec_plot .- tv))
    t_actual = t_grid[idx]
    phi_dot_row = phi_dot_grid[idx, :]
    plot!(p1, phi_vec_plot, real.(phi_dot_row), label="t=$(round(t_actual,digits=2))", color=col)
    plot!(p2, phi_vec_plot, imag.(phi_dot_row), label="t=$(round(t_actual,digits=2))", color=col)
end

savefig(plot(p1, p2, layout=(1,2), size=(900,400)), joinpath(@__DIR__, "data", "02_velocity_field.png"))
println("Plot saved to data/02_velocity_field.png")
