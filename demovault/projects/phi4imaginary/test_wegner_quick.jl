using QuadGK, JLD2, Plots

# ── Load parameters ───────────────────────────────────────────────────────────
data = load(joinpath(@__DIR__, "data", "action_check.jld2"))
const M2  = data["M2"]
const LAM = data["LAM"]
println("Loaded M2=$(M2), LAM=$(LAM)")

# ── Action and helpers (same as 01) ──────────────────────────────────────────
S_hat(φ, t) = M2 * exp(-im*π/2*t)/2 * φ^2 + LAM * exp(-im*π/2*t)/4 * φ^4

function integration_bound(t)
    c = cos(π * t / 2)
    c > 0.01 && return clamp(sqrt(56.0 / c), 8.0, 50.0)
    return 50.0
end

function action_mean(t; L=integration_bound(t))
    num, _ = quadgk(φ -> S_hat(φ, t) * exp(-S_hat(φ, t)), -L, L, rtol=1e-8)
    den, _ = quadgk(φ -> exp(-S_hat(φ, t)), -L, L, rtol=1e-8)
    return num / den
end

function phi_dot(φ_real::Real, t, Sbar; L=integration_bound(t))
    integrand(x) = (S_hat(x, t) - Sbar) * exp(-S_hat(x, t))
    integral, _ = quadgk(integrand, φ_real, L, rtol=1e-8)  # L approximates +∞
    return (im*π/2) * exp(S_hat(φ_real, t)) * integral
end

# Quick test with smaller grid
const N_phi = 50
const N_t   = 5
phi_grid = range(-5.0, 5.0, length=N_phi)
t_grid   = range(0.0, 0.99, length=N_t)

phi_dot_grid = zeros(ComplexF64, N_t, N_phi)

println("Computing velocity field on $(N_t) × $(N_phi) grid...")
for (i, t) in enumerate(t_grid)
    Sbar = action_mean(t)
    for (j, φ) in enumerate(phi_grid)
        phi_dot_grid[i, j] = phi_dot(φ, t, Sbar)
    end
    println("  t-slice $i / $N_t done")
end
println("Done.")

# Test point (at φ=0, t=0)
idx_phi_zero = argmin(abs.(collect(phi_grid) .- 0.0))
idx_t_zero = 1
test_val = phi_dot_grid[idx_t_zero, idx_phi_zero]
println("Test: phi_dot(0, 0) = $(round(test_val, digits=6))")
println("Expected: close to 0 (symmetry)")
println("Status: $(abs(test_val) < 1e-6 ? "PASS" : "CONCERN")")
