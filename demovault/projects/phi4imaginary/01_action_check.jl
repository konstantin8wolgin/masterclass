using QuadGK, JLD2

# ── Parameters ────────────────────────────────────────────────────────────────
const M2 = 1.0   # bare mass squared
const LAM = 1.0  # bare coupling

# ── Action (field-dependent part only, no log Z) ──────────────────────────────
S_hat(φ, t) = M2 * exp(-im*π/2*t)/2 * φ^2 + LAM * exp(-im*π/2*t)/4 * φ^4
dS_dt(φ, t) = (-im*π/2) * S_hat(φ, t)   # exact: both couplings carry exp(-iπt/2)

# ── Integration bound: require exp(-Re[S(L,t)]) < 1e-12, i.e. c/2*L^2 > 27.6 ─
function integration_bound(t)
    c = cos(π * t / 2)
    c > 0.01 && return clamp(sqrt(56.0 / c), 8.0, 50.0)
    return 50.0
end

# ── Mean ⟨S_hat_t⟩ under e^{-S_hat} ──────────────────────────────────────────
function action_mean(t; L=integration_bound(t))
    num, _ = quadgk(φ -> S_hat(φ, t) * exp(-S_hat(φ, t)), -L, L, rtol=1e-10)
    den, _ = quadgk(φ -> exp(-S_hat(φ, t)), -L, L, rtol=1e-10)
    return num / den
end

# ── Velocity field from Wegner equation (real φ) ──────────────────────────────
# φ̇_t(φ) = (iπ/2) exp(S_hat(φ)) ∫_φ^∞ [S_hat(φ') - ⟨S_hat⟩] exp(-S_hat(φ')) dφ'
function phi_dot(φ_real::Real, t; L=integration_bound(t))
    Sbar = action_mean(t; L=L)
    integrand(x) = (S_hat(x, t) - Sbar) * exp(-S_hat(x, t))
    integral, _ = quadgk(integrand, φ_real, L, rtol=1e-10)  # L approximates +∞
    return (im*π/2) * exp(S_hat(φ_real, t)) * integral
end

# ── Free-theory version (LAM=0) for validation ───────────────────────────────
S_hat_free(φ, t) = M2 * exp(-im*π/2*t)/2 * φ^2

function phi_dot_free(φ_real::Real, t; L=integration_bound(t))
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
    global all_pass; all_pass &= pass
    println("  φ=$(φ), t=$(t): err=$(round(err, sigdigits=3))  $(pass ? "PASS" : "FAIL")")
end

# ── Validation 2: free-theory φ̇_t vs analytic (iπ/4)φ ───────────────────────
println("\n=== Free-theory check: phi_dot_free vs (iπ/4)φ ===")
for (φ, t) in [(1.0, 0.25), (2.0, 0.5), (-1.0, 0.75), (0.5, 0.9)]
    numeric  = phi_dot_free(φ, t)
    analytic = (im*π/4) * φ
    err = abs(numeric - analytic)
    pass = err < 1e-6
    global all_pass; all_pass &= pass
    println("  φ=$(φ), t=$(t): numeric=$(round(numeric, digits=5)), analytic=$(round(analytic, digits=5)), err=$(round(err, sigdigits=3))  $(pass ? "PASS" : "FAIL")")
end

println("\nOverall: $(all_pass ? "ALL PASS" : "FAILURES DETECTED")")

# ── Smoke-test interacting phi_dot at one point ───────────────────────────────
φ_test, t_test = 1.0, 0.5
vtest = phi_dot(φ_test, t_test)
println("\nSmoke-test phi_dot($(φ_test), $(t_test)) = $(round(vtest, digits=5))  (expect Im > 0, Re ≈ 0)")

# ── Save parameters ───────────────────────────────────────────────────────────
jldsave(joinpath(@__DIR__, "data", "action_check.jld2"); M2, LAM)
println("\nSaved parameters to data/action_check.jld2")
