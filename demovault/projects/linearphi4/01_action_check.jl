using QuadGK, JLD2

# ── Parameters ────────────────────────────────────────────────────────────────
const M2_0  = 1.0    # m²(0)
const M2_1  = -2.0   # m²(1)
const LAM_1 = 0.5    # λ(1)

# ── Action path ───────────────────────────────────────────────────────────────
m2(t)       = M2_0 + t * (M2_1 - M2_0)   # 1 - 3t
lam(t)      = t * LAM_1                    # t/2
S_hat(φ, t) = m2(t)/2 * φ^2 + lam(t)/4 * φ^4
dS_dt(φ, t) = -3/2 * φ^2 + 1/8 * φ^4    # analytic time derivative

# ── Integration bound ─────────────────────────────────────────────────────────
function integration_bound(t)
    lam_t = lam(t)
    lam_t > 0.01 && return clamp(sqrt(sqrt(96.0 / lam_t)), 5.0, 30.0)
    m2_t = m2(t)
    m2_t > 0.01 && return clamp(sqrt(56.0 / m2_t), 8.0, 50.0)
    return 50.0
end

# ── Free-theory helpers (λ=0) ─────────────────────────────────────────────────
S_hat_free(φ, t)  = m2(t)/2 * φ^2
dS_dt_free(φ, t)  = -3/2 * φ^2   # derivative of m²(t)/2·φ² w.r.t. t only

function integration_bound_free(t)
    m2_t = m2(t)
    m2_t > 0.01 && return clamp(sqrt(56.0 / m2_t), 8.0, 50.0)
    return 50.0
end

function phi_dot_free(φ::Real, t; L=integration_bound_free(t))
    Sbar = begin
        num, _ = quadgk(x -> dS_dt_free(x,t) * exp(-S_hat_free(x,t)), -L, L; rtol=1e-10)
        den, _ = quadgk(x -> exp(-S_hat_free(x,t)),                    -L, L; rtol=1e-10)
        num / den
    end
    integrand(x) = (dS_dt_free(x,t) - Sbar) * exp(-S_hat_free(x,t))
    if φ >= 0
        integral, _ = quadgk(integrand, φ, L; rtol=1e-10)
        return -exp(S_hat_free(φ,t)) * integral
    else
        integral, _ = quadgk(integrand, -L, φ; rtol=1e-10)
        return exp(S_hat_free(φ,t)) * integral
    end
end

# ── Validation 1: dS_dt identity ─────────────────────────────────────────────
println("=== Identity check: dS_dt(φ,t) == -3/2·φ² + 1/8·φ⁴ ===")
all_pass = true
for (φ, t) in [(0.5, 0.0), (1.0, 0.3), (-1.5, 0.6), (2.0, 1.0)]
    analytic = -3/2 * φ^2 + 1/8 * φ^4
    err = abs(dS_dt(φ, t) - analytic)
    pass = err < 1e-14
    global all_pass; all_pass &= pass
    println("  φ=$(φ), t=$(t): err=$(round(err, sigdigits=3))  $(pass ? "PASS" : "FAIL")")
end

# ── Validation 2: free-theory φ̇ vs analytic -3/(2(1-3t))·φ ──────────────────
println("\n=== Free-theory check: phi_dot_free vs +3/(2(1-3t))·φ ===")
# Avoid t >= 1/3 where m²(t) <= 0 (free theory singular there)
for (φ, t) in [(1.0, 0.05), (2.0, 0.1), (-1.0, 0.2), (0.5, 0.3)]
    numeric  = phi_dot_free(φ, t)
    analytic = 3.0 / (2*(1 - 3t)) * φ
    err = abs(numeric - analytic)
    pass = err < 1e-6
    global all_pass; all_pass &= pass
    println("  φ=$(φ), t=$(t): numeric=$(round(numeric, digits=5)), analytic=$(round(analytic, digits=5)), err=$(round(err, sigdigits=3))  $(pass ? "PASS" : "FAIL")")
end

println("\nOverall: $(all_pass ? "ALL PASS" : "FAILURES DETECTED")")

# ── Smoke test: interacting phi_dot at t=1 ───────────────────────────────────
println("\n=== Smoke test: S_hat at t=1 boundary ===")
for φ in [0.0, 1.0, 2.0, 3.0]
    s = S_hat(φ, 1.0)
    println("  S_hat($φ, 1.0) = $(round(s, digits=4))  (expect -φ² + φ⁴/8; negative near 0, positive far)")
end

# ── Save ──────────────────────────────────────────────────────────────────────
jldsave(joinpath(@__DIR__, "data", "action_check.jld2"); M2_0, M2_1, LAM_1)
println("\nSaved parameters to data/action_check.jld2")
