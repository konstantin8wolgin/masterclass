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
L0 = domain_limit(0.0)
phi2_real_direct, _ = quadgk(φ -> φ^2 * exp(-S_hat(φ, 0.0)), -L0, L0, rtol=1e-10)
norm_direct,      _ = quadgk(φ ->       exp(-S_hat(φ, 0.0)), -L0, L0, rtol=1e-10)
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
plot!(p1, t_vec, imag.(phi2_arr),
    label="Im⟨φ²⟩", lw=2, color=:blue, linestyle=:dash)
plot!(p1, t_vec, imag.(absphi_arr),
    label="Im⟨|φ|⟩", lw=2, color=:red, linestyle=:dash)
savefig(p1, "data/04_observables.png")

p2 = plot(t_vec, abs.(avgsign_arr),
    xlabel="t", ylabel="|avg sign|",
    title="Sign problem severity: |⟨e^{-i Im S}⟩_Q|",
    lw=2, color=:black, ylims=(0, 1.05), legend=false)
hline!(p2, [1.0], linestyle=:dash, color=:gray)
savefig(p2, "data/04_sign_severity.png")

println("Plots saved to data/04_observables.png and data/04_sign_severity.png")
