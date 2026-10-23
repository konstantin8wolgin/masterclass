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

# Smooth domain limit — no abrupt switch at t→1
domain_limit(t) = Regulators.L_smooth(t; ε=EPS_REG)

function compute_observables(t)
    L = domain_limit(t)
    # Use the regularised action everywhere so the t→1 limit is smooth.
    phase(φ)    = exp(-im * imag(S_hat_reg(φ, t)))
    quenched(φ) = exp(-real(S_hat_reg(φ, t)))
    weight(φ)   = phase(φ) * quenched(φ)

    norm_Q,  _ = quadgk(quenched,                  -L, L, rtol=1e-8)
    denom,   _ = quadgk(weight,                    -L, L, rtol=1e-8)
    num_phi2,_ = quadgk(φ -> φ^2   * weight(φ),    -L, L, rtol=1e-8)
    num_phi4,_ = quadgk(φ -> φ^4   * weight(φ),    -L, L, rtol=1e-8)
    num_abs, _ = quadgk(φ -> abs(φ) * weight(φ),   -L, L, rtol=1e-8)
    num_S,   _ = quadgk(φ -> S_hat(φ,t) * weight(φ), -L, L, rtol=1e-8)

    phi2     = num_phi2 / denom
    phi4     = num_phi4 / denom
    absphi   = num_abs  / denom
    Smean    = num_S    / denom
    avg_sign = denom    / norm_Q
    return phi2, phi4, absphi, Smean, avg_sign
end

# ── Compute over t grid ───────────────────────────────────────────────────────
const N_t = 100
# t exactly = 1 is the genuine Minkowski sign singularity (cos(πt/2)=0). The
# Tikhonov regulator gives smooth values for t up to ~0.99; t=1.0 itself
# is reported via free-theory analytic and via the thimble regulator in 07.
t_grid = range(0.0, 0.99, length=N_t)

phi2_arr     = zeros(ComplexF64, N_t)
phi4_arr     = zeros(ComplexF64, N_t)
absphi_arr   = zeros(ComplexF64, N_t)
Smean_arr    = zeros(ComplexF64, N_t)
avgsign_arr  = zeros(ComplexF64, N_t)

println("Computing observables at $N_t t-values (reg=$REG_METHOD, ε=$EPS_REG)...")
for (i, t) in enumerate(t_grid)
    phi2_arr[i], phi4_arr[i], absphi_arr[i], Smean_arr[i], avgsign_arr[i] = compute_observables(t)
    i % 10 == 0 && println("  t=$(round(t,digits=3))  ⟨φ²⟩=$(round(phi2_arr[i],digits=4))  |avg_sign|=$(round(abs(avgsign_arr[i]),sigdigits=3))")
end
println("Done.")

# ── Sanity checks at t=0 ─────────────────────────────────────────────────────
println("\n=== Sanity checks at t=0 ===")
println("  avg_sign(0) should be 1: $(round(avgsign_arr[1], digits=6))")
L0 = domain_limit(0.0)
phi2_real_direct, _ = quadgk(φ -> φ^2 * exp(-S_hat_reg(φ, 0.0)), -L0, L0, rtol=1e-10)
norm_direct,      _ = quadgk(φ ->       exp(-S_hat_reg(φ, 0.0)), -L0, L0, rtol=1e-10)
phi2_euclidean = phi2_real_direct / norm_direct
println("  ⟨φ²⟩ from reweighting: $(round(phi2_arr[1], digits=6))")
println("  ⟨φ²⟩ direct Euclidean: $(round(phi2_euclidean, digits=6))")
println("  Difference: $(round(abs(phi2_arr[1] - phi2_euclidean), sigdigits=3))")

# Free-theory analytic comparison at t=1: ⟨φ²⟩_free = exp(+iπt/2)/m²
println("\n=== Free-theory rotation comparison ===")
for tv in [0.0, 0.5, 1.0]
    idx = argmin(abs.(t_grid .- tv))
    analytic = exp(im*π/2*tv) / M2
    println("  t=$(round(t_grid[idx],digits=2))  ⟨φ²⟩_num=$(round(phi2_arr[idx],digits=4))  free⟨φ²⟩=$(round(analytic,digits=4))")
end

# ── Save ──────────────────────────────────────────────────────────────────────
jldsave("data/observables.jld2";
    t_grid    = collect(t_grid),
    phi2      = phi2_arr,
    phi4      = phi4_arr,
    absphi    = absphi_arr,
    Smean     = Smean_arr,
    avg_sign  = avgsign_arr,
    reg_method = String(REG_METHOD),
    eps_reg    = EPS_REG)
println("\nSaved to data/observables.jld2")

# ── Plots ─────────────────────────────────────────────────────────────────────
t_vec = collect(t_grid)
phi2_free = @. exp(im*π/2*t_vec) / M2   # free-theory analytic ⟨φ²⟩

p1 = plot(t_vec, real.(phi2_arr),
    label="Re⟨φ²⟩", xlabel="t", ylabel="Observable",
    title="Observables vs t  (reg=$REG_METHOD, ε=$EPS_REG)",
    lw=2, color=:blue)
plot!(p1, t_vec, real.(absphi_arr), label="Re⟨|φ|⟩", lw=2, color=:red)
plot!(p1, t_vec, imag.(phi2_arr),   label="Im⟨φ²⟩",  lw=2, color=:blue, linestyle=:dash)
plot!(p1, t_vec, imag.(absphi_arr), label="Im⟨|φ|⟩", lw=2, color=:red, linestyle=:dash)
# free-theory overlay
plot!(p1, t_vec, real.(phi2_free), label="Re⟨φ²⟩ free", lw=1.5, color=:black, linestyle=:dot)
plot!(p1, t_vec, imag.(phi2_free), label="Im⟨φ²⟩ free", lw=1.5, color=:gray,  linestyle=:dot)
savefig(p1, "data/04_observables.png")

p2 = plot(t_vec, abs.(avgsign_arr),
    xlabel="t", ylabel="|avg sign|",
    title="Sign problem severity: |⟨e^{-i Im S}⟩_Q|",
    lw=2, color=:black, ylims=(0, 1.05), legend=false)
hline!(p2, [1.0], linestyle=:dash, color=:gray)
savefig(p2, "data/04_sign_severity.png")

println("Plots saved to data/04_observables.png and data/04_sign_severity.png")
