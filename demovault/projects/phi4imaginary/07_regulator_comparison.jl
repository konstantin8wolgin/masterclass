using QuadGK, JLD2, Plots, Printf

include("reg.jl")
using .Regulators

data = load("data/action_check.jld2")
const M2  = data["M2"]
const LAM = data["LAM"]

S_bare(φ, t) = Regulators.S_hat_bare(φ, t; m2=M2, lam=LAM)

function observables_along_contour(t; method::Symbol, ε::Float64)
    if method === :thimble
        # ── Contour-deformed integral: φ_c(φ₀, t) = φ₀ · exp(i π t/4)
        # This rotation makes the quadratic term on the contour purely real-positive
        # (M2·φ₀²/2 — Gaussian σ=1 damping for all t). The quartic becomes
        # LAM·φ₀⁴·exp(iπt/2)/4 — Re ≥ 0 for t ≤ 1. Tractable everywhere on [0,1].
        # Free-theory limit: ⟨φ²⟩ = exp(iπt/2) · 1 → i at t=1 (matches analytic).
        # On-contour quadratic damping is M²·φ₀²/2 (σ=1) for all t → fixed L=12 suffices.
        L = 12.0
        J = exp(im*π/4*t)
        φc(φ0) = φ0 * J
        Sval(φ0) = S_bare(φc(φ0), t)
        w(φ0)    = exp(-Sval(φ0)) * J
        # Quenched (|w|) for sign normalisation
        wq(φ0)   = abs(w(φ0))

        denom, _    = quadgk(w,                          -L, L, rtol=1e-8)
        num_phi2, _ = quadgk(φ0 -> φc(φ0)^2 * w(φ0),     -L, L, rtol=1e-8)
        num_abs, _  = quadgk(φ0 -> abs(φc(φ0)) * w(φ0),  -L, L, rtol=1e-8)
        norm_Q, _   = quadgk(wq,                         -L, L, rtol=1e-8)
        return num_phi2/denom, num_abs/denom, denom/norm_Q
    else
        L = Regulators.L_smooth(t; ε=ε)
        Sr(φ) = Regulators.S_reg(φ, t; method=method, ε=ε, m2=M2, lam=LAM)
        phase(φ)   = exp(-im * imag(Sr(φ)))
        quenched(φ)= exp(-real(Sr(φ)))
        weight(φ)  = phase(φ) * quenched(φ)
        denom, _    = quadgk(weight,                    -L, L, rtol=1e-8)
        norm_Q, _   = quadgk(quenched,                  -L, L, rtol=1e-8)
        num_phi2, _ = quadgk(φ -> φ^2  * weight(φ),     -L, L, rtol=1e-8)
        num_abs, _  = quadgk(φ -> abs(φ) * weight(φ),   -L, L, rtol=1e-8)
        return num_phi2/denom, num_abs/denom, denom/norm_Q
    end
end

t_grid = range(0.0, 1.0, length=50)

regs = [
    (:tikhonov, 1e-3, :blue,   "Tikhonov ε=1e-3"),
    (:epsmass,  1e-3, :green,  "iε-prescription δ=1e-3"),
    (:thimble,  1e-3, :red,    "free-theory thimble contour"),
]

results = Dict{Symbol, NamedTuple}()
for (m, ε, _, _) in regs
    println("Running regulator $m (ε=$ε)...")
    phi2   = zeros(ComplexF64, length(t_grid))
    absphi = zeros(ComplexF64, length(t_grid))
    sgn    = zeros(ComplexF64, length(t_grid))
    for (i, t) in enumerate(t_grid)
        phi2[i], absphi[i], sgn[i] = observables_along_contour(t; method=m, ε=ε)
    end
    results[m] = (phi2=phi2, absphi=absphi, sign=sgn)
end

# Baseline (no regulator) for t ≤ 0.95
println("Running baseline (unregulated) for t ≤ 0.95...")
t_base = filter(t -> t ≤ 0.95, collect(t_grid))
phi2_base   = zeros(ComplexF64, length(t_base))
absphi_base = zeros(ComplexF64, length(t_base))
sign_base   = zeros(ComplexF64, length(t_base))
for (i, t) in enumerate(t_base)
    phi2_base[i], absphi_base[i], sign_base[i] =
        observables_along_contour(t; method=:none, ε=0.0)
end

# Plot 2×2
p_re  = plot(title="Re⟨φ²⟩",   xlabel="t", ylabel="Re⟨φ²⟩",   legend=:topright)
p_im  = plot(title="Im⟨φ²⟩",   xlabel="t", ylabel="Im⟨φ²⟩",   legend=:topright)
p_abs = plot(title="Re⟨|φ|⟩",  xlabel="t", ylabel="Re⟨|φ|⟩",  legend=:topright)
p_sgn = plot(title="|avg sign|", xlabel="t", ylabel="|avg sign|",
             ylims=(0,1.05), legend=:topright)

for (m, ε, c, label) in regs
    r = results[m]
    plot!(p_re,  t_grid, real.(r.phi2),   lw=2, color=c, label=label)
    plot!(p_im,  t_grid, imag.(r.phi2),   lw=2, color=c, label=label)
    plot!(p_abs, t_grid, real.(r.absphi), lw=2, color=c, label=label)
    plot!(p_sgn, t_grid, abs.(r.sign),    lw=2, color=c, label=label)
end
plot!(p_re,  t_base, real.(phi2_base),   lw=1.5, color=:black, linestyle=:dash, label="baseline")
plot!(p_im,  t_base, imag.(phi2_base),   lw=1.5, color=:black, linestyle=:dash, label="baseline")
plot!(p_abs, t_base, real.(absphi_base), lw=1.5, color=:black, linestyle=:dash, label="baseline")
plot!(p_sgn, t_base, abs.(sign_base),    lw=1.5, color=:black, linestyle=:dash, label="baseline")

# free-theory overlay
plot!(p_re, t_grid, [real(exp(im*π/2*t)/M2) for t in t_grid],
      lw=1.5, color=:gray, linestyle=:dot, label="free")
plot!(p_im, t_grid, [imag(exp(im*π/2*t)/M2) for t in t_grid],
      lw=1.5, color=:gray, linestyle=:dot, label="free")

savefig(plot(p_re, p_im, p_abs, p_sgn, layout=(2,2), size=(1300, 900)),
        "data/07_regulator_comparison.png")
println("Saved data/07_regulator_comparison.png")

jldsave("data/regulator_comparison.jld2";
    t_grid = collect(t_grid),
    tikhonov_phi2 = results[:tikhonov].phi2,
    epsmass_phi2  = results[:epsmass].phi2,
    thimble_phi2  = results[:thimble].phi2,
    tikhonov_sign = results[:tikhonov].sign,
    epsmass_sign  = results[:epsmass].sign,
    thimble_sign  = results[:thimble].sign)

# Table at t = 1
println("\n=== Values at t = 1.0 ===")
@printf("%-32s  %-18s  %-18s  %-10s\n", "Regulator", "Re⟨φ²⟩", "Im⟨φ²⟩", "|sign|")
for (m, ε, _, label) in regs
    r = results[m]
    p2 = r.phi2[end]; sg = r.sign[end]
    @printf("%-32s  %+.6f       %+.6f       %.4f\n", label, real(p2), imag(p2), abs(sg))
end
analytic = exp(im*π/2)/M2
@printf("%-32s  %+.6f       %+.6f\n", "free-theory analytic", real(analytic), imag(analytic))
