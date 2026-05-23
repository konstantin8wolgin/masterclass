# # Fermionic one-site Hubbard sampling via PIK
#
# Reference: `DESIGN.md` in this directory. Reproduces the analytic
# three-Gaussian decomposition of e^{-S} for the one-site spin-1/2
# Hubbard model and samples observables via the constant-velocity PIK
# $\dot\phi = iU$.

# %% Phase 0: Setup, parameters, action
using Pkg
Pkg.activate(@__DIR__)
using Random, Distributions, QuadGK, JLD2, Plots, StatsBase
Random.seed!(20260522)

"""
    S(phi, U, beta, mu)

One-site spin-1/2 Hubbard action with HS auxiliary field `phi`:

    S(phi) = beta/(2U) * phi^2 - 2 * log(1 + exp(beta * (i*phi + mu + U/2)))

Returns a `ComplexF64`. `phi` may be real or complex.
"""
function S(phi, U::Float64, beta::Float64, mu::Float64)
    return beta/(2U) * phi^2 - 2 * log(1 + exp(beta * (im*phi + mu + U/2)))
end

"""
    boltzmann(phi, U, beta, mu)

Returns `exp(-S(phi, U, beta, mu))` as `ComplexF64`.
"""
boltzmann(phi, U, beta, mu) = exp(-S(phi, U, beta, mu))

# Headline parameters used for all single-case plots:
const U_HL, BETA_HL, MU_HL = 1.0, 10.0, 1.0

# Sanity check the action: at phi = 0,
# S(0) = -2 * log(1 + exp(beta * (mu + U/2))).
let
    U, beta, mu = U_HL, BETA_HL, MU_HL
    s0 = S(0.0, U, beta, mu)
    expected = -2 * log(1 + exp(beta * (mu + U/2)))
    @assert isapprox(s0, expected; atol=1e-12) "S(0) mismatch: got $s0, expected $expected"
    println("S(0) = $s0  ✓")
end

# Plot Re S, Im S, |e^{-S}| for the headline parameters over phi in [-6, 6].
let
    U, beta, mu = U_HL, BETA_HL, MU_HL
    phi_grid = range(-6.0, 6.0; length=801)
    s_vals = S.(phi_grid, U, beta, mu)
    p1 = plot(phi_grid, real.(s_vals); xlabel="φ", ylabel="Re S(φ)",
              title="Re S", label="", lw=2)
    p2 = plot(phi_grid, imag.(s_vals); xlabel="φ", ylabel="Im S(φ)",
              title="Im S", label="", lw=2)
    p3 = plot(phi_grid, abs.(exp.(-s_vals)); xlabel="φ", ylabel="|e^{-S}|",
              title="|e^{-S}|", label="", lw=2)
    pcombined = plot(p1, p2, p3; layout=(1,3), size=(1200, 350))
    savefig(pcombined, joinpath(@__DIR__, "data", "01_action_check.png"))
    println("Saved data/01_action_check.png")
end

# %% Phase 1: Analytic three-Gaussian decomposition
#
# Identity:
#   e^{-S(φ)} = Σ_{j=0}^{2}  w_j · exp( -(φ - i*m_j)^2 / (2σ²) )
# with σ² = U/β, m_j = j·U, w = (1, 2 e^{βμ}, e^{β(2μ - U)}).

"""
    mixture_params(U, beta, mu) -> (weights, centers, sigma2)

Returns the analytic three-Gaussian decomposition of `exp(-S(φ; U, β, μ))`
for the one-site spin-1/2 Hubbard model.
- `weights[j+1] = w_j` for j = 0, 1, 2
- `centers[j+1] = m_j` (real numbers; the IMAGINARY shifts are i*m_j)
- `sigma2 = U/β`
"""
function mixture_params(U::Float64, beta::Float64, mu::Float64)
    weights = [1.0, 2*exp(beta*mu), exp(beta*(2mu - U))]
    centers = [0.0, U, 2U]
    sigma2  = U / beta
    return weights, centers, sigma2
end

# Sanity check on the headline case.
let
    w, m, σ2 = mixture_params(U_HL, BETA_HL, MU_HL)
    @assert w ≈ [1.0, 2*exp(BETA_HL*MU_HL), exp(BETA_HL*(2*MU_HL - U_HL))]
    @assert m == [0.0, U_HL, 2*U_HL]
    @assert σ2 ≈ U_HL / BETA_HL
    println("mixture_params headline: w=$w, m=$m, σ²=$σ2  ✓")
end

"""
    p_mix(phi, U, beta, mu) -> ComplexF64

Returns the analytic mixture density evaluated at `phi`.
Equal to `exp(-S(phi, U, beta, mu))` up to floating-point round-off.
"""
function p_mix(phi, U::Float64, beta::Float64, mu::Float64)
    w, m, σ2 = mixture_params(U, beta, mu)
    s = zero(ComplexF64)
    for j in 1:3
        s += w[j] * exp(-(phi - im*m[j])^2 / (2σ2))
    end
    return s
end

# Verification gate: RELATIVE pointwise agreement < 1e-12 over φ in [-6, 6].
# Absolute differences are not meaningful here — at headline (U=1, β=10, μ=1)
# the Boltzmann weight reaches ~1e13 near its peak, so absolute Δ is bounded
# below by ~1e-3 even when the algebra is exact. Use relative error instead.
let
    U, beta, mu = U_HL, BETA_HL, MU_HL
    grid = range(-6.0, 6.0; length=601)
    b   = boltzmann.(grid, U, beta, mu)
    pm  = p_mix.(grid, U, beta, mu)
    rel = abs.(pm .- b) ./ max.(abs.(b), 1.0)   # denominator clamp avoids /0 in tails
    err = maximum(rel)
    @assert err < 1e-12 "mixture vs Boltzmann relative mismatch: max|Δ|/|e^{-S}| = $err"
    println("max relative |p_mix - e^{-S}|/|e^{-S}| over φ∈[-6,6] = $err  ✓")
end

# Plot Re and Im of both functions overlaid.
let
    U, beta, mu = U_HL, BETA_HL, MU_HL
    grid = range(-6.0, 6.0; length=601)
    pm = p_mix.(grid, U, beta, mu)
    b  = boltzmann.(grid, U, beta, mu)
    pre = plot(grid, real.(b); label="Re e^{-S}", lw=2)
    plot!(pre, grid, real.(pm); label="Re p_mix", lw=2, ls=:dash)
    pim = plot(grid, imag.(b); label="Im e^{-S}", lw=2)
    plot!(pim, grid, imag.(pm); label="Im p_mix", lw=2, ls=:dash)
    pcomb = plot(pre, pim; layout=(1,2), size=(1100, 400),
                 xlabel="φ", title="Mixture vs e^{-S}, headline (U,β,μ)=(1,10,1)")
    savefig(pcomb, joinpath(@__DIR__, "data", "01_mixture_vs_target.png"))
    println("Saved data/01_mixture_vs_target.png")
end

# %% Phase 2: Partition function and analytic observables
#
# Direct Fock-trace formulas:
#   Z          = 1 + 2 e^{βμ} + e^{β(2μ - U)}
#   ⟨n⟩       = (2 e^{βμ} + 2 e^{β(2μ - U)}) / Z
#   ⟨n²⟩       = (2 e^{βμ} + 4 e^{β(2μ - U)}) / Z
#   ⟨n↑ n↓⟩  = e^{β(2μ - U)} / Z

"""
    fock_trace_observables(U, beta, mu) -> NamedTuple

Closed-form `Z, ⟨n⟩, ⟨n²⟩, ⟨n↑n↓⟩` from summing the 4 Fock states.
"""
function fock_trace_observables(U::Float64, beta::Float64, mu::Float64)
    z1 = 2 * exp(beta * mu)
    z2 = exp(beta * (2mu - U))
    Z = 1.0 + z1 + z2
    n_mean   = (z1 + 2 * z2) / Z
    n2_mean  = (z1 + 4 * z2) / Z
    nup_ndown = z2 / Z
    return (; Z, n_mean, n2_mean, nup_ndown)
end

"""
    gaussian_integral_Z(U, beta, mu) -> Float64

Z computed from the three-Gaussian decomposition by closed-form
integration of each component:
    ∫ e^{-(φ - i m_j)^2 / (2σ²)} dφ = sqrt(2π σ²)
Multiplied by the HS normalisation N = sqrt(β/(2π U)) the result must
match the Fock-trace `Z`.
"""
function gaussian_integral_Z(U::Float64, beta::Float64, mu::Float64)
    w, _, σ2 = mixture_params(U, beta, mu)
    N = sqrt(beta / (2π * U))
    return N * sum(w) * sqrt(2π * σ2)
end

# Verification gate: closed-form Z from the mixture equals the Fock trace.
let
    cases = [(U, β, μ) for U in (1.0, 4.0), β in (2.0, 10.0), μ in (0.0, 1.0)]
    for (U, β, μ) in cases
        z_trace = fock_trace_observables(U, β, μ).Z
        z_gauss = gaussian_integral_Z(U, β, μ)
        @assert isapprox(z_trace, z_gauss; atol=1e-12, rtol=1e-12) "Z mismatch at (U,β,μ)=($U,$β,$μ): trace=$z_trace gauss=$z_gauss"
    end
    println("Z(Fock trace) == Z(Gaussian integral) for all 8 cases  ✓")
end

# Build the 8-case grid and tabulate analytic observables.
const CASES = [(U=U, β=β, μ=μ)
               for U in (1.0, 4.0), β in (2.0, 10.0), μ in (0.0, 1.0)]
const CASES_VEC = vec(CASES)

analytic_rows = map(CASES_VEC) do c
    obs = fock_trace_observables(c.U, c.β, c.μ)
    (; c.U, c.β, c.μ, obs.Z, obs.n_mean, obs.n2_mean, obs.nup_ndown)
end

println("Analytic observables per case:")
println(" U    β     μ     Z            ⟨n⟩      ⟨n²⟩     ⟨n↑n↓⟩")
for r in analytic_rows
    println(rpad(string(r.U), 5), rpad(string(r.β), 6), rpad(string(r.μ), 6),
            rpad(string(round(r.Z;        digits=6)), 13),
            rpad(string(round(r.n_mean;    digits=4)),  9),
            rpad(string(round(r.n2_mean;   digits=4)),  9),
            rpad(string(round(r.nup_ndown; digits=4)),  9))
end

# Persist the analytic table and the mixture parameters per case.
let
    mixture_per_case = map(CASES_VEC) do c
        w, m, σ2 = mixture_params(c.U, c.β, c.μ)
        (; c.U, c.β, c.μ, weights=w, centers=m, sigma2=σ2)
    end
    jldsave(joinpath(@__DIR__, "data", "observables.jld2");
            analytic = analytic_rows)
    jldsave(joinpath(@__DIR__, "data", "mixture_params.jld2");
            mixtures = mixture_per_case)
    println("Saved data/observables.jld2 and data/mixture_params.jld2")
end

# %% Phase 3: PIK sampler
#
# Estimator (Fukushima 2026 eq. 20):
#   ⟨O⟩ = Σ_j w_j · (1/M) Σ_i O(φ₀^{(i)} + i m_j)  /  Σ_j w_j
# with φ₀ ∼ 𝒩(0, σ² = U/β) drawn M times.

"""
    pik_sample(O, U, beta, mu; M, rng) -> NamedTuple

Returns the PIK estimate of <O> from `M` real-Gaussian samples as
`(; est, se_re, se_im, phi0, contrib)`:

- `est::ComplexF64` — Monte Carlo mean of the per-sample contributions.
- `se_re`, `se_im::Float64` — standard error of `real`/`imag` parts.
  The imaginary part of `<O>` is zero *in expectation*; at finite `M`
  it has a fluctuation of the same order as `se_im`. Gate Im with a
  statistical bound, not a hard tolerance.
- `phi0::Vector{Float64}` — base samples (kept for histograms).
- `contrib::Vector{ComplexF64}` — per-sample weighted contributions
  (their mean is `est`; their std/√M is the SE).
"""
function pik_sample(O, U::Float64, beta::Float64, mu::Float64;
                    M::Int, rng::AbstractRNG)
    w, m, σ2 = mixture_params(U, beta, mu)
    σ = sqrt(σ2)
    phi0 = randn(rng, M) .* σ      # 𝒩(0, σ²)
    Wtot = sum(w)
    contrib = zeros(ComplexF64, M)
    for j in 1:3
        @inbounds for i in 1:M
            contrib[i] += w[j] * O(phi0[i] + im*m[j])
        end
    end
    contrib ./= Wtot
    est = sum(contrib) / M
    se_re = std(real.(contrib); corrected=true) / sqrt(M)
    se_im = std(imag.(contrib); corrected=true) / sqrt(M)
    return (; est, se_re, se_im, phi0, contrib)
end

# PIK observables for the Hubbard one-site model.
#
# Critical subtlety: the per-φ trace-ratio operators 2·σ(β(iφ+μ+U/2)),
# 2s+2s², s² (where s = σ(β(iφ+μ+U/2))) have logarithmic-derivative poles
# at Im(φ) = μ + U/2 — inside the contour-shift strip [0, 2U] for 7 of
# the 8 parameter cases. Contour-deforming the integral to R + ij·U
# crosses these poles and gives the WRONG answer.
#
# The fix is the HOLOMORPHIC representation from Fukushima eq. 21:
# the auxiliary field IS the density, n = -i⟨φ⟩/U. Polynomial in φ,
# so every contour shift is exact. Derivation by direct Gaussian
# integration of e^{-S(φ)}·φ^k against the (1+ξ)² expansion:
#
#   ⟨-iφ⟩_phys    = U · ⟨n⟩
#   ⟨φ²⟩_phys    = U/β − U² · ⟨n²⟩
#   n↑n↓         = (n² − n)/2   (fermion algebra)
#
# Yielding the holomorphic PIK observables:
n_op(phi, U, beta, mu)       = -im * phi / U
n2_op(phi, U, beta, mu)      = 1/(beta*U) - phi^2 / U^2
nupdown_op(phi, U, beta, mu) = (n2_op(phi, U, beta, mu) - n_op(phi, U, beta, mu)) / 2

# Sanity check the sampler.
# ⟨n⟩ via O = -iφ/U has identically constant real part per sample
# (Re(-iφ_i/U) = 0 for real φ_i, plus the constant Fock value), so the
# estimator is exact to machine precision and SE_Re ≈ 1e-18. Gate it
# with an absolute tolerance instead of |Δ|/SE.
# ⟨n²⟩ has nontrivial variance from the φ² term; gate statistically.
let
    U, beta, mu = U_HL, BETA_HL, MU_HL
    rng = MersenneTwister(20260522)
    r_n  = pik_sample(phi -> n_op(phi,  U, beta, mu), U, beta, mu; M = 100_000, rng = rng)
    r_n2 = pik_sample(phi -> n2_op(phi, U, beta, mu), U, beta, mu; M = 100_000, rng = rng)
    ana = fock_trace_observables(U, beta, mu)
    @assert abs(real(r_n.est)  - ana.n_mean)  < 1e-10 "⟨n⟩ off by $(abs(real(r_n.est) - ana.n_mean))"
    @assert abs(real(r_n2.est) - ana.n2_mean) < 5*r_n2.se_re "⟨n²⟩ off by $(abs(real(r_n2.est) - ana.n2_mean)) > 5·SE=$(5*r_n2.se_re)"
    println("Headline ⟨n⟩  = $(real(r_n.est))   (analytic $(ana.n_mean), exact to 1e-10)  ✓")
    println("Headline ⟨n²⟩ = $(real(r_n2.est)) ± $(r_n2.se_re)   (analytic $(ana.n2_mean))  ✓")
end

# Full 8-case sampling table at M = 10^5. Gate per observable:
#   ⟨n⟩      — absolute (zero-variance estimator, machine-precision)
#   ⟨n²⟩     — statistical, |Δ| < 5·SE
#   ⟨n↑n↓⟩   — statistical, |Δ| < 5·SE
const M_TABLE = 100_000

sampled_rows = let
    rng = MersenneTwister(20260522)
    rows = []
    for c in CASES_VEC
        U, β, μ = c.U, c.β, c.μ
        analytic = fock_trace_observables(U, β, μ)
        rn   = pik_sample(phi -> n_op(phi,       U, β, μ), U, β, μ; M = M_TABLE, rng = rng)
        rn2  = pik_sample(phi -> n2_op(phi,      U, β, μ), U, β, μ; M = M_TABLE, rng = rng)
        rnup = pik_sample(phi -> nupdown_op(phi, U, β, μ), U, β, μ; M = M_TABLE, rng = rng)
        @assert abs(real(rn.est) - analytic.n_mean) < 1e-10 "Case ($U,$β,$μ): ⟨n⟩ off by $(abs(real(rn.est) - analytic.n_mean))"
        for (name, r, ana) in (("⟨n²⟩", rn2, analytic.n2_mean),
                               ("⟨n↑n↓⟩", rnup, analytic.nup_ndown))
            diff = abs(real(r.est) - ana)
            @assert diff < 5 * r.se_re "Case ($U,$β,$μ): $name diff=$diff > 5·SE=$(5*r.se_re)"
        end
        push!(rows, (; U, β, μ,
                     n_sampled   = real(rn.est),   n_se   = rn.se_re,   n_analytic   = analytic.n_mean,
                     n2_sampled  = real(rn2.est),  n2_se  = rn2.se_re,  n2_analytic  = analytic.n2_mean,
                     nud_sampled = real(rnup.est), nud_se = rnup.se_re, nud_analytic = analytic.nup_ndown))
    end
    rows
end

println("Sampled vs analytic (M = $M_TABLE, all 8 cases):")
println(" U    β     μ    ⟨n⟩ (exact)    ⟨n²⟩ sampled (SE)      analytic    |Δ|/SE")
for r in sampled_rows
    z = abs(r.n2_sampled - r.n2_analytic) / r.n2_se
    println(rpad(string(r.U), 5), rpad(string(r.β), 6), rpad(string(r.μ), 5),
            rpad(string(round(r.n_sampled; digits=5)), 15),
            rpad("$(round(r.n2_sampled; digits=5)) ± $(round(r.n2_se; sigdigits=2))", 25),
            rpad(string(round(r.n2_analytic; digits=5)), 12),
            round(z; digits=2))
end

# Headline convergence study: ⟨n²⟩ at M ∈ {10³, 10⁴, 10⁵, 10⁶}. The SE
# should scale as 1/√M. (⟨n⟩ is zero-variance under polynomial PIK at one
# site, so there is nothing to converge — ⟨n²⟩ is the natural choice.)
const M_CONV = (10^3, 10^4, 10^5, 10^6)

conv = let
    U, beta, mu = U_HL, BETA_HL, MU_HL
    rng = MersenneTwister(20260601)
    Ms = Int[]
    ests = Float64[]
    ses = Float64[]
    for M in M_CONV
        r = pik_sample(phi -> n2_op(phi, U, beta, mu), U, beta, mu; M = M, rng = rng)
        push!(Ms, M); push!(ests, real(r.est)); push!(ses, r.se_re)
        println("M = $(lpad(M, 7)):  ⟨n²⟩ = $(round(real(r.est); digits=6)) ± $(round(r.se_re; sigdigits=3))")
    end
    (; Ms, ests, ses)
end

let
    U, beta, mu = U_HL, BETA_HL, MU_HL
    analytic = fock_trace_observables(U, beta, mu).n2_mean
    Mvec = collect(conv.Ms)
    ref = conv.ses[1] .* sqrt.(Mvec[1] ./ Mvec)
    p = plot(Mvec, conv.ses; xscale=:log10, yscale=:log10,
             marker=:circle, lw=2, label="measured SE",
             xlabel="M (samples)", ylabel="SE[Re ⟨n²⟩]",
             title="Headline (U,β,μ)=($U,$beta,$mu): 1/√M convergence of ⟨n²⟩")
    plot!(p, Mvec, ref; ls=:dash, lw=2, label="1/√M reference")
    savefig(p, joinpath(@__DIR__, "data", "03_convergence.png"))
    println("Saved data/03_convergence.png  (analytic ⟨n²⟩ = $analytic)")
end

# Persist headline base samples + per-observable contributions, and
# extend the analytic table with sampled rows + convergence data.
let
    rng = MersenneTwister(20260522)
    U, beta, mu = U_HL, BETA_HL, MU_HL
    rn   = pik_sample(phi -> n_op(phi,       U, beta, mu), U, beta, mu; M = M_TABLE, rng = rng)
    rn2  = pik_sample(phi -> n2_op(phi,      U, beta, mu), U, beta, mu; M = M_TABLE, rng = rng)
    rnup = pik_sample(phi -> nupdown_op(phi, U, beta, mu), U, beta, mu; M = M_TABLE, rng = rng)
    jldsave(joinpath(@__DIR__, "data", "samples_headline.jld2");
            phi0      = rn.phi0,
            contrib_n = rn.contrib,
            contrib_n2  = rn2.contrib,
            contrib_nud = rnup.contrib,
            U = U, beta = beta, mu = mu, M = M_TABLE)
    f = jldopen(joinpath(@__DIR__, "data", "observables.jld2"), "r")
    analytic = f["analytic"]; close(f)
    jldsave(joinpath(@__DIR__, "data", "observables.jld2");
            analytic = analytic, sampled = sampled_rows, convergence = conv)
    println("Saved data/samples_headline.jld2 and extended data/observables.jld2")
end

# %% Phase 4: Distribution comparison plot
# (filled in by Task 10-11)
