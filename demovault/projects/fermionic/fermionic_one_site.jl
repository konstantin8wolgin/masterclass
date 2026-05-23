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
    pik_eval(O, U, beta, mu, phi0) -> NamedTuple

Evaluate the holomorphic-PIK estimator of `<O>` on a *pre-drawn* vector
of base samples `phi0 ∼ 𝒩(0, σ² = U/β)`. Returns
`(; est, se_re, se_im, contrib)`:

- `est::ComplexF64` — Monte Carlo mean of the per-sample contributions.
- `se_re`, `se_im::Float64` — standard error of `Re`/`Im`.
- `contrib::Vector{ComplexF64}` — per-sample contributions (mean = est).

Use this directly when you need to evaluate several operators on the
*same* base draw (e.g. so saved `phi0` reproduces every saved
`contrib` byte-for-byte).
"""
function pik_eval(O, U::Float64, beta::Float64, mu::Float64,
                  phi0::AbstractVector{<:Real})
    w, m, _ = mixture_params(U, beta, mu)
    M = length(phi0)
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
    return (; est, se_re, se_im, contrib)
end

"""
    pik_sample(O, U, beta, mu; M, rng) -> NamedTuple

Convenience wrapper around `pik_eval` that draws `M` fresh base samples
`phi0 ∼ 𝒩(0, σ² = U/β)` from `rng` and returns
`(; est, se_re, se_im, phi0, contrib)`.

Each invocation advances `rng` by `M` standard normals — back-to-back
calls produce *different* `phi0` vectors. If you need shared samples
across operators, draw `phi0` once and call `pik_eval` directly.
"""
function pik_sample(O, U::Float64, beta::Float64, mu::Float64;
                    M::Int, rng::AbstractRNG)
    σ = sqrt(U / beta)
    phi0 = randn(rng, M) .* σ
    r = pik_eval(O, U, beta, mu, phi0)
    return (; r.est, r.se_re, r.se_im, phi0, r.contrib)
end

# PIK observables for the Hubbard one-site model.
#
# Critical subtlety: the per-φ trace-ratio operators 2·σ(β(iφ+μ+U/2)),
# 2s+2s², s² (where s = σ(β(iφ+μ+U/2))) have logarithmic-derivative poles
# at Im(φ) = μ + U/2. The strict-interior condition for the contour-shift
# strip is 0 < μ + U/2 < 2U, i.e. -U/2 < μ < 3U/2. On the (U,μ) grid
# {1,4} × {0,1} this holds for ALL 8 parameter cases (μ + U/2 ∈
# {0.5, 1.5, 2.0, 3.0}; strip upper edges 2U ∈ {2, 8}). Contour-deforming
# the integral to R + ij·U therefore crosses these poles in every case
# and gives the WRONG answer.
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
#
# IMPORTANT: all three operators are evaluated on the *same* base draw
# `phi0`. Earlier versions of this notebook used three sequential calls
# to `pik_sample`, which advanced the RNG between operators and saved
# only the first call's `phi0` — making the stored `contrib_n2` and
# `contrib_nud` unreproducible from the saved `phi0`. We now draw once
# and reuse via `pik_eval`.
let
    rng = MersenneTwister(20260522)
    U, beta, mu = U_HL, BETA_HL, MU_HL
    σ = sqrt(U / beta)
    phi0 = randn(rng, M_TABLE) .* σ
    rn   = pik_eval(phi -> n_op(phi,       U, beta, mu), U, beta, mu, phi0)
    rn2  = pik_eval(phi -> n2_op(phi,      U, beta, mu), U, beta, mu, phi0)
    rnup = pik_eval(phi -> nupdown_op(phi, U, beta, mu), U, beta, mu, phi0)
    # Three contour images (Fukushima eq. 17): the deterministic shifts
    # φ₀ → φ₀ + i·j·U for j ∈ {0, 1, 2}. Stored explicitly so consumers
    # can reconstruct the weighted-sample picture without recomputing.
    phi_image_j0 = ComplexF64.(phi0)
    phi_image_j1 = phi0 .+ im * U
    phi_image_j2 = phi0 .+ im * 2U
    jldsave(joinpath(@__DIR__, "data", "samples_headline.jld2");
            phi0         = phi0,
            phi_image_j0 = phi_image_j0,
            phi_image_j1 = phi_image_j1,
            phi_image_j2 = phi_image_j2,
            contrib_n    = rn.contrib,
            contrib_n2   = rn2.contrib,
            contrib_nud  = rnup.contrib,
            U = U, beta = beta, mu = mu, M = M_TABLE)
    f = jldopen(joinpath(@__DIR__, "data", "observables.jld2"), "r")
    analytic = f["analytic"]; close(f)
    jldsave(joinpath(@__DIR__, "data", "observables.jld2");
            analytic = analytic, sampled = sampled_rows, convergence = conv)
    println("Saved data/samples_headline.jld2 (phi0 + 3 shifted images + 3 contribs from a single base draw) and extended data/observables.jld2")
end

# %% Phase 4: Distribution comparison and per-case observable agreement
#
# Two figures:
#  (a) `04_t1_comparison.png` — marquee for the headline case, TWO panels:
#       (left)  histogram of base samples φ₀ vs the analytic
#               Normal(0, σ²=U/β). Each mixture component is a complex
#               Gaussian *shifted* by i·j·U; the contour pull-back to
#               φ₀ ∈ ℝ collapses all three components onto the same real
#               Gaussian, which is what we histogram. This panel is the
#               KS-gated quantity.
#       (right) complex-plane "weighted-samples vs analytic-mixture"
#               picture: the three deterministic contour images
#               φ₀ + i·j·U (j ∈ {0,1,2}) scattered as horizontal bands
#               annotated by their mixture weights w_j = (1, 2e^{βμ},
#               e^{β(2μ−U)}), with the pole line Im(φ) = μ + U/2 and
#               the contour-deformation strip [0, 2U] overlaid.
#  (b) `04_observables_per_case.png` — 2×4 grid, one panel per case,
#       analytic vs sampled bars for ⟨n⟩, ⟨n²⟩, ⟨n↑n↓⟩ (5·SE error bars).
#
# Verification gates:
#   - KS: formal two-sided 1-sample Kolmogorov–Smirnov distance of φ₀
#     vs Normal(0, σ²), i.e. D_M = sup_x |F_M(x) − F(x)| computed as
#     max_i max(|i/M − F(xᵢ)|, |F(xᵢ) − (i−1)/M|). Must be below
#     5×(α=0.05 critical = 1.36/√M).
#   - Reproducibility: the persisted `phi0` reconstructs `contrib_n2`
#     byte-for-byte via `pik_eval` (gate against Finding 1).
#   - Per-case observable agreement was already gated in Phase 3.

# Load the headline base samples and the three shifted contour images.
let
    f = jldopen(joinpath(@__DIR__, "data", "samples_headline.jld2"), "r")
    global PHI0_HL    = f["phi0"]
    global PHI_IMG_J0 = f["phi_image_j0"]
    global PHI_IMG_J1 = f["phi_image_j1"]
    global PHI_IMG_J2 = f["phi_image_j2"]
    global CONTRIB_N2 = f["contrib_n2"]
    @assert f["U"] == U_HL && f["beta"] == BETA_HL && f["mu"] == MU_HL "headline parameters mismatch in samples_headline.jld2"
    close(f)
    @assert PHI_IMG_J0 ≈ ComplexF64.(PHI0_HL)
    @assert PHI_IMG_J1 ≈ PHI0_HL .+ im * U_HL
    @assert PHI_IMG_J2 ≈ PHI0_HL .+ im * 2 * U_HL
    println("Loaded $(length(PHI0_HL)) headline φ₀ samples + 3 shifted images")
end

# Reproducibility gate: re-evaluating ⟨n²⟩ on the saved `phi0` must
# reproduce the saved `contrib_n2` exactly. This catches Finding 1
# (mixed sample streams) if it ever regresses.
let
    rep = pik_eval(phi -> n2_op(phi, U_HL, BETA_HL, MU_HL),
                   U_HL, BETA_HL, MU_HL, PHI0_HL)
    maxdiff = maximum(abs.(rep.contrib .- CONTRIB_N2))
    @assert maxdiff < 1e-12 "saved contrib_n2 not reproducible from saved phi0: max |Δ| = $maxdiff"
    println("Reproducibility: pik_eval(n²) on saved φ₀ reproduces contrib_n2 to $(round(maxdiff; sigdigits=2))  ✓")
end

# Formal two-sided 1-sample KS gate: φ₀ samples vs Normal(0, σ²).
let
    σ = sqrt(U_HL / BETA_HL)
    M = length(PHI0_HL)
    sorted_φ = sort(PHI0_HL)
    F_ana    = cdf.(Normal(0.0, σ), sorted_φ)
    D_plus   = maximum((1:M) ./ M .- F_ana)             # max over upper jumps
    D_minus  = maximum(F_ana .- (0:M-1) ./ M)           # max over lower jumps
    D        = max(D_plus, D_minus)
    crit     = 1.36 / sqrt(M)
    @assert D < 5 * crit "KS distance φ₀ vs Normal(0, σ²) = $D exceeds 5·crit = $(5*crit)"
    println("Formal KS distance φ₀ vs Normal(0, σ²) = $(round(D; sigdigits=3))   (D₊=$(round(D_plus;sigdigits=3)), D₋=$(round(D_minus;sigdigits=3));  5·crit_{α=0.05} = $(round(5*crit; sigdigits=3)))  ✓")
end

# (a) Two-panel marquee: histogram (left) + complex-plane weighted-sample
#     picture (right).
let
    U, beta, mu = U_HL, BETA_HL, MU_HL
    σ = sqrt(U / beta)
    w, m_centers, _ = mixture_params(U, beta, mu)
    Wsum = sum(w)

    # ─ Left panel: KS-gated histogram.
    p_left = plot(title="Base samples φ₀ vs analytic Normal(0, σ²=U/β)",
                  xlabel="φ₀", ylabel="density",
                  legend=:topright, titlefontsize=10)
    histogram!(p_left, PHI0_HL;
               bins=70, normalize=:pdf,
               label="PIK samples (M=$(length(PHI0_HL)))",
               fillalpha=0.55, color=:steelblue, linewidth=0)
    gridφ = range(-4σ, 4σ; length=400)
    plot!(p_left, gridφ, pdf.(Normal(0, σ), gridφ);
          label="Normal(0, σ²=$(round(σ^2; sigdigits=3)))",
          color=:orange, lw=2.5)

    # ─ Right panel: complex-plane scatter of the three weighted contour
    #   images, with the pole line Im(φ)=μ+U/2, the deformation strip
    #   [0, 2U], and the analytic mixture density |e^{−S(φ)}|/Z (real
    #   axis) sketched as a peak-normalised silhouette at Im(φ) = 0.
    nshow = 2000
    idx   = 1:max(1, length(PHI0_HL) ÷ nshow):length(PHI0_HL)
    p_right = plot(title="Weighted contour images φ₀ + i·j·U  vs  mixture density",
                   xlabel="Re φ", ylabel="Im φ",
                   xlims=(-1.2, 1.2),
                   ylims=(-0.5, max(2U, mu + U/2) + 0.5),
                   legend=:topright, titlefontsize=10)
    # Strip of the contour deformation [0, 2U]
    plot!(p_right, [-1.2, 1.2, 1.2, -1.2, -1.2], [0, 0, 2U, 2U, 0];
          seriestype=:shape, color=:gray, fillalpha=0.07, linewidth=0,
          label="deformation strip [0, 2U]")
    # Pole line Im(φ) = μ + U/2 and the n=0,±1 poles
    pole_im = mu + U/2
    plot!(p_right, [-1.2, 1.2], [pole_im, pole_im];
          ls=:dash, color=:red, lw=1.5,
          label="pole line Im φ = μ + U/2")
    for n in -2:2
        pole_re = (2n + 1) * π / beta
        if -1.2 ≤ pole_re ≤ 1.2
            scatter!(p_right, [pole_re], [pole_im];
                     marker=:xcross, color=:red, markersize=7, label="")
        end
    end
    # Three contour bands at Im φ = j·U, weighted by w_j (label shows
    # raw weight; visual marker uses log-scaled alpha so the suppressed
    # j=0 band is still visible).
    rel_logw = log10.(w ./ minimum(w))
    αs = 0.15 .+ 0.7 .* rel_logw ./ maximum(rel_logw)
    for (j, (mj, wj, αj)) in enumerate(zip(m_centers, w, αs))
        scatter!(p_right, real.(PHI0_HL[idx]), fill(mj, length(idx));
                 markersize=2, markerstrokewidth=0,
                 color=[:steelblue, :seagreen, :purple][j], alpha=αj,
                 label="j=$(j-1), Im φ=$(mj),  w=$(round(wj; sigdigits=3))")
    end
    # Silhouette of |e^{-S(φ)}|/max on the real line for shape context.
    gridφ_wide = range(-1.2, 1.2; length=400)
    boltz_abs  = abs.(boltzmann.(gridφ_wide, U, beta, mu))
    silhouette = 0.35 .* boltz_abs ./ maximum(boltz_abs)
    plot!(p_right, gridφ_wide, silhouette .- 0.5;
          color=:black, lw=1.5, label="|e^{−S(φ)}| (real axis, peak-normalised)")

    pmarquee = plot(p_left, p_right; layout=(1, 2), size=(1500, 500),
                    plot_title="Headline (U,β,μ)=($U,$beta,$mu): weighted samples vs analytic mixture density")
    savefig(pmarquee, joinpath(@__DIR__, "data", "04_t1_comparison.png"))
    println("Saved data/04_t1_comparison.png")
end

# (b) 2×4 per-case observable grid.
let
    panels = Any[]
    for c in CASES_VEC
        # Find the matching sampled row.
        sr = nothing
        for row in sampled_rows
            if row.U == c.U && row.β == c.β && row.μ == c.μ
                sr = row; break
            end
        end
        @assert sr !== nothing "missing sampled row for case $c"
        labels = ["⟨n⟩", "⟨n²⟩", "⟨n↑n↓⟩"]
        x      = [1.0, 2.0, 3.0]
        ana    = [sr.n_analytic, sr.n2_analytic, sr.nud_analytic]
        smp    = [sr.n_sampled,  sr.n2_sampled,  sr.nud_sampled]
        ses    = [sr.n_se,       sr.n2_se,       sr.nud_se]
        p = bar(x .- 0.18, ana; bar_width=0.32,
                color=:lightsteelblue, linewidth=0, label="analytic",
                xticks=(x, labels),
                title="U=$(c.U), β=$(c.β), μ=$(c.μ)",
                ylabel="value",
                ylims=(0, max(maximum(ana), maximum(smp)) * 1.25 + 1e-3),
                legend=:topleft, titlefontsize=10, tickfontsize=8,
                legendfontsize=7)
        scatter!(p, x .+ 0.18, smp;
                 yerror=5 .* ses, color=:darkorange,
                 markersize=5, label="sampled (±5·SE)", linewidth=2)
        push!(panels, p)
    end
    pgrid = plot(panels...; layout=(2,4), size=(1700, 750),
                 plot_title="PIK observables vs analytic, all 8 cases (M=$M_TABLE)")
    savefig(pgrid, joinpath(@__DIR__, "data", "04_observables_per_case.png"))
    println("Saved data/04_observables_per_case.png")
end
