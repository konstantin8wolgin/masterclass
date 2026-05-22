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
# (filled in by Task 8-9)

# %% Phase 4: Distribution comparison plot
# (filled in by Task 10-11)
