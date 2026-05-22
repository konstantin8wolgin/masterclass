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
# (filled in by Task 4-5)

# %% Phase 2: Partition function and analytic observables
# (filled in by Task 6-7)

# %% Phase 3: PIK sampler
# (filled in by Task 8-9)

# %% Phase 4: Distribution comparison plot
# (filled in by Task 10-11)
