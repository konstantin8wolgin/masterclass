module Regulators

export S_reg, weight_reg, L_smooth, S_hat_bare

const M2_DEFAULT  = 1.0
const LAM_DEFAULT = 1.0

# Bare (unregularised) action — same formula as 01_action_check.jl
S_hat_bare(φ, t; m2=M2_DEFAULT, lam=LAM_DEFAULT) =
    m2 * exp(-im*π/2*t)/2 * φ^2 + lam * exp(-im*π/2*t)/4 * φ^4

"""
    S_reg(φ, t; method, ε, m2, lam)

Regularised action. `method` ∈ {:tikhonov, :epsmass, :thimble, :none}.

- :tikhonov  →  S_hat(φ,t) + ε·φ²/2          (real positive Tikhonov damping)
- :epsmass   →  S_hat(φ, t·(1 - iε))         (Feynman iε prescription via t)
- :thimble   →  same as :tikhonov for the *evaluation* (contour handled in caller)
- :none      →  S_hat(φ, t)
"""
function S_reg(φ, t; method::Symbol=:tikhonov, ε::Float64=1e-3,
               m2=M2_DEFAULT, lam=LAM_DEFAULT)
    if method === :tikhonov
        return S_hat_bare(φ, t; m2=m2, lam=lam) + ε * φ^2 / 2
    elseif method === :epsmass
        # Feynman iε via mass shift m² → m² + iε.
        # At t=1: extra contribution = iε·(-i)·φ²/2 = εφ²/2 → real damping.
        return S_hat_bare(φ, t; m2=m2 + im*ε, lam=lam)
    elseif method === :thimble
        return S_hat_bare(φ, t; m2=m2, lam=lam) + ε * φ^2 / 2
    elseif method === :none
        return S_hat_bare(φ, t; m2=m2, lam=lam)
    else
        error("unknown regulator method $method")
    end
end

"""
    weight_reg(φ, t; method, ε)  →  exp(-S_reg(φ,t))
"""
weight_reg(φ, t; method::Symbol=:tikhonov, ε::Float64=1e-3,
           m2=M2_DEFAULT, lam=LAM_DEFAULT) =
    exp(-S_reg(φ, t; method=method, ε=ε, m2=m2, lam=lam))

"""
    L_smooth(t; ε=1e-3)

Smooth integration bound that does not switch abruptly at t→1.
"""
L_smooth(t; ε::Float64=1e-3) = clamp(sqrt(48.0 / (cos(π*t/2) + ε)), 8.0, 220.0)

end # module
