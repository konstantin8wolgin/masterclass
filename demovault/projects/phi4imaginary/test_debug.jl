using QuadGK

const M2 = 1.0
const LAM = 1.0

S_hat_free(φ, t) = M2 * exp(-im*π/2*t)/2 * φ^2

function phi_dot_free_debug(φ_real::Real, t; L=8.0)
    # Compute Sbar_free
    Sbar_free = begin
        num, _ = quadgk(φ -> S_hat_free(φ,t)*exp(-S_hat_free(φ,t)), -L, L, rtol=1e-10)
        den, _ = quadgk(φ -> exp(-S_hat_free(φ,t)), -L, L, rtol=1e-10)
        num/den
    end
    
    println("t=$t, φ_real=$φ_real: Sbar_free = $Sbar_free")
    println("  Re[Sbar_free] = $(real(Sbar_free)), Im[Sbar_free] = $(imag(Sbar_free))")
    
    integrand(x) = (S_hat_free(x, t) - Sbar_free) * exp(-S_hat_free(x, t))
    
    # Check integrand at boundary
    println("  S_hat_free($φ_real, $t) = $(S_hat_free(φ_real, t))")
    println("  integrand at L=$L: $(integrand(L))")
    
    integral, err = quadgk(integrand, φ_real, L, rtol=1e-10)
    println("  integral = $integral, error estimate = $err")
    
    result = (im*π/2) * exp(S_hat_free(φ_real, t)) * integral
    println("  result = $result")
    return result
end

# Test one case
phi_dot_free_debug(1.0, 0.75)
