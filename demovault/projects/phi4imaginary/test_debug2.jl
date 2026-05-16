using QuadGK

const M2 = 1.0

S_hat_free(φ, t) = M2 * exp(-im*π/2*t)/2 * φ^2

function test_with_bound(φ_real::Real, t, L)
    Sbar_free = begin
        num, _ = quadgk(φ -> S_hat_free(φ,t)*exp(-S_hat_free(φ,t)), -L, L, rtol=1e-10)
        den, _ = quadgk(φ -> exp(-S_hat_free(φ,t)), -L, L, rtol=1e-10)
        num/den
    end
    
    integrand(x) = (S_hat_free(x, t) - Sbar_free) * exp(-S_hat_free(x, t))
    integral, _ = quadgk(integrand, φ_real, L, rtol=1e-10)
    result = (im*π/2) * exp(S_hat_free(φ_real, t)) * integral
    
    analytic = (im*π/4) * φ_real
    err = abs(result - analytic)
    
    println("L=$L: result=$(round(result, digits=5)), analytic=$(round(analytic, digits=5)), err=$(round(err, sigdigits=3))")
    return err
end

println("Testing φ=1.0, t=0.75:")
for L in [8.0, 12.0, 20.0, 50.0]
    test_with_bound(1.0, 0.75, L)
end

println("\nTesting φ=-1.0, t=0.75:")
for L in [8.0, 12.0, 20.0, 50.0]
    test_with_bound(-1.0, 0.75, L)
end
