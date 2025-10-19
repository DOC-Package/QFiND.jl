"""
    ChebyshevExpansion

Structure to store Chebyshev expansion coefficients for bath correlation functions.
"""
struct ChebyshevExpansion
    ω_min::Float64
    ω_max::Float64
    ω_bar::Float64  
    Ω::Float64      
    coeffs::Vector{ComplexF64} 
    basis::Vector{Function}
    n_terms::Int
end

"""
    chebyshev_expansion(sd::SpectralDensity, Temp::Real, ω_min::Real, ω_max::Real, n_terms::Int; 
                       scale::Float64=1.0, rtol::Real=1e-8, atol::Real=1e-12)

Compute Chebyshev expansion coefficients for a bath correlation function using QuadGK.
"""
function chebyshev_expansion(sd::SpectralDensity, Temp::Real, ω_min::Real, ω_max::Real, 
                           n_terms::Int; scale::Float64=icm2ifs, rtol::Real=1e-8)
    ω_bar = (ω_max + ω_min) / 2 * scale
    Ω = (ω_max - ω_min) / 2 * scale
    β = ħ * 1e15 / (kb * Temp) 
    coeffs = zeros(ComplexF64, n_terms)
    
    # Construct basis functions
    basis = Vector{Function}(undef, n_terms)
    for k in 0:(n_terms-1)
        basis[k+1] = t -> besselj(k, Ω * t) * exp(-1im * ω_bar * t)
    end
    
    for k in 0:(n_terms-1)
        integrand = x -> begin
            ω = Ω * x + ω_bar
            sbeta = BosonicQNSD(sd, Temp)(ω; scale=scale)
            # Chebyshev polynomial T_k
            T_k = cos(k * acos(x))
            prefactor = (k == 0) ? 1.0 : 2.0
            weight = 1 / sqrt(1 - x^2)
            #return prefactor * (-1im)^k * T_k * J_val * thermal_factor * weight
            return prefactor * (-1im)^k * T_k * sbeta * weight
        end
        result1, _ = quadgk(integrand, -1, 0; rtol=rtol)
        result2, _ = quadgk(integrand, 0, 1; rtol=rtol)
        result = result1 + result2
        coeffs[k+1] = (Ω / π) * result 
    end
    return ChebyshevExpansion(ω_min, ω_max, ω_bar, Ω, coeffs, basis, n_terms)
end

"""
    evaluate_correlation(cheb::ChebyshevExpansion, t::Real)

Evaluate the correlation function at time t using Chebyshev expansion.
"""
function evaluate_correlation(cheb::ChebyshevExpansion, t::Real)
    result = 0.0 + 0.0im
    for k in 1:cheb.n_terms
        result += cheb.coeffs[k] * cheb.basis[k](t)
    end
    return result
end

(cheb::ChebyshevExpansion)(t::Real) = evaluate_correlation(cheb, t)

"""
    chebyshev_bcf(sd::SpectralDensity, Temp::Real, ω_min::Real, ω_max::Real, 
                  n_terms::Int; kwargs...)

Create a bath correlation function using Chebyshev expansion.
"""
function chebyshev_bcf(sd::SpectralDensity, Temp::Real, ω_min::Real, ω_max::Real, 
                      n_terms::Int; kwargs...)
    cheb = chebyshev_expansion(sd, Temp, ω_min, ω_max, n_terms; kwargs...)
    return t -> evaluate_correlation(cheb, t)
end

"""
    chebyshev_bcf(cheb::ChebyshevExpansion)

Create a bath correlation function using Chebyshev expansion.
"""
function chebyshev_bcf(cheb::ChebyshevExpansion)
    return t -> evaluate_correlation(cheb, t)
end