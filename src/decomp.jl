function bath_correlation_exponentials_lower(sd::SpectralDensity, beta::Real; npade::Int=5)
    exponents = ComplexF64[]
    coeffs    = ComplexF64[]

    # Contributions from the spectral density
    poles_J  = nodes(sd)
    res_J    = residues(sd)
    @assert length(poles_J) == length(res_J) "Number of poles and residues must match"
    for (ω_p, r_p) in zip(poles_J, res_J)
        if imag(ω_p) < 0
            A = r_p * (coth((beta*ω_p)/2) + 1)
            push!(exponents, ω_p)
            push!(coeffs,    A)
        end
    end

    # Contributions from the Pade approximation 
    pade_poles, η = pade_coth_terms(beta, npade)
    for (ω_p, η_j) in zip(pade_poles, η)
        A = (2η_j / beta) * sd(ω_p)
        push!(exponents, ω_p)
        push!(coeffs,    A)
    end

    return exponents, coeffs
end
