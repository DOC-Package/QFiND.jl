function spectral_decomposition(sd::SpectralDensity, Temp::Real; npade::Int=5)
    β = ħ * 1e15 / (kb * Temp)
    exponents = ComplexF64[]
    coeffs    = ComplexF64[]

    # Contributions from the spectral density
    poles_J  = sd_poles(sd; scale=icm2ifs)
    res_J    = sd_residues(sd; scale=icm2ifs)
    @assert length(poles_J) == length(res_J) "Number of poles and residues must match"
    for (ω_p, r_p) in zip(poles_J, res_J)
        if imag(ω_p) < 0
            A = r_p * (coth((β*ω_p)/2) + 1)
            push!(exponents, ω_p*1.0im)
            push!(coeffs,    A*(-1.0im))
        end
    end

    # Contributions from the Pade approximation
    if npade > 1
        ξ, η = padeN_Nm1(npade)
        poles_pade = ξ ./ β
        for (ω_p, η_j) in zip(poles_pade, η)
            A = (2η_j / β) * sd(ω_p)
            push!(exponents, ω_p)
            push!(coeffs,    A)
        end
    end

    return exponents, coeffs
end
