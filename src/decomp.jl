function spectral_decomposition(poles_J, residues_J, poles_f, residues_f, f_func, J_func)

    in_contour(ω) = (imag(ω) < 0) 
    #in_contour(w) = (t >= 0 && imag(w) < 0) || (t < 0 && imag(w) > 0)

    # Dictionary for composing contributions from each pole
    terms = Dict{ComplexF64, ComplexF64}()

    # Contribution from poles of J(ω): resJ * f_func(ω)
    for (ω, resJ) in zip(poles_J, residues_J)
        if in_contour(ω)
            terms[ω] = get(terms, ω, 0.0 + 0im) + resJ * f_func(ω)
        end
    end

    # Contribution from poles of f(ω): J_func(ω) * resf
    for (ω, resf) in zip(poles_f, residues_f)
        if in_contour(ω)
            terms[ω] = get(terms, ω, 0.0 + 0im) + J_func(ω) * resf
        end
    end

    coeff = ComplexF64[]
    expon = ComplexF64[]
    for (ω, r) in terms
         push!(expon, ω)
         push!(coeff, 2π * im * r)
    end

    return (coeff=coeff, expon=expon)
end

function esprit_decomp(bcf::BosonicBathCorrelationFunction, N_t::Integer, tc::Real, eps::Real) 
    ef = esprit(bcf, 0.0, tc, N_t, eps)
    return (expon=ef.expon, coeff=ef.coeff)
end 
