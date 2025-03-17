function esprit_decomp(bcf::BosonicBathCorrelationFunction, N_t::Integer, tc::Real, eps::Real) 
    ef = esprit(bcf, 0.0, tc, N_t, eps)
    return (expon=ef.expon, coeff=ef.coeff)
end 