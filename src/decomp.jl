function time_domain_decomp(bcf::BosonicBathCorrelationFunction, N_t::Integer, tc::Real)
    t = collect(range(0, tc, length=N_t))
    cc = bcf.(t)
    
end