using Test
using QFiND

@testset "bo.jl" begin
    
    Ω = 1400.0
    Γ = 200.0
    λ = 600.0
    Temp = 300.0
    ub = 6000.0
    Ω_min = -3000.0
    Ω_max = 3000.0
    N_w = 2000
    T_c = 500.0
    N_t = 1000
    eps = 1e-2
    rank = 30

    sdens = BrownianSD(Ω, Γ, λ)
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-6) 
    dataset, _ = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)

    res = id_discr(dataset, eps)
    @test !isnothing(res)
    freq = res.freq
    coeff = res.coeff
    evaluate_error(freq, coeff, bcf, T_c, N_t)

    Temp = 0.0
    eps = 1e-12
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-6) 
    dataset, _ = InitialData(DiscrID(), sbeta, bcf, 0.0, Ω_max, T_c; n_freq=N_w, n_time=N_t)

    res = id_discr(dataset, eps)
    @test !isnothing(res)
    freq = res.freq
    weight = res.weight
    coeff = res.coeff
    Er = reorganization_energy(freq, coeff)
    # error is within 10percent 
    @test abs(Er - λ) / λ < 0.1
    println("Reorganization energy: ", Er)

    res = bsdo_discr(sbeta, 1e-15, Ω_max, 1000; n_lanczos=N_w)
    wk = res.freq
    gk = res.coeff
    Er = reorganization_energy(freq, coeff)
    @test abs(Er - λ) / λ < 0.1
    println("Reorganization energy: ", Er)

end
