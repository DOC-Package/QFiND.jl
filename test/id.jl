include("plot.jl")
using Test
using QFiND

@testset "id.jl" begin
    
    s = 1.0
    alpha = 1.0
    gamc = 50.0
    Temp = 300.0
    Ω_c = 1000.0
    Ω_min = -400.0
    Ω_max = 400.0
    N_ω = 2000
    T_max = 500.0
    N_t = 200
    eps = 1e-4
    rank = 27

    sdens = PowerLawExpSD(s, alpha, gamc)
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp)
    dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t)

    res = id_discr(dataset, eps)
    @test !isnothing(res)
    freq = res.freq
    coeff = res.coeff
    evaluate_error(freq, coeff, bcf, T_max, N_t)

    res = id_discr(sbeta, bcf, Ω_min, Ω_max, T_max, N_ω, N_t, eps)
    @test !isnothing(res)
    freq = res.freq
    coeff = res.coeff
    evaluate_error(freq, coeff, bcf, T_max, N_t)

    res = id_discr(sbeta, bcf, Ω_min, Ω_max, T_max, N_ω, N_t, rank)
    @test !isnothing(res)
    freq = res.freq
    coeff = res.coeff
    evaluate_error(freq, coeff, bcf, T_max, N_t)
    plot_bcf(freq, coeff, bcf, T_max, N_t, "bcf_ohmic_id.png")

end
