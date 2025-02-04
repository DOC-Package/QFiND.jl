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
    N_w = 2000
    T_c = 500.0
    N_t = 200
    eps = 1e-4
    rank = 27

    sdens = PowerLawExpSD(s, alpha, gamc)
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp, Ω_c)

    res = id_discr(sbeta, bcf, N_t, N_w, T_c, Ω_min, Ω_max, eps)
    @test !isnothing(res)
    freq = res.freq
    coef = res.coef
    calc_error(freq, coef, bcf, T_c, N_t)

    res = id_discr(sbeta, bcf, N_t, N_w, T_c, Ω_min, Ω_max, rank)
    @test !isnothing(res)
    freq = res.freq
    coef = res.coef
    calc_error(freq, coef, bcf, T_c, N_t)
    plot_bcf(freq, coef, bcf, T_c, N_t)

end
