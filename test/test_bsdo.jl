include("plot.jl")
using Test
using CairoMakie
using QFiND



@testset "bsdo.jl" begin
    
    s = 1.0
    alpha = 1.0
    gamc = 50.0
    Temp = 300.0
    Ω_c = 1000.0
    Ω_min = -400.0
    Ω_max = 400.0
    N_ω = 2000
    M_sp = 60
    Tc = 500.0
    N_t = 200

    sdens = PowerLawExpSD(s, alpha, gamc)
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp, Ω_c)
    #plot_qnsd(sbeta, Ω_min, Ω_max, N_ω)

    res = bsdo_discr(sbeta, Ω_min, Ω_max, N_ω, M_sp)
    freq = res.freq
    coef = res.coef
    @test !isnothing(res)

    calc_error(freq, coef, bcf, Tc, N_t)
    plot_bcf(freq, coef, bcf, Tc, N_t)

end
