using Test
using CairoMakie
using QFiND

@testset "bsdo.jl" begin
    
    s = 1.0
    α = 1.0
    γ_c = 50.0
    Temp = 300.0
    Ω_c = 500.0
    Ω_min = -300.0
    Ω_max = 300.0
    N_freq = 4000
    M_sample = 50

    ple = PowerLawExpSD(s, α, γ_c)
    sdens = make_sdens(ple; scale=icm2ifs)
    sbeta = make_sbeta(ple, Temp)

    bcf = make_bcf(sdens, Temp, Ω_c)
    res = bsdo_discr(sbeta, Omega_min, Omega_max, N_w, Msp)
    wk = res.freq
    zk = res.coef
    calc_error(wk, zk, bcf, 500.0, 200)

end
