using Test
using QFiND

@testset "id.jl" begin
    
    s = 1.0
    alpha = 1.0
    gamc = 50.0
    Temp = 300.0
    Ω_c = 1000.0
    Omega_min = -400.0
    Omega_max = 400.0
    N_w = 2000
    Msp = 50
    Tc = 500.0
    N_t = 200
    eps = 1e-4
    rank = 27

    ple = PowerLawExpSD(s, alpha, gamc)
    sdens = make_sdens(ple; scale=icm2ifs)
    sbeta = make_sbeta(ple, Temp; scale=icm2ifs)

    bcf = make_bcf(sdens, Temp, Ω_c)
    res = id_discr(sbeta, bcf, N_t, N_w, Tc, Omega_min, Omega_max, eps)
    wk = res.freq
    gk = res.coef
    calc_error(wk, gk, bcf, 500.0, 200)

    res = id_discr(sbeta, bcf, N_t, N_w, Tc, Omega_min, Omega_max, rank)
    wk = res.freq
    gk = res.coef
    calc_error(wk, gk, bcf, 500.0, 200)

end
