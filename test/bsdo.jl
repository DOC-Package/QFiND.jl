using Test
using QFiND

@testset "bsdo.jl" begin
    
    s = 1.0
    alpha = 50.0
    gamc = 50.0
    Temp = 300.0
    ub = 1000.0
    Ω_min = -400.0
    Ω_max = 400.0
    N_w = 2000
    M_sp = 100
    T_c = 500.0
    N_t = 200

    sdens = PowerLawExpSD(s, gamc; alpha=alpha)
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-6)
    dataset, _ = InitialData(DiscrBSDO(), sbeta, bcf, Ω_min, Ω_max, T_c; n_lanczos=N_w, n_time=N_t)
    res = bsdo_discr(dataset, M_sp)
    wk = res.freq
    gk = res.coeff
    @test !isnothing(res)
    evaluate_error(wk, gk, bcf, T_c, N_t)
end
