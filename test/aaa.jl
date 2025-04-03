using Test
using QFiND
using RationalFunctionApproximation

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
    w = collect(range(1e-12, Ω_max, length=N_w))
    Jw = sdens.(w)
    bary = aaa(w, Jw; tol=1e-10)

    sdens = AAAfittedSD(bary)
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp, Ω_c)
    dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)

    res = id_discr(sbeta, bcf, N_t, N_w, T_c, Ω_min, Ω_max, rank)
    @test !isnothing(res)
    freq = res.freq
    coef = res.coef
    evaluate_error(freq, coef, bcf, T_c, N_t)
    plot_bcf(freq, coef, bcf, T_c, N_t)

end
