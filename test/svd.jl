include("plot.jl")
using Test
using QFiND

@testset "svd.jl" begin
    
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
    eps = 5e-3
    rank = 27

    sdens = PowerLawExpSD(s, alpha, gamc)
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp, Ω_c)
    dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)

    res = svd_basis_decomp(dataset, eps)
    @test !isnothing(res)
    U = res.basis
    zk = res.coef
    evaluate_error(U, zk, dataset.bcf)
    plot_bcf(U, zk, dataset.bcf, dataset.time)

end
