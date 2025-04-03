using Test
using QFiND

@testset "svd.jl" begin
    
    s = 1.0
    alpha = 1.0
    gamc = 50.0
    Temp = 300.0
    Ω_c = 1000.0
    Ω_min = -500.0
    Ω_max = 500.0
    N_w = 10000
    T_c = 1000.0
    N_t = 200
    eps = 1e-4
    rank = 0

    sdens = PowerLawExpSD(s, alpha, gamc)
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp, Ω_c)
    dC = BosonicBCF_dt(sdens, Temp, Ω_c)
    dataset = InitialData(DecompSVD(), sbeta, bcf, dC, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)

    res = svd_intermed_decomp(dataset, eps)
    @test !isnothing(res)
    U = res.basis_time
    V = res.basis_freq
    zk = res.coeff
    Dt = res.Dt
    alpha = Dt * zk
    evaluate_error(U, zk, dataset.bcf)
    #plot_bcf(U, zk, dataset.bcf, dataset.time, "bcf_svd.png")

    evaluate_error(U, alpha, dataset.dC)
    #plot_bcf(U, alpha, dataset.dC, dataset.time, "bcf_dt.png")

    #plot_basis_time(U, dataset.time)
    #plot_basis_freq(V, dataset.freq ./ icm2ifs)

end
