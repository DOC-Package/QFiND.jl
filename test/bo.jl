include("plot.jl")
using Test
using QFiND

@testset "id.jl" begin
    
    Ω = 1400.0
    Γ = 200.0
    λ = 600.0
    Temp = 300.0
    Ω_c = 6000.0
    Ω_min = -3000.0
    Ω_max = 3000.0
    N_w = 2000
    T_c = 500.0
    N_t = 1000
    eps = 1e-2
    rank = 27

    sdens = BrownianSD(Ω, Γ, λ)
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp, Ω_c)
    dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)

    res = id_discr(dataset, eps)
    @test !isnothing(res)
    freq = res.freq
    coef = res.coef
    evaluate_error(freq, coef, bcf, T_c, N_t)
    plot_bcf(freq, coef, bcf, T_c, N_t, "bcf_bo.png")

end
