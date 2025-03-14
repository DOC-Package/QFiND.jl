include("plot.jl")
using Test
using QFiND

@testset "specdens.jl" begin
    
    s = 1.0
    alpha = 1.0
    gamc = 50.0
    Temp = 300.0
    Ω_c = 1000.0
    Ω_min = -500.0
    Ω_max = 500.0
    N_w = 2000
    T_c = 500.0
    N_t = 200
    eps = 1e-4
    rank = 27

    sdens = PowerLawExpSD(s, alpha, gamc)
    sbeta = BosonicQNSD(sdens, Temp)
    bcf = BosonicBCF(sdens, Temp, Ω_c)

    #plot_qnsd(sdens, 0.0, Ω_max, N_w)
    plot_qnsd(sbeta, Ω_min, Ω_max, N_w)

end
