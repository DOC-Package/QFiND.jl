using Test
using CairoMakie
using QFiND

@testset "bsdo.jl" begin
    
    s = 1.0
    alpha = 1.0
    gamc = 50.0
    Temp = 300.0
    Ω_c = 500.0
    Omega_min = -300.0
    Omega_max = 300.0
    N_w = 4000
    Msp = 50

    ple = PowerLawExpSD(s, alpha, gamc)
    sdens = make_sdens(ple; scale=icm2ifs)
    sbeta = make_sbeta(ple, Temp)

    w = range(Omega_min, Omega_max, length=N_w) |> collect
    j = sbeta.(w)
    # Plot the spectral density and save the plot
    fig = Figure()
    ax = Axis(fig[1, 1],
    xlabel = "Frequency (cm^-1)",
    ylabel = "Spectral Density"
    )
    lines!(ax, w, j, label = "Spectral Density")
    axislegend(ax)
    save("spectral_density.png", fig)

    bcf = make_bcf(sdens, Temp, Ω_c)
    res = bsdo_discr(sbeta, Omega_min, Omega_max, N_w, Msp)
    wk = res.freq
    zk = res.coef
    calc_error(wk, zk, bcf, 500.0, 200)

end
