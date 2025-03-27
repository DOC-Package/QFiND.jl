include("../plot.jl")
using LinearAlgebra
using KissSmoothing
using DelimitedFiles

function smoothing(ω, J)
    idx = findfirst(ω .> 150)
    J1 = J[1:idx]
    J2 = J[idx+1:end]
    plot_qnsd(ω, J, "fmo.png")
    #
    S1, N1 = denoise(J1; factor=0.95)
    S2, N2 = denoise(J2; factor=0.5)
    println("N1: ", norm(N1))
    println("N2: ", norm(N2))
    # 
    S = vcat(S1, S2)
    plot_qnsd(ω, S, "fmo_smoothed.png")
    return S
end


