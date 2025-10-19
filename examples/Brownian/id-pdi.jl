include("../plot.jl")
using QFiND
using DelimitedFiles

A = readdlm("pdi-bo.txt")
Ω = A[:,1]
λ = A[:,2]
Γ = A[:,3]

Temp = 300.0
Ω_min = -600.0
Ω_max = 2000.0
N_ω = 2000
T_max = 500.0
N_t = 400
eps = 5e-3

sdens = BrownianSD(Ω, Γ, λ)
sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp; ub=6000.0, rtol=1e-8) # atol=1e-9

elapsed_time = @elapsed begin
    println("Started the Computing of the initial data...")
    dataset, errbcf = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t)
    println("Estimated maximum error in the BCF: ", errbcf)
end
println("Elapsed time for initial data preparation: ", elapsed_time, " seconds")
elapsed_time = @elapsed begin
    println("Computing the ID...")
    res = id_discr(dataset, eps)
end
println("Elapsed time for ID computation: ", elapsed_time, " seconds")
ω = res.freq
g = res.coeff
t = dataset.time
approx = bcf_approx.(t, Ref(ω), Ref(g))
evaluate_error(t, approx, dataset.bcf)
plot_bcf(t, approx, dataset.bcf, "./figure/bcf_id.png")
save_freq_coeff(ω, g, "freq_coeff_bo_pdi.txt")

plot_freq_coeff(sbeta, ω, g, Ω_min, Ω_max, N_ω, "./figure/brownian_id.png")


