include("../plot.jl")
using QFiND

Ω = [30.0, 50.0, 100.0, 200.0]
Γ = [10.0, 20.0, 30.0, 40.0]                                                
λ = [30.0, 50.0, 100.0, 200.0]
Temp = 300.0
Ω_min = -600.0
Ω_max = 1100.0
N_ω = 2000
T_max = 200.0
N_t = 400
eps = 1e-3

sdens = BrownianSD(Ω, Γ, λ)
sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t)
res = id_discr(dataset, eps)
ω = res.freq
g = res.coeff
t = dataset.time
approx = bcf_approx.(t, Ref(ω), Ref(g))
evaluate_error(t, approx, dataset.bcf)
plot_bcf(t, approx, dataset.bcf, "./figure/bcf_id.png")
save_freq_coeff(ω, g, "freq_coeff_bo50.txt")

plot_freq_coeff(sbeta, ω, g, Ω_min, Ω_max, N_ω, "./figure/brownian_id.png")
