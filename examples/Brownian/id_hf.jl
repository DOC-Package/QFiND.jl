include("../plot.jl")
using QFiND

Ω = 1400.0
Γ = 200.0
λ = 600.0
Temp = 300.0
Ω_min = -2000.0
Ω_max = 3000.0
N_ω = 2000
T_max = 2000.0
N_t = 500
eps = 1e-2

sdens = BrownianSD(Ω, Γ, λ)
sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, 6000.0; rtol=1e-6) # atol=1e-9
dataset, _ = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t)
res = id_discr(dataset, eps)
ω = res.freq
g = res.coeff
t = dataset.time
approx = bcf_approx.(t, Ref(ω), Ref(g))
evaluate_error(t, approx, dataset.bcf)
plot_bcf(t, approx, dataset.bcf, "./figure/bcf_id.png")
save_freq_coeff(ω, g, "freq_coeff_bo50.txt")

plot_freq_coeff(sbeta, ω, g, Ω_min, Ω_max, N_ω, "./figure/brownian_id.png")
