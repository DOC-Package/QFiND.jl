include("../plot.jl")
using Test
using QFiND

# spectral density
s = 1.0
alpha = 50.0
gamc = 50.0
sdens = PowerLawExpSD(s, gamc; alpha=alpha)

# 300 K
Temp = 300.0
Ω_min = -400.0
Ω_max = 400.0
N_ω = 2000
T_max = 500.0
N_t = 200
eps = 1e-3

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t)
res = id_discr(dataset, eps)
ω = res.freq
g = res.coeff
evaluate_error(ω, g, bcf, T_max, N_t)
plot_bcf(ω, g, bcf, dataset.time, "bcf_id.png")
save_freq_coeff(ω, g, "freq_coeff_id.txt")

# 0 K
Temp = 0.0
Ω_min = 0.0
Ω_max = 600.0
N_ω = 2000
T_max = 500.0
N_t = 200
eps = 1e-3

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t)
res = id_discr(dataset, eps)
ω = res.freq
g = res.coeff
evaluate_error(ω, g, bcf, T_max, N_t)
plot_bcf(ω, g, bcf, dataset.time, "bcf_0K_id.png")
save_freq_coeff(ω, g, "freq_coeff_0K_id.txt")

E_reorg = reorganization_energy(ω, g)
println("Effective reorganization energy: ", E_reorg)