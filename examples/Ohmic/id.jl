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
Ω_min = -500.0
Ω_max = 500.0
N_ω = 2000
T_max = 1000.0
N_t = 500
eps = 1e-6

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, 1000)
dataset, _ = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t)

res = id_discr(dataset, eps)
ω = res.freq
g = res.coeff
t = dataset.time
approx = bcf_approx.(t, Ref(ω), Ref(g))
evaluate_error(t, approx, dataset.bcf)
plot_bcf(t, approx, dataset.bcf, "bcf_id.png")
save_freq_coeff(ω, g, "freq_coeff_id.txt")
plot_freq_coeff(sbeta, ω, g, Ω_min, Ω_max, N_ω, "./figure/ohmic_id.png")

# 0 K
Temp = 0.0
Ω_min = 0.0
Ω_max = 500.0
N_ω = 2000
T_max = 1000.0
N_t = 200
eps = 1e-6

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp)
dataset, _ = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t)
res = id_discr(dataset, eps)
ω = res.freq
g = res.coeff
t = dataset.time
approx = bcf_approx.(t, Ref(ω), Ref(g))
save_freq_coeff(ω, g, "freq_coeff_0K_id.txt")
evaluate_error(t, approx, dataset.bcf)
plot_bcf(t, approx, dataset.bcf, "bcf_0K_id.png")
plot_freq_coeff(sbeta, ω, g, Ω_min, Ω_max, N_ω, "./figure/ohmic_id_0K.png")

E_reorg = reorganization_energy(ω, g)
println("Effective reorganization energy: ", E_reorg)

