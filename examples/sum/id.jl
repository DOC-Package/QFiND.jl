include("../plot.jl")
using Test
using QFiND

# spectral density 
# sdens1
s1 = 1.0
alpha1 = 50.0
gamc1 = 50.0
sdens1 = PowerLawExpSD(s1, gamc1; alpha=alpha1)
#sdens2
s2 = 2.0
alpha2 = 100.0
gamc2 = 100.0
sdens2 = PowerLawExpSD(s2, gamc2; alpha=alpha2)

# combined spectral density
sdens = sdens1 + sdens2         # sdens = SumSD(sdens1, sdens2) 

# Parameters
Temp = 300.0
Ω_min = -700.0
Ω_max = 1000.0
N_ω = 2000
T_max = 500.0
N_t = 500
eps = 1e-2

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, 1000)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t)
res = id_discr(dataset, eps)
ω = res.freq
g = res.coeff
t = dataset.time
approx = bcf_approx.(t, Ref(ω), Ref(g))
evaluate_error(t, approx, dataset.bcf)
save_freq_coeff(ω, g, "freq_coeff_id.txt")
plot_bcf(t, approx, dataset.bcf, "./figure/bcf_id.png")
plot_freq_coeff(sbeta, ω, g, Ω_min, Ω_max, N_ω, "./figure/sum_id.png")