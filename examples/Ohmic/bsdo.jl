include("../plot.jl")
using Test
using QFiND
  
# spectral density
s = 1.0
alpha = 50.0
gamc = 50.0
sdens = PowerLawExpSD(s, gamc; alpha=alpha)

# 300K 
# frequency range
Ω_min = -300.0
Ω_max = 400.0
# degree
degree = 60
# time range
T_max = 500.0
N_t = 200

Temp = 300.0
sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp)
res = bsdo_discr(sbeta, Ω_min, Ω_max, degree)
ω = res.freq
g = res.coeff
@test !isnothing(res)
evaluate_error(ω, g, bcf, T_max, N_t)
plot_bcf(ω, g, bcf, T_max, N_t, "bcf_bsdo.png")
save_freq_coeff(ω, g, "freq_coeff_bsdo.txt")

# 0K 
Temp = 0.0
Ω_min = 0.0
Ω_max = 500.0
# degree
degree = 60
# time range
T_max = 500.0
N_t = 200

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp)
res = bsdo_discr(sbeta, Ω_min, Ω_max, degree)
ω = res.freq
g = res.coeff
@test !isnothing(res)
evaluate_error(ω, g, bcf, T_max, N_t)
plot_bcf(ω, g, bcf, T_max, N_t, "bcf_0K_bsdo.png")
save_freq_coeff(ω, g, "freq_coeff_0K_bsdo.txt")

E_reorg = reorganization_energy(ω, g)
println("Effective reorganization energy: ", E_reorg)
