include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Statistics
using Serialization

Temp = 300.0
Ω_c = 4000.0
Ω_min = -1000.0
Ω_max = 4000.0
N_ω = 5000
T_max = 250.0
N_t = 500
eps = 1e-3
lam = 1.34 / icm2ev
r = open("r_atcry.bin", "r") do io
    deserialize(io)
end

sdens = AAAfittedSD(r, lam; ub=Ω_c)
E0 = reorganization_energy(sdens; ub=Ω_c)
println("Reorganization energy (eV): ", E0 * icm2ev)

sbeta = BosonicQNSD(sdens, Temp)
plot_qnsd(sbeta, Ω_min, Ω_max, N_ω, "./figure/sbeta_atcry_300K.png")
bcf = BosonicBCF(sdens, Temp, Ω_c)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t)

res = id_discr(dataset, eps)
ω = res.freq
g = res.coeff
evaluate_error(ω, g, dataset.bcf, dataset.time)
plot_bcf(ω, g, dataset.bcf, dataset.time, "./figure/bcf_atcry_300K_id.png")
