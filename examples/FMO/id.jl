include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Statistics
using Serialization

Temp = 300.0
Ω_c = 495.0
Ω_min = -500.0
Ω_max = 500.0
N_ω = 2000
T_c = 1500.0
N_t = 2000
eps = 7e-2

r = open("r_fmo.bin", "r") do io
    deserialize(io)
end

E_reorg = 35.0
sdens = AAAfittedSD(r, E_reorg)
println("Reorganization energy: ", sdens.reorgene)

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_ω, n_time=N_t)

res = id_discr(dataset, eps)
ω = res.freq
g = res.coeff
evaluate_error(ω, g, dataset.bcf, dataset.time)
plot_bcf(ω, g, dataset.bcf, dataset.time, "./figure/bcf_fmo_id.png")
save_freq_coeff(ω, g, "fmo_id.txt")

plot_freq_coeff(sbeta, ω, g, Ω_min, Ω_max, N_ω, "./figure/fmo_id.png")
