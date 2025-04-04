include("../plot.jl")
include("../save.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Statistics
using Serialization

Temp = 0.0
Ω_c = 500.0
Ω_min = -500.0
Ω_max = 500.0
N_w = 2000
T_c = 1500.0
N_t = 2000
E_reorg = 35.0
eps = 1e-2

r = open("r_fmo.bin", "r") do io
    deserialize(io)
end

sdens = AAAfittedSD(r, E_reorg)
println("Reorganization energy: ", sdens.λ)

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)

res = id_discr(dataset, eps)
freq = res.freq
coeff = res.coeff
evaluate_error(freq, coeff, dataset.bcf, dataset.time)
plot_bcf(freq, coeff, dataset.bcf, dataset.time, "bcf_fmo_id.png")
save_freq_coeff(freq, sqrt.(coeff), "fmo_id.txt")
