include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Statistics
using Serialization

Temp = 0.0
Ω_c = 4000.0
Ω_min = 0.0
Ω_max = 4000.0
N_w = 4000
T_c = 250.0
N_t = 500
eps = 2e-2

r = open("r_atcry.bin", "r") do io
    deserialize(io)
end

sdens = AAAfittedSD(r)
E0 = reorganization_energy(sdens; uplim=Ω_c)
println("Reorganization energy: ", E0)
E_r = 50.0

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)

res = id_discr(dataset, eps)
freq = res.freq
coeff = res.coeff
evaluate_error(freq, coeff, dataset.bcf, dataset.time)
plot_bcf(freq, coeff, dataset.bcf, dataset.time, "bcf_atcry_id.png")
