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
N_w = 2000
T_c = 1500.0
N_t = 2000
eps = 6e-2

r = open("r_fmo.bin", "r") do io
    deserialize(io)
end

sdens = AAAfittedSD(r)
E0 = reorganization_energy(sdens; ub=Ω_c)
println("Reorganization energy: ", E0)
E_r = 35.0
factor = E_r / E0

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)

res = id_discr(dataset, eps)
freq = res.freq
coeff = res.coeff
evaluate_error(freq, coeff, dataset.bcf, dataset.time)
plot_bcf(freq, coeff, dataset.bcf, dataset.time, "bcf_fmo_id.png")
