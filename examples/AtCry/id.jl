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
eps = 1e-5
lam = 1.34 / icm2ev
r = open("r_atcry.bin", "r") do io
    deserialize(io)
end

sdens = AAAfittedSD(r, lam; ub=Ω_c)
E0 = reorganization_energy(sdens; ub=Ω_c)
println("Reorganization energy (eV): ", E0 * icm2ev)

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)
println("Preparing initial data")
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)

println("ID started")
res = id_discr(dataset, eps)
ω = res.freq
g = res.coeff
evaluate_error(ω, g, dataset.bcf, dataset.time)
plot_bcf(ω, g, dataset.bcf, dataset.time, "./figure/bcf_atcry_id.png")

E_reorg = reorganization_energy(ω, g)
println("Effective reorganization energy: ", E_reorg * icm2ev)
