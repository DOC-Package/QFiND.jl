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
T_max = 250.0
N_t = 500
degree = 50

lam = 1.34 / icm2ev
r = open("r_atcry.bin", "r") do io
    deserialize(io)
end

sdens = AAAfittedSD(r, lam; ub=Ω_c)
E0 = reorganization_energy(sdens; ub=Ω_c)
println("Reorganization energy (eV): ", E0 * icm2ev)

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)
res = bsdo_discr(sbeta, Ω_min, Ω_max, degree; n_lanczos=N_w)
ω = res.freq
g = res.coeff

# time range
t = collect(range(0.0, T_max, length=N_t))
bcf_t = bcf.(t)
evaluate_error(ω, g, bcf_t, t)
plot_bcf(ω, g, bcf_t, t, "./figure/bcf_atcry_bsdo.png")

E_reorg = reorganization_energy(ω, g)
println("Effective reorganization energy: ", E_reorg * icm2ev)
