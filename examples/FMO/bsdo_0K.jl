include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Serialization

Temp = 0.0
Ω_c = 500.0
Ω_min = 0.0
Ω_max = 500.0
N_ω = 4000

# reorganization energy
E_reorg = 35.0
# degree
degree = 100  

r = open("r_fmo.bin", "r") do io
    deserialize(io)
end

sdens = AAAfittedSD(r, E_reorg; ub=Ω_c)

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)

res = bsdo_discr(sbeta, Ω_min, Ω_max, degree; n_lanczos=N_ω)
ω = res.freq
g = res.coeff

# time range
T_max = 1500.0
N_t = 2000
t = collect(range(0.0, T_max, length=N_t))
reference = bcf.(t)
approx = bcf_approx.(t, Ref(ω), Ref(g))
evaluate_error(t, approx, reference)
plot_bcf(t, approx, reference, "./figure/bcf_fmo_bsdo_0K.png")


E_reorg = reorganization_energy(ω, g)
println("Effective reorganization energy: ", E_reorg)

