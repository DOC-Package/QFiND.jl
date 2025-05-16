include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Serialization

Temp = 300.0
Ω_c = 500.0
Ω_min = -350.0
Ω_max = 490.0
N_ω = 2000
T_max = 1500.0
N_t = 2000
M_sp = 64
E_reorg = 35.0  

r = open("r_fmo.bin", "r") do io
    deserialize(io)
end

sdens = AAAfittedSD(r, E_reorg; ub=Ω_c)
plot_qnsd(sdens, 0.0, Ω_max, 2000, "fmo_aaa.png")
E0 = reorganization_energy(sdens; ub=Ω_c)
println("Reorganization energy: ", E0)

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)
dataset = InitialData(DiscrBSDO(), sbeta, Ω_min, Ω_max; n_lanczos=N_ω)

res = bsdo_discr(dataset, M_sp)
ω = res.freq
g = res.coeff

t = collect(range(0.0, T_max, length=N_t))
bcf_t = bcf.(t)
approx = bcf_approx.(t, Ref(ω), Ref(g))
evaluate_error(t, approx, bcf_t)
plot_bcf(t, approx, bcf_t, "./figure/bcf_fmo_bsdo.png")
save_freq_coeff(ω, g, "fmo_bsdo.txt")
plot_freq_coeff(sbeta, ω, g, Ω_min, Ω_max, N_ω, "./figure/fmo_bsdo.png")

