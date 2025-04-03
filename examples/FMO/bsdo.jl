include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Serialization


Temp = 300.0
Ω_c = 495.0
Ω_min = -400.0
Ω_max = 500.0
N_w = 2000
T_c = 1500.0
N_t = 2000
M_sp = 70

r = open("r_fmo.bin", "r") do io
    deserialize(io)
end

sdens = AAAfittedSD(r)
plot_qnsd(sdens, 0.0, Ω_max, 2000, "fmo_aaa.png")
E0 = reorganization_energy(sdens; uplim=Ω_c)
println("Reorganization energy: ", E0)
E_r = 50.0

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)

res = bsdo_discr(sbeta, Ω_min, Ω_max, M_sp; n_lanczos=N_w)
freq = res.freq
coeff = res.coeff * π / 2.0
evaluate_error(freq, coeff, bcf, T_c, N_t)
plot_bcf(freq, coeff, bcf, T_c, N_t, "bcf_fmo_bsdo.png")


