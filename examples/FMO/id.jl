include("smooth.jl")
include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation

Temp = 300.0
Ω_c = 500.0
Ω_min = -500.0
Ω_max = 500.0
N_w = 2000
T_c = 1000.0
N_t = 1000
eps = 1e-4

data = readdlm("fmo_smoothed.txt")
col1 = data[:, 1]
col2 = data[:, 2] #.* data[:, 1] .^ 2.0
ω = Float64.(col1) 
J = Float64.(col2) 
#S = smoothing(ω, J)

bary = aaa(ω, J; tol=1e-4*maximum(J))
check(bary)
sdens = AAAfittedSD(bary)
plot_qnsd(sdens, 0.0, Ω_max * icm2ifs, 2000, "fmo_aaa.png")
E0 = reorganization_energy(sdens; uplim=Ω_c)
println("Reorganization energy: ", E0)
E_r = 50.0

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)

res = id_discr(sbeta, bcf, Ω_min, Ω_max, T_c, N_w, N_t, eps)
freq = res.freq
coeff = res.coeff
evaluate_error(freq, coeff, bcf, T_c, N_t)
plot_bcf(freq, coeff, bcf, T_c, N_t, "bcf_fmo_id.png")

