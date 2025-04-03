include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Statistics

Temp = 300.0
Ω_c = 495.0
Ω_min = -500.0
Ω_max = 500.0
N_w = 2000
T_c = 1500.0
N_t = 2000
eps = 3e-2

data = readdlm("fmo_smoothed.txt")
col1 = data[:, 1]
col2 = data[:, 2] 
ω = Float64.(col1) 
J = Float64.(col2)
J[J .< 0] .= 0.0 

println("AAA started")
r = aaa(ω, J; tol=1e-14, max_degree=1200, lookahead=200)
#r = mylawson(ω, J, r, 10)
println("degree: ", length(r.nodes))
err = norm(r.(ω) - J)
println("error: ", err)

sdens = AAAfittedSD(r)
plot_qnsd(sdens, 0.0, Ω_max, 2000, "fmo_aaa.png")
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
plot_bcf(freq, coeff, dataset.bcf, dataset.time, "bcf_fmo_id.png")
