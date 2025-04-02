include("smooth.jl")
include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation

Temp = 300.0
Ω_c = 490.0
Ω_min = -500.0
Ω_max = 500.0
N_w = 2000
T_c = 1000.0
N_t = 500
eps = 1e-4

data = readdlm("jw-fmo.txt")

residues = data[:, 1] + 1im * data[:, 2]
poles = data[:, 3] + 1im * data[:, 4]
residues = residues .* 35.0 ./ 42.734770
sdens = RationalSD(poles, residues) 
plot_qnsd(sdens, 0.0, Ω_max, 2000, "fmo_aaa.png")

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)
"finished preparing initial data"

res = id_discr(dataset, eps)
freq = res.freq
coeff = res.coeff
evaluate_error(freq, coeff, bcf, T_c, N_t)
plot_bcf(freq, coeff, bcf, T_c, N_t, "bcf_fmo_id.png")

