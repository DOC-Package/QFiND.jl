include("../plot.jl")
using LinearAlgebra
using QFiND
using ExpFit
using DelimitedFiles

data = readdlm("DATA_S0S1.txt")
col1 = data[:, 1]
col2 = data[:, 2] 
ω = Float64.(col1) 
J = Float64.(col2)
J[J .< 0] .= 0.0 

Temp = 300.0
freq, coeff = BosonicQNSD_Discrete(ω, J, Temp)
println("frequencies: ", freq)
println("coefficients: ", coeff)
bcf = t -> bcf_discrete(t, freq, coeff)

T_max = 10000.0      # maximum time in fs
N_t = 5000          # number of time points
eps = 1e-2          # fitting error tolerance or set degree directly

t = collect(range(0.0, T_max, length=N_t))
dt = t[2] - t[1]
c = bcf.(t)

ef = esprit(c, dt, eps)
println("dgree: ", size(ef.expon))

evaluate_error(t, ef.(t), c)
plot_bcf(t, ef.(t), c,  "./figure/bcf_s0s1.png")
