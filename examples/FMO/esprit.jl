include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Statistics
using Serialization
using ExpFit

Temp = 300.0
Ω_c = 495.0
T_max = 1500.0
N_t = 2000
eps = 1e-1

r = open("r_fmo.bin", "r") do io
    deserialize(io)
end

E_reorg = 35.0
sdens = AAAfittedSD(r, E_reorg)
println("Reorganization energy: ", sdens.reorgene)
sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)

t = collect(range(0.0, T_max, length=N_t))
dt = t[2] - t[1]
c = bcf.(t)

ef = esprit(c, dt, eps)
println("dgree: ", size(ef.expon))

evaluate_error(t, ef.(t), c)
plot_bcf(t, ef.(t), c,  "./figure/bcf_fmo_esprit.png")
