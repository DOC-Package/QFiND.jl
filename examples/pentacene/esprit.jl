include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Statistics
using Serialization
using ExpFit

Temp = 300.0
ub = 2000.0
T_max = 500.0
N_t = 1000
eps = 1e-1

r = open("r_pen.bin", "r") do io
    deserialize(io)
end

E_reorg = 600.0
sdens = AAAfittedSD(r, E_reorg)
println("Reorganization energy: ", sdens.reorgene)
sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp; ub=ub)

t = collect(range(0.0, T_max, length=N_t))
dt = t[2] - t[1]
c = bcf.(t)

ef = esprit(c, dt, eps)
println("dgree: ", size(ef.expon))

evaluate_error(t, ef.(t), c)
plot_bcf(t, ef.(t), c,  "./figure/bcf_pen_esprit.png")
