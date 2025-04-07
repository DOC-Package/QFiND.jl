include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Statistics
using Serialization
using ExpFit

Temp = 300.0
Ω_c = 4000.0
T_max = 250.0
N_t = 250
eps = 5e-2

r = open("r_atcry.bin", "r") do io
    deserialize(io)
end

E_reorg = 1.34 / icm2ev
sdens = AAAfittedSD(r, E_reorg)
println("Reorganization energy (eV): ", sdens.reorgene * icm2ev)
sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)

println("Preparing initial data")
t = collect(range(0.0, T_max, length=N_t))
dt = t[2] - t[1]
c = bcf.(t)

println("ESPRIT started")
ef = esprit(c, dt, eps)
println("dgree: ", size(ef.expon))

evaluate_error(t, ef.(t), c)
plot_bcf(t, ef.(t), c,  "./figure/bcf_atcry_esprit.png")
