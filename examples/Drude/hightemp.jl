include("../plot.jl")
include("../save.jl")
using QFiND

γ = 50.0
λ = 100.0
Temp = 300.0
Ω_c = 6000.0
T_c = 100.0
N_t = 1000
eps = 5e-3

sdens = DrudeSD(γ, λ)
dht = Drude_HighT(sdens, Temp; scale=icm2ifs)
bcf = t -> dht.coeff * exp(-dht.expon * t)


