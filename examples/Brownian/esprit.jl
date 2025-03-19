include("../plot.jl")
include("../save.jl")
using QFiND

Ω = 1400.0
Γ = 200.0
λ = 600.0
Temp = 300.0
Ω_c = 6000.0
T_c = 500.0
N_t = 1000
eps = 1e-2

sdens = BrownianSD(Ω, Γ, λ)
bcf = BosonicBCF(sdens, Temp, Ω_c)

res = esprit_decomp(bcf, N_t, T_c, eps) 
a = res.expon
c = res.coeff
println("degree of the exponential: ", size(a,1))
evaluate_error(a, c, bcf, T_c, N_t)
plot_bcf(a, c, bcf, T_c, N_t, "esprit.png")
save_array("bo.txt", c, a)

