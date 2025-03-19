include("../plot.jl")
include("../save.jl")
using QFiND

Ω = 100.0
Γ = 100.0
λ = 200.0
Temp = 300.0
Ω_c = 6000.0
T_c = 100.0
N_t = 1000
eps = 5e-3

sdens = BrownianSD(Ω, Γ, λ)
bcf = BosonicBCF(sdens, Temp, Ω_c)

res = esprit_decomp(bcf, N_t, T_c, 2) 
a = res.expon
c = res.coeff
println("degree of the exponential: ", size(a,1))
evaluate_error(a, c, bcf, T_c, N_t)
plot_bcf(a, c, bcf, T_c, N_t, "esprit.png")
save_array_union("bo-lf-esprit.txt", c, a)

