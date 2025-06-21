include("../plot.jl")
using QFiND

Ω = 50.0
Γ = 20.0
λ = 30.0
Temp = 10.0
T_max = 4000.0
N_t = 4000
eps = 1e-2
ub = 6000.0

sdens = BrownianSD(Ω, Γ, λ)
bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-6) # atol=1e-9

t = collect(range(0.0, T_max, length=N_t))
bcf = bcf.(t)

npade = 1
a, c = psd(sdens, Temp, npade)
println("a: ", a)
println("c: ", c)
approx = bcf_approx(t, a, c) 
evaluate_error(t, approx, bcf)
plot_bcf(t, approx, bcf,  "./figure/bcf_psd.png")
save_expon_coeff(a, c, "expon_coeff_bo50_l30_10K.txt")