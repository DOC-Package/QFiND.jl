include("../plot.jl")
using QFiND
using ExpFit

Ω = 400.0
Γ = 80.0
λ = 150.0
Temp = 10.0
T_max = 200.0
N_t = 400
eps = 1e-2
ub = 6000.0

sdens = BrownianSD(Ω, Γ, λ)
bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-6) # atol=1e-9

t = collect(range(0.0, T_max, length=N_t))
bcf = bcf.(t)

a, c = spectral_decomposition(sdens, Temp; npade=5)
asd = a[1:2]
csd = c[1:2]
alt = a[3:end] 
clt = c[3:end]
er = expred(alt, clt, eps)
alt = er.expon
clt = er.coeff
a = [asd; alt]
c = [csd; clt]
approx = bcf_approx(t, a, c) 
println("degree: ", size(a))
evaluate_error(t, approx, bcf)
plot_bcf(t, approx, bcf,  "./figure/bcf_psd.png")
save_expon_coeff_union(a, c, "expon_coeff_bo1400_l600_300K.txt")