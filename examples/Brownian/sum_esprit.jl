include("../plot.jl")
using QFiND
using ExpFit

Ω = 1400.0
Γ = 200.0
λ = 450.0
sdens1 = BrownianSD(Ω, Γ, λ)

Ω = 50.0
Γ = 10.0
λ = 150.0
sdens2 = BrownianSD(Ω, Γ, λ)

Temp = 300.0
ub = 10000.0
T_max = 2000.0
N_t = 2000
eps = 1e-2

sdens = sdens1 + sdens2
bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-6) # atol=1e-9

t = collect(range(0.0, T_max, length=N_t))
dt = t[2] - t[1]
c = bcf.(t)
ef = expfit(c, dt, eps)
println("dgree: ", size(ef.expon))
evaluate_error(t, ef.(t), c)
plot_bcf(t, ef.(t), c,  "./figure/bcf_esprit.png")
save_expon_coeff_union(ef.expon ./ icm2ifs, ef.coeff ./ icm2ifs^2.0, "expon_coeff_bo1400_l600_300K.txt")
