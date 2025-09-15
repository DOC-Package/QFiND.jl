include("../plot.jl")
using QFiND

Ω = 1400.0
Γ = 200.0
λ = 600.0
Temp = 50.0
T_max = 400.0
N_t = 400
eps = 1e-2
ub = 6000.0

sdens = BrownianSD(Ω, Γ, λ)
bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-6) # atol=1e-9

t = collect(range(0.0, T_max, length=N_t))
bcf = bcf.(t)

npade=4
a, c = tpsd(sdens, Temp, npade, 1; pade_type=:Nm1)
approx = bcf_approx(t, a, c) 
println("degree: ", size(a))
evaluate_error(t, approx, bcf)
plot_bcf(t, approx, bcf,  "./figure/bcf_psd.png")
save_expon_coeff_union(a ./ icm2ifs, c ./ icm2ifs^2.0, "expon_coeff_bo1400_l600_50K.txt")