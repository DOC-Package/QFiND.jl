include("../plot.jl")
using QFiND

Temp = 50.0
eps = 1e-2

Ω = 1400.0
Γ = 200.0
λ = 570.0
sd1 = BrownianSD(Ω, Γ, λ)

γ = 50
λ = 30.0
sd2 = DrudeSD(γ, λ)
sdens = sd1 + sd2

T_max = 2000.0
N_t = 2000
t = collect(range(0.0, T_max, length=N_t))

npade=11
a, c = psd(sdens, Temp, npade; pade_type=:Nm1)
at, ct = tpsd(sdens, Temp, npade, 1; pade_type=:Nm1)
bcf = bcf_approx(t, a, c) 
bcf_t = bcf_approx(t, at, ct)
println("degree: ", size(at))
evaluate_error(t, bcf_t, bcf)
plot_bcf(t, bcf_t, bcf,  "./figure/bcf_psd.png")
#save_expon_coeff(a ./ icm2ifs, c ./ icm2ifs^2.0, "expon_coeff_d50_l150_50K_psd11.txt")
save_expon_coeff(at ./ icm2ifs, ct ./ icm2ifs^2.0, "expon_coeff_bo1400_l570_d50_l30_50K_tpsd1_Nm1.txt")