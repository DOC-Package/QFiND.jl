include("../plot.jl")
using QFiND
using ExpFit

γ = 50.0
λ = 30.0
Temp = 400.0
eps = 1e-2

sdens = DrudeSD(γ, λ)

T_max = 2000.0
N_t = 2000
t = collect(range(0.0, T_max, length=N_t))

npade=11
a, c = psd(sdens, Temp, npade; pade_type=:N)
at, ct = tpsd(sdens, Temp, npade, 1; pade_type=:N)
bcf = bcf_approx(t, a, c) 
bcf_t = bcf_approx(t, at, ct)
println("degree: ", size(at))
evaluate_error(t, bcf_t, bcf)
plot_bcf(t, bcf_t, bcf,  "./figure/bcf_psd.png")
save_expon_coeff(a ./ icm2ifs, c ./ icm2ifs^2.0, "expon_coeff_d50_l150_50K_psd11.txt")
save_expon_coeff(at ./ icm2ifs, ct ./ icm2ifs^2.0, "expon_coeff_d50_l30_50K_psd1.txt")