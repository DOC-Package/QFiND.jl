include("../plot.jl")
using QFiND
using ExpFit

γ = 50.0
λ = 1200.0
Temp = 300.0
eps = 1e-2

sdens = DrudeSD(γ, λ)

T_max = 2000.0
N_t = 2000
t = collect(range(0.0, T_max, length=N_t))

npade = 5
ntpsd = 1
a, c = psd(sdens, Temp, npade; pade_type=:N)
at, ct = tpsd(sdens, Temp, npade, ntpsd; pade_type=:N)
bcf = bcf_approx(t, a, c)
bcf_t = bcf_approx(t, at, ct)
println("degree: ", size(at))
evaluate_error(t, bcf_t, bcf)
plot_bcf(t, bcf_t, bcf,  "./figure/bcf_psd.png")
#save_expon_coeff(a ./ icm2ifs, c ./ icm2ifs^2.0, "expon_coeff_d50_l150_50K_psd4.txt")
save_expon_coeff(at ./ icm2ifs, ct ./ icm2ifs^2.0, "expon_coeff_d50_l1200_300K_tpsd5-1.txt")