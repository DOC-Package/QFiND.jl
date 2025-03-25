include("../plot.jl")
include("../save.jl")
using QFiND

Ω = 400.0
Γ = 200.0
λ = 600.0
Temp = 300.0
Ω_c = 6000.0
Ω_min = -3000.0
Ω_max = 3000.0
N_w = 4000
T_c = 200.0
N_t = 500
eps = 1e-3

sdens = BrownianSD(Ω, Γ, λ)
sbeta = BosonicQNSD_HighT(sdens, Temp)
bht = BrownianHighT(sdens, Temp; scale=icm2ifs)
bcf = t -> sumexp(t, bht.expon, bht.coeff)

#println("coeff", bht.coeff ./ icm2ifs^2.0)
#println("expon", bht.expon ./ icm2ifs) 

save_array_unite("bo_ht_hf.txt", bht.coeff ./ icm2ifs^2.0, bht.expon ./ icm2ifs)

dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)
res = id_discr(dataset, eps)

freq = res.freq
coef = res.coef
evaluate_error(freq, coef, bcf, T_c, N_t)
plot_bcf(freq, coef, bcf, T_c, N_t, "id_ht.png")

save_array("omega_g_ht_hf.txt", freq, sqrt.(coef))

