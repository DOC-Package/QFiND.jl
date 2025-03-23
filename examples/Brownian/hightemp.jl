include("../plot.jl")
include("../save.jl")
using QFiND

Ω = 100.0
Γ = 100.0
λ = 50.0
Temp = 300.0
Ω_c = 6000.0
Ω_min = -600.0
Ω_max = 600.0
N_w = 2000
T_c = 150.0
N_t = 100
eps = 1e-8

sdens = BrownianSD(Ω, Γ, λ)
sbeta = BosonicQNSD_HighT(sdens, Temp)
bht = BrownianHighT(sdens, Temp; scale=icm2ifs)
bcf = t -> sumexp(t, bht.expon, bht.coeff)

dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)
res = id_discr(dataset, eps)

freq = res.freq
coef = res.coef
evaluate_error(freq, coef, bcf, T_c, N_t)
plot_bcf(freq, coef, bcf, T_c, N_t, "id_ht.png")

save_array("omega_g_ht.txt", freq, sqrt.(coef))

