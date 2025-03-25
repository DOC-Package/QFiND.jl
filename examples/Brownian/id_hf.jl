include("../plot.jl")
include("../save.jl")
using QFiND

Ω = 300.0
Γ = 100.0
λ = 100.0
Temp = 300.0
Ω_c = 6000.0
Ω_min = -1000.0
Ω_max = 1000.0
N_w = 4000
T_c = 200.0
N_t = 500
eps = 5e-3

sdens = BrownianSD(Ω, Γ, λ)
sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp, Ω_c)
dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)
res = id_discr(dataset, eps)

freq = res.freq
coef = res.coef
evaluate_error(freq, coef, bcf, T_c, N_t)
plot_bcf(freq, coef, bcf, T_c, N_t, "id_hf.png")

save_array("omega_g_bo_hf.txt", freq, sqrt.(coef))

