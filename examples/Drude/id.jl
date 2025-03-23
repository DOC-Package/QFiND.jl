include("../plot.jl")
include("../save.jl")
using QFiND

    
γ = 50.0
λ = 100.0
Temp = 300.0
Ω_c = 6000.0
Ω_min = -1000.0
Ω_max = 1000.0
N_w = 2000
T_c = 150.0
N_t = 100
eps = 1e-8

sdens = BrownianSD(γ, λ)
sbeta = BosonicQNSD_HighT(sdens, Temp)
dht = Drude_HighT(sdens, Temp; scale=icm2ifs)
bcf = t -> dht.coeff * exp(-dht.expon * t)

dataset = InitialData(DiscrID(), sbeta, bcf, Ω_min, Ω_max, T_c; n_freq=N_w, n_time=N_t)
res = id_discr(dataset, eps)

freq = res.freq
coef = res.coef
evaluate_error(freq, coef, bcf, T_c, N_t)
plot_bcf(freq, coef, bcf, T_c, N_t, "id.png")

save_array("omega_g_bo_l100.txt", freq, sqrt.(coef))

