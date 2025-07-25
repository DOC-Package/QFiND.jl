include("../plot.jl")
using QFiND
using ExpFit

function esprit_decomp(bcf::Function, N_t::Integer, tc::Real, eps::Real) 
    ef = esprit(bcf, 0.0, tc, N_t, eps)
    return (expon = ef.expon / icm2ifs, coeff = ef.coeff / icm2ifs^2.0)
end

# spectral density
s = 1.0
alpha = 50.0
gamc = 50.0
sdens = PowerLawExpSD(s, gamc; alpha=alpha)

Temp = 300.0
# time range
T_max = 4000.0
N_t = 200
eps = 1e-9

sbeta = BosonicQNSD(sdens, Temp)
bcf = BosonicBCF(sdens, Temp; ub=1000, rtol=1e-8) # atol=1e-9

t = collect(range(0.0, T_max, length=N_t))
dt = t[2] - t[1]
c = bcf.(t)
ef = esprit(c, dt, eps)
println("dgree: ", size(ef.expon))
evaluate_error(t, ef.(t), c)
plot_bcf(t, ef.(t), c,  "./figure/bcf_esprit.png")

ω = collect(range(-400.0, 400.0, length=1000))
expon = ef.expon ./ icm2ifs
coeff = ef.coeff ./ icm2ifs^2.0
sbeta_eff = EffectiveBosonicQNSD(expon, coeff)
plot_qnsd(ω, [sbeta.(ω), sbeta_eff.(ω)], "./figure/sbeta.png")
