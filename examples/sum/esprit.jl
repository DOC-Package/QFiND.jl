include("../plot.jl")
using QFiND
using ExpFit

function esprit_decomp(bcf::Function, N_t::Integer, tc::Real, eps::Real) 
    ef = esprit(bcf, 0.0, tc, N_t, eps)
    return (expon = ef.expon / icm2ifs, coeff = ef.coeff / icm2ifs^2.0)
end

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

ub = 6000.0
eps = 5e-3

bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-6) # atol=1e-9

t = collect(range(0.0, T_max, length=N_t))
dt = t[2] - t[1]
c = bcf.(t)
ef = esprit(c, dt, eps)
println("dgree: ", size(ef.expon))
evaluate_error(t, ef.(t), c)
plot_bcf(t, ef.(t), c,  "./figure/bcf_esprit.png")
save_expon_coeff_union(ef.expon ./ icm2ifs, ef.coeff ./ icm2ifs^2.0, "expon_coeff_bo1400_l600_300K.txt")