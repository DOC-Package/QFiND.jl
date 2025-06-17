include("../plot.jl")
using QFiND
using ExpFit

function esprit_decomp(bcf::Function, N_t::Integer, tc::Real, eps::Real) 
    ef = esprit(bcf, 0.0, tc, N_t, eps)
    return (expon = ef.expon / icm2ifs, coeff = ef.coeff / icm2ifs^2.0)
end

Ω = 1400.0
Γ = 200.0
λ = 600.0
Temp = 300.0
ub = 6000.0
T_max = 1000.0
N_t = 500
eps = 5e-3

sdens = BrownianSD(Ω, Γ, λ)
bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-6) # atol=1e-9

t = collect(range(0.0, T_max, length=N_t))
dt = t[2] - t[1]
c = bcf.(t)
ef = esprit(c, dt, eps)
println("dgree: ", size(ef.expon))
evaluate_error(t, ef.(t), c)
plot_bcf(t, ef.(t), c,  "./figure/bcf_esprit.png")
save_expon_coeff_union(ef.expon ./ icm2ifs, ef.coeff ./ icm2ifs^2.0, "expon_coeff_bo1400_l600_300K.txt")