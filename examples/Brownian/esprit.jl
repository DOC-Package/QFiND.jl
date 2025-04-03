include("../plot.jl")
include("../save.jl")
using QFiND

function esprit_decomp(bcf::Function, N_t::Integer, tc::Real, eps::Real) 
    ef = esprit(bcf, 0.0, tc, N_t, eps)
    return (expon = ef.expon / icm2ifs, coeff = ef.coeff / icm2ifs^2.0)
end

Ω = 300.0
Γ = 100.0
λ = 50.0
Temp = 300.0
Ω_c = 6000.0
T_c = 100.0
N_t = 1000
eps = 5e-3

sdens = BrownianSD(Ω, Γ, λ)
bcf = BosonicBCF(sdens, Temp, Ω_c)

res = esprit_decomp(bcf, N_t, T_c, 2) 
a = res.expon
c = res.coeff
println("degree of the exponential: ", size(a,1))
evaluate_error(a, c, bcf, T_c, N_t)
plot_bcf(a, c, bcf, T_c, N_t, "esprit.png")
save_array_union("bo_esprit_300_2.txt", c, a)

