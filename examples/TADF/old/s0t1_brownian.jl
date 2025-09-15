include("../plot.jl")
using LinearAlgebra
using QFiND
using ExpFit
using DelimitedFiles

data = readdlm("DATA_S0T1.txt")
col1 = data[:, 1]
col2 = data[:, 2] 
Ω = Float64.(col1)           # frequencies in cm^-1
λ = Float64.(col2)           # reorganization energies in cm^-1
λ[λ .< 0] .= 0.0 
Γ = ones(length(Ω)) .* 50.0  # width

# Define Brownian spectral density
sdens = BrownianSD(Ω, Γ, λ)

# Plot spectral density
omega = range(0.0, stop=3300.0, length=1000)
Jw = sdens.(omega)
plot_qnsd(omega, Jw, "./figure/Jw_s0t1_brownian.png")

# Parameters for BCF calculation and fitting
Temp = 300.0        # temperature in K
ub = 3500.0         # upper bound of integration frequency in cm^-1
T_max = 1500.0      # maximum time in fs
N_t = 500          # number of time points
eps = 5e-2          # fitting error tolerance or set degree directly

bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-7)         # function for BCF
t = collect(range(0.0, T_max, length=N_t))              # time points
dt = t[2] - t[1]                                        # time step
println("Start BCF preparation")
c = bcf.(t)                                     
println("finished BCF preparation")

# ESPRIT fitting
ef = esprit(c, dt, eps)                       
println("dgree: ", size(ef.expon))
save_expon_coeff(ef.expon ./ icm2ifs, ef.coeff ./ icm2ifs^2.0, "expon_coeff_s0t1.txt")  

evaluate_error(t, ef.(t), c)
plot_bcf(t, ef.(t), c,  "./figure/bcf_s0t1.png")
