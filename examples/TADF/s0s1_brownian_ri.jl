include("../plot.jl")
using LinearAlgebra
using QFiND
using ExpFit
using DelimitedFiles

data = readdlm("DATA_S0S1.txt")
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
plot_qnsd(omega, Jw, "./figure/Jw_s0s1_brownian.png")

# Parameters for BCF calculation and fitting
Temp = 300.0        # temperature in K
ub = 3500.0         # upper bound of integration frequency in cm^-1
T_max = 1500.0      # maximum time in fs
N_t = 1000           # number of time points
eps = 5e-2          # fitting error tolerance or set degree directly

bcf_real = BosonicBCF_Real(sdens, Temp; ub=ub, rtol=1e-7)         # function for BCF
bcf_imag = BosonicBCF_Imag(sdens, Temp; ub=ub, rtol=1e-7)         # function for BCF
t = collect(range(0.0, T_max, length=N_t))              # time points
dt = t[2] - t[1]                                        # time step
println("Start BCF preparation")
cr = bcf_real.(t) 
ci = bcf_imag.(t)   
c = cr + im * ci                                 
println("finished BCF preparation")

# ESPRIT fitting
efr = esprit(cr, dt, eps)
efi = esprit(ci, dt, eps)    
println("dgree (real): ", size(efr.expon))
println("dgree (imag): ", size(efi.expon))
save_expon_coeff(efr.expon ./ icm2ifs, efr.coeff ./ icm2ifs^2.0, "expon_coeff_s0s1_real.txt")  
save_expon_coeff(efi.expon ./ icm2ifs, efi.coeff ./ icm2ifs^2.0, "expon_coeff_s0s1_imag.txt")

approx = efr.(t) + im * efi.(t)
evaluate_error(t, approx, c)
plot_bcf(t, approx, c,  "./figure/bcf_s0s1_real.png")
