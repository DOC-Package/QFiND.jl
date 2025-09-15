include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Statistics
using Serialization
using ExpFit

# Super-Ohmic spectral density
s = 3.0                  # super-Ohmic
lam = 20.0               # reorganization energy in cm^-1
gamc = 2.0 / icm2ev      # cutoff frequency in eV, converted to cm^-1
sdens = PowerLawExpSD(s, gamc; reorgene=lam)
println("Reorganization energy (cm^-1): ", sdens.reorgene)

Temp = 6000.0   # temperature in K

# fitting parameters
ub = 1e7          # upper bound of integration frequency in cm^-1
T_max = 5.0       # maximum time in fs
N_t = 200         # number of time points
eps = 7e-2        # fitting error tolerance or set degree directly
deg = 3

bcf_real = BosonicBCF_Real(sdens, Temp; ub=ub)  # function for real part of BCF
bcf_imag = BosonicBCF_Imag(sdens, Temp; ub=ub)  # function for imaginary part of BCF

t = collect(range(0.0, T_max, length=N_t))      # time points
dt = t[2] - t[1]                                # time step                               
cr = bcf_real.(t)
ci = bcf_imag.(t)
c = cr + im * ci

efr = esprit(cr, dt, eps)
efi = esprit(ci, dt, eps)
#efr = esprit(cr, dt, deg)
#efi = esprit(ci, dt, deg)

println("degree (real): ", size(efr.expon))
println("degree (imag): ", size(efi.expon))
save_expon_coeff(efr.expon ./ icm2ifs, efr.coeff ./ icm2ifs^2.0, "expon_coeff_radiation_real.txt")  
save_expon_coeff(efi.expon ./ icm2ifs, efi.coeff ./ icm2ifs^2.0, "expon_coeff_radiation_imag.txt")

ef = efr.(t) + im * efi.(t)
evaluate_error(t, ef, c)
plot_bcf(t, ef, c,  "./figure/bcf_radiation.png")
