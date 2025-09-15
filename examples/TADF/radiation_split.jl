include("../plot.jl")
include("split_ri.jl")
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

bcf = BosonicBCF(sdens, Temp; ub=ub, rtol=1e-7)         # function for BCF
t = collect(range(0.0, T_max, length=N_t))              # time points
dt = t[2] - t[1]                                        # time step
println("Start BCF preparation")
c = bcf.(t)                                     
println("finished BCF preparation")

# ESPRIT fitting
ef = esprit(c, dt, eps)                       
println("dgree: ", size(ef.expon))

freq_real, coeff_real, freq_imag, coeff_imag = split_real_imag(ef.expon, ef.coeff)
save_expon_coeff(freq_real ./ icm2ifs, coeff_real ./ icm2ifs^2.0, "expon_coeff_radiation_real.txt")  
save_expon_coeff(freq_imag ./ icm2ifs, coeff_imag ./ icm2ifs^2.0, "expon_coeff_radiation_imag.txt")

bcf_real = t -> sum(coeff_real .* exp.(-freq_real .* t))
bcf_imag = t -> sum(coeff_imag .* exp.(-freq_imag .* t))
bcf_approx = t -> bcf_real(t) + im * bcf_imag(t)
evaluate_error(t, bcf_approx.(t), c)
plot_bcf(t, bcf_approx.(t), c,  "./figure/bcf_radiation.png")
