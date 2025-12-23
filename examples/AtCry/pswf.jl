include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using ProlateSpheroidalWaveFunctions
using Statistics
using Serialization

Temp = 0.0
Omega_c = 4000.0
Omega_min = 0.0
Omega_max = 3700.0
N_w = 4000
T_c = 300.0
N_t = 500
eps = 1e-5
lam = 1.34 / icm2ev
r = open("r_atcry.bin", "r") do io
    deserialize(io)
end

sdens = AAAfittedSD(r, lam; ub=Omega_c)
E0 = reorganization_energy(sdens; ub=Omega_c)
println("Reorganization energy (eV): ", E0 * icm2ev)

Temp = 0.0
sbeta = BosonicQNSD(sdens, Temp)
f = w -> sbeta(w; scale=icm2ifs) / pi

# Correlation function via PSWF expansion
n_terms = 50
pswfft = pswf_expansion_fourier(f, Omega_min * icm2ifs, Omega_max * icm2ifs, T_c, n_terms; pm=-1.0)

sbeta_fit = pswfft.pswf

omega = range(Omega_min, Omega_max, length=2000)
fw = [f(w * icm2ifs) for w in omega]
fit = [sbeta_fit(w * icm2ifs) for w in omega]
plot_qnsd(omega, [fw, fit], "figure/pswf_ohmic_sbeta.png"; labels=["Exact", "PSWF fit"])

# Correlation function from analytical solution (for comparison)
bcf_exact = BosonicBCF(sdens, Temp; ub=2000.0, rtol=1e-6)

# Set time range
t_max = 250.0  # fs
n_points = 500
t_vals = range(0.0, t_max, length=n_points)

# Compute correlation functions
C_pswf = [pswfft(t) for t in t_vals]
C_exact = [bcf_exact(t) for t in t_vals]

# Compute errors
normalized_absolute_error = abs.(C_pswf - C_exact) ./ abs(C_exact[1])

# Display statistics
max_abs_error = maximum(normalized_absolute_error)
mean_abs_error = sum(normalized_absolute_error) / length(normalized_absolute_error)

println("Normalized maximum absolute error: $(max_abs_error)")
println("Normalized mean absolute error:    $(mean_abs_error)")
plot_bcf(t_vals, C_pswf, C_exact, "figure/pswf_ohmic_bcf.png")