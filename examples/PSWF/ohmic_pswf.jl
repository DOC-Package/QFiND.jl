include("../plot.jl")
using QFiND
using ProlateSpheroidalWaveFunctions

# Define spectral density (Ohmic)
s = 1.0
gam = 50.0
lam = 50.0
Omega_min = -180.0
Omega_max = 250.0

sd = PowerLawExpSD(s, gam; reorgene=lam)
Temp = 300.0
sbeta = BosonicQNSD(sd, Temp)
f = w -> sbeta(w; scale=icm2ifs) / pi

# Correlation function via PSWF expansion
n_terms = 15
T = 1200.0  # Time duration parameter
pswfft = pswf_expansion_fourier(f, Omega_min * icm2ifs, Omega_max * icm2ifs, T, n_terms; pm=-1.0)

sbeta_fit = pswfft.pswf

omega = range(Omega_min, Omega_max, length=2000)
fw = [f(w * icm2ifs) for w in omega]
fit = [sbeta_fit(w * icm2ifs) for w in omega]
plot_qnsd(omega, [fw, fit], "figure/pswf_ohmic_sbeta.png"; labels=["Exact", "PSWF fit"])

# Correlation function from analytical solution (for comparison)
bcf_exact = BosonicBCF(sd, Temp; ub=6000.0)

# Set time range
t_max = 1500.0  # fs
n_points = 200
t_vals = range(0.01, t_max, length=n_points)

# Compute correlation functions
C_pswf = [pswfft(t) for t in t_vals]
C_exact = [bcf_exact(t) for t in t_vals]

# Compute errors
absolute_error = abs.(C_pswf - C_exact)
normalized_absolute_error = absolute_error ./ abs.(C_exact[1])

# Display statistics
max_abs_error = maximum(normalized_absolute_error)
mean_abs_error = sum(normalized_absolute_error) / length(normalized_absolute_error)

println("Normalized maximum absolute error: $(max_abs_error)")
println("Normalized mean absolute error:    $(mean_abs_error)")
plot_bcf(t_vals, C_pswf, C_exact, "figure/pswf_ohmic_bcf.png")