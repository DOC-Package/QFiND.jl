using QFiND
include("../plot.jl")

# Ohmic spectral density
s = 1.0
alpha = 50.0
gamc = 50.0
sd = PowerLawExpSD(s, gamc; alpha=alpha)
Temp = 300.0
bcf_exact = BosonicBCF(sd, Temp)

# Chebyshev
n_terms = 50
cheb = chebyshev_expansion(sd, Temp, -300.0, 300.0, n_terms)
bcf_cheb = chebyshev_bcf(cheb)

t_max = 1000.0  # fs
n_points = 1000
t_vals = range(0.01, t_max, length=n_points)

C_cheb = [bcf_cheb(t) for t in t_vals]
C_exact = [bcf_exact(t) for t in t_vals]

# error
evaluate_error(t_vals, C_cheb, C_exact)
plot_bcf(t_vals, C_cheb, C_exact, "figure/bcf_chebyshev.png")
