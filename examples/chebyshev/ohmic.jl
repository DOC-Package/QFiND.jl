using QFiND
include("../plot.jl")
#include("../matio.jl")

s = 1.0
alpha = 5.0
gamc = 50.0
sd = PowerLawExpSD(s, gamc; alpha=alpha)
reorgene = reorganization_energy(sd; ub=6000.0)
println("reorganization energy: ", reorgene)
Temp = 300.0
bcf_exact = BosonicBCF(sd, Temp; ub=6000.0)
bcf_dt = BosonicBCF_dt(sd, Temp; ub=6000.0)

# Chebyshev
n_terms = 40
cheb = chebyshev_expansion(sd, Temp, -400.0, 400.0, n_terms)
bcf_cheb = chebyshev_bcf(cheb)


t_max = 2000.0  # fs
n_points = 1000
t_vals = range(0.01, t_max, length=n_points)

C_cheb = [bcf_cheb(t) for t in t_vals]
C_exact = [bcf_exact(t) for t in t_vals]

evaluate_error(t_vals, C_cheb, C_exact)
plot_bcf(t_vals, C_cheb, C_exact, "figure/chebyshev_bcf_ohmic.png")

dC = [bcf_dt(t) for t in t_vals]
dC_cheb = [chebyshev_expansion_dt(t, cheb) for t in t_vals]

evaluate_error(t_vals, dC_cheb, dC)
plot_bcf(t_vals, dC_cheb, dC, "figure/chebyshev_bcf_dt_ohmic.png")

D = chebyshev_derivative_matrix(cheb)
#D = D ./ icm2ifs
#write_complex_matrix_txt("D_ohmic.txt", D; delimiter = ",")

#phi0 = [cheb.basis[i](0.0) for i in 1:length(cheb.coeffs)]
#save_coeff_phi0(cheb.coeffs ./ icm2ifs^2.0, phi0, "ohmic_coeffs_phi0.txt")


