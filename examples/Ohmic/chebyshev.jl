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
n_terms = 30
cheb = chebyshev_expansion(sd, Temp, -300.0, 300.0, n_terms)
bcf_cheb = chebyshev_bcf(cheb)

t_max = 1000.0  # fs
n_points = 1000
t_vals = range(0.01, t_max, length=n_points)

# 相関関数を計算
C_cheb = [bcf_cheb(t) for t in t_vals]
C_exact = [bcf_exact(t) for t in t_vals]

# 誤差を計算
absolute_error = abs.(C_cheb - C_exact)
relative_error = absolute_error ./ abs.(C_exact[1])

# 統計情報を表示
max_abs_error = maximum(absolute_error)
max_rel_error = maximum(relative_error)
mean_abs_error = sum(absolute_error) / length(absolute_error)
mean_rel_error = sum(relative_error) / length(relative_error)

println("Chebyshev expansion analysis (n_terms = $n_terms)")
println("=" ^ 50)
println("Maximum absolute error: $(max_abs_error)")
println("Maximum relative error: $(max_rel_error) ($(max_rel_error * 100)%)")
println("Mean absolute error:    $(mean_abs_error)")
println("Mean relative error:    $(mean_rel_error) ($(mean_rel_error * 100)%)")

plot_bcf(t_vals, C_cheb, C_exact, "figure/chebyshev_bcf_ohmic.png")
