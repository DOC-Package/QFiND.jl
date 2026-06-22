include("../plot.jl")
using QFiND
using LinearAlgebra
using RationalFunctionApproximation
  
# spectral density
s = 1.0
alpha = 50.0
gamc = 50.0
Temp = 300.0
Ω_min = -600.0
Ω_max = 600.0
N_ω = 2000
eps = 1e-2

sdens = PowerLawExpSD(s, gamc; alpha=alpha)
sbeta = BosonicQNSD(sdens, Temp)

ω = collect(range(Ω_min, Ω_max, length=N_ω)) 
Sbeta = [sbeta(ωi) for ωi in ω]

# AAA approximation of S_beta
r = aaa(ω, Sbeta; tol=eps, max_degree=100, lookahead=50)
println("degree: ", length(r.nodes))
if length(r.nodes) > 2
    r = lawson(ω, Sbeta, r, 20)
end
err = norm(r.(ω) - Sbeta)
println("error: ", err)

pols, res, rinf = AAAPartialFraction(r)
Srat = RationalSD(pols, res, rinf)
err = norm(Srat.(ω) - Sbeta)
println("error (RationalSD): ", err)

Srat_no_rinf = RationalSD(pols, res)
err = norm(Srat_no_rinf.(ω) - Sbeta)
println("error (RationalSD wo rinf): ", err)

# poles_to_exponents for correlation function
result = poles_to_exponents(pols, res)
c_k = result.c .* icm2ifs^2.0*2.0
a_k = result.a .* icm2ifs
println("\nNumber of exponential terms: ", length(c_k))
println("Poles in lower half plane: ", pols[findall(z -> imag(z) < 0, complex.(pols))])
println("r_inf: ", rinf)

# Exponential sum representation of C(t): C(t) = Σ c_k exp(-a_k t)
function C_exp_sum(t::Real, c::Vector{ComplexF64}, a::Vector{ComplexF64})
    if isempty(c)
        return 0.0 + 0.0im
    end
    return sum(c .* exp.(-a .* t))
end

# Compare with BosonicBCF
bcf = BosonicBCF(sdens, Temp; ub=Ω_max)
t_max = 1000.0  # fs
N_t = 100
t = collect(range(0.0, t_max, length=N_t))

# Compute both correlation functions
C_exact = [bcf(ti) for ti in t]
# a_k is in cm⁻¹, multiply by icm2ifs to get fs⁻¹
C_approx = [C_exp_sum(ti, c_k, a_k) for ti in t]

# Numerical Fourier transform of AAA approximation r(ω)
# C_FT(t) = 1/(2π) ∫ r(ω) e^{-iωt} dω
# ω is in cm⁻¹, t is in fs, so we need to convert units properly
using QuadGK
function C_fourier_transform(ti::Real, r_func, ω_min::Real, ω_max::Real)
    # Convert ω from cm⁻¹ to fs⁻¹ inside the integral
    # dω also needs to be converted: dω[cm⁻¹] * icm2ifs = dω[fs⁻¹]
    integrand = ω -> r_func(ω) * exp(-im * ω * icm2ifs * ti)
    result, _ = quadgk(integrand, ω_min, ω_max; rtol=1e-6)
    return result * icm2ifs^2.0 / π  # icm2ifs for dω conversion
end
println("Computing numerical Fourier transform of AAA approximation...")
C_ft = [C_fourier_transform(ti, r, Ω_min, Ω_max) for ti in t]
println("Done.")

# Numerical FT of RationalSD (sum of poles) with rinf
println("Computing numerical Fourier transform of RationalSD (with r_inf)...")
C_ft_srat = [C_fourier_transform(ti, Srat, Ω_min, Ω_max) for ti in t]
println("Done.")

# Numerical FT of RationalSD (sum of poles) without rinf
println("Computing numerical Fourier transform of RationalSD (without r_inf)...")
C_ft_srat_no_rinf = [C_fourier_transform(ti, Srat_no_rinf, Ω_min, Ω_max) for ti in t]
println("Done.")

# ========== Visualization ==========
fig = Figure(size = (1200, 1000))

# Plot 1: S_beta(ω) comparison
ax1 = Axis(fig[1, 1],
    xlabel = L"\omega \, (\mathrm{cm}^{-1})",
    ylabel = L"S_\beta(\omega)",
    title = "Quantum Noise Spectral Density"
)
lines!(ax1, ω, Sbeta, label = "Original", linewidth = 2)
lines!(ax1, ω, real.(r.(ω)), label = "AAA", linewidth = 2, linestyle = :dash)
axislegend(ax1, position = :rt)

# Plot 2: Poles in complex plane
ax2 = Axis(fig[1, 2],
    xlabel = L"\mathrm{Re}(z_k)",
    ylabel = L"\mathrm{Im}(z_k)",
    title = "Poles of S_β(ω)"
)
scatter!(ax2, real.(pols), imag.(pols), markersize = 10)
hlines!(ax2, [0.0], color = :gray, linestyle = :dash)

# Plot 3: Correlation function comparison (Real part)
ax3 = Axis(fig[2, 1],
    xlabel = L"t \, (\mathrm{fs})",
    ylabel = L"\mathrm{Re}[C(t)]",
    title = "Correlation Function (Real)"
)
lines!(ax3, t, real.(C_exact), label = "BosonicBCF", linewidth = 2)
lines!(ax3, t, real.(C_ft), label = "FT of AAA", linewidth = 2, linestyle = :dot)
lines!(ax3, t, real.(C_ft_srat), label = "FT of Srat (w/ r∞)", linewidth = 2, linestyle = :dashdot)
lines!(ax3, t, real.(C_ft_srat_no_rinf), label = "FT of Srat (w/o r∞)", linewidth = 2, linestyle = :dashdotdot)
if !isempty(c_k)
    lines!(ax3, t, real.(C_approx), label = "Exp. Sum", linewidth = 2, linestyle = :dash)
end
axislegend(ax3, position = :rt)

# Plot 4: Correlation function comparison (Imaginary part)
ax4 = Axis(fig[2, 2],
    xlabel = L"t \, (\mathrm{fs})",
    ylabel = L"\mathrm{Im}[C(t)]",
    title = "Correlation Function (Imaginary)"
)
lines!(ax4, t, imag.(C_exact), label = "BosonicBCF", linewidth = 2)
lines!(ax4, t, imag.(C_ft), label = "FT of AAA", linewidth = 2, linestyle = :dot)
lines!(ax4, t, imag.(C_ft_srat), label = "FT of Srat (w/ r∞)", linewidth = 2, linestyle = :dashdot)
lines!(ax4, t, imag.(C_ft_srat_no_rinf), label = "FT of Srat (w/o r∞)", linewidth = 2, linestyle = :dashdotdot)
if !isempty(c_k)
    lines!(ax4, t, imag.(C_approx), label = "Exp. Sum", linewidth = 2, linestyle = :dash)
end
axislegend(ax4, position = :rt)

# Plot 5: Error in correlation function
ax5 = Axis(fig[3, 1:2],
    xlabel = L"t \, (\mathrm{fs})",
    ylabel = L"|C_\mathrm{exact} - C_\mathrm{approx}|",
    title = "Correlation Function Error",
    yscale = log10
)
lines!(ax5, t[2:end], abs.(C_exact[2:end] .- C_ft[2:end]), label = "FT of AAA", linewidth = 2)
lines!(ax5, t[2:end], abs.(C_exact[2:end] .- C_ft_srat[2:end]), label = "FT of Srat (w/ r∞)", linewidth = 2)
lines!(ax5, t[2:end], abs.(C_exact[2:end] .- C_ft_srat_no_rinf[2:end]), label = "FT of Srat (w/o r∞)", linewidth = 2)
if !isempty(c_k)
    lines!(ax5, t[2:end], abs.(C_exact[2:end] .- C_approx[2:end]), label = "Exp. Sum", linewidth = 2)
end
axislegend(ax5, position = :rt)

save("figure/aaa_visualization.png", fig)
println("Figure saved to figure/aaa_visualization.png")