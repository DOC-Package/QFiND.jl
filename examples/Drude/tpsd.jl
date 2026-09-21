include("../plot.jl")
using QFiND
using ExpFit

γ = 50.0
λ = 30.0
eps = 1e-6
npade = 11
ntpsd = 1
temps = [50.0, 300.0]
figure_dir = joinpath(@__DIR__, "figure")

sdens = DrudeSD(γ, λ)

T_max = 500.0
N_t = 2000
t = collect(range(0.0, T_max, length=N_t))

function compare_temperatures_plot(
	t::AbstractVector{<:Real},
	references::Dict{Float64, Vector{ComplexF64}},
	approximations::Dict{Float64, Vector{ComplexF64}},
	filename::String,
)
	fig = Figure(size = (1000, 700))
	ax_re = Axis(fig[1, 1],
		xlabel = "",
		ylabel = "Re C(t)",
		xlabelsize = 20,
		ylabelsize = 20,
	)
	ax_im = Axis(fig[2, 1],
		xlabel = L"t (\mathrm{fs})",
		ylabel = "Im C(t)",
		xlabelsize = 20,
		ylabelsize = 20,
	)

	colors = Dict(50.0 => :royalblue, 300.0 => :orangered)
	for temp in sort(collect(keys(references)))
		color = colors[temp]
		label_prefix = "$(Int(round(temp))) K"
		lines!(ax_re, t, real.(references[temp]), color = color, linestyle = :dash, linewidth = 2.5,
			label = label_prefix * " PSD")
		lines!(ax_re, t, real.(approximations[temp]), color = color, linewidth = 2.5,
			label = label_prefix * " TPSD")
		lines!(ax_im, t, imag.(references[temp]), color = color, linestyle = :dash, linewidth = 2.5,
			label = label_prefix * " PSD")
		lines!(ax_im, t, imag.(approximations[temp]), color = color, linewidth = 2.5,
			label = label_prefix * " TPSD")
	end

	axislegend(ax_re, position = :rt, labelsize = 16)
	xlims!(ax_re, (first(t), last(t)))
	xlims!(ax_im, (first(t), last(t)))
	save(filename, fig)
	return fig
end

function compare_temperatures_psd_real_plot(
	t::AbstractVector{<:Real},
	references::Dict{Float64, Vector{ComplexF64}},
	filename::String,
)
	fig = Figure(size = (850, 450))
	ax = Axis(fig[1, 1],
		xlabel = L"t (\mathrm{fs})",
		ylabel = "Re C(t)",
		xlabelsize = 20,
		ylabelsize = 20,
	)

	colors = Dict(50.0 => :royalblue, 300.0 => :orangered)
	for temp in sort(collect(keys(references)))
		color = colors[temp]
		label = "$(Int(round(temp))) K"
		lines!(ax, t, real.(references[temp]), color = color, linewidth = 3, label = label)
	end

	axislegend(ax, position = :rt, labelsize = 16)
	xlims!(ax, (first(t), last(t)))
	save(filename, fig)
	return fig
end

references = Dict{Float64, Vector{ComplexF64}}()
approximations = Dict{Float64, Vector{ComplexF64}}()

for temp in temps
	a, c = psd(sdens, temp, npade; pade_type=:N)
	at, ct = tpsd(sdens, temp, npade, ntpsd; pade_type=:N)
	bcf = bcf_approx(t, a, c)
	bcf_t = bcf_approx(t, at, ct)

	references[temp] = bcf
	approximations[temp] = bcf_t

	println("temperature: ", temp, " K")
	println("degree: ", size(at))
	evaluate_error(t, bcf_t, bcf)
	plot_bcf(t, bcf_t, bcf, joinpath(figure_dir, "bcf_psd_$(Int(round(temp)))K.png"))
	save_expon_coeff(a ./ icm2ifs, c ./ icm2ifs^2.0,
		joinpath(@__DIR__, "expon_coeff_d50_l30_$(Int(round(temp)))K_psd11.txt"))
	save_expon_coeff(at ./ icm2ifs, ct ./ icm2ifs^2.0,
		joinpath(@__DIR__, "expon_coeff_d50_l30_$(Int(round(temp)))K_tpsd1.txt"))
end

compare_temperatures_plot(t, references, approximations, joinpath(figure_dir, "bcf_temp_compare.png"))
compare_temperatures_psd_real_plot(t, references, joinpath(figure_dir, "bcf_temp_compare_psd_real.png"))