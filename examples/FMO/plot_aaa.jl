include("../plot.jl")
using RationalFunctionApproximation
using Serialization
using DelimitedFiles
using CairoMakie
using LaTeXStrings

function plot_aaa_comparison(
    w::AbstractVector{<:Real}, 
    j_approx::AbstractVector{<:Real}, 
    j_ref::AbstractVector{<:Real}, 
    filename::String;
    labels::Tuple{String,String}=("AAA", "Reference"))

    error = j_approx .- j_ref

    ls1 = 20      # label font size
    ls2 = 15      # tick label font size
    lw1 = 2.5     # line width for approximation
    lw2 = 2.0     # line width for reference (dashed)
    color1 = :red
    color2 = :black
    color3 = :red

    fig = Figure(size = (800, 700))

    # Upper panel: spectral density comparison
    ax1 = Axis(fig[1, 1],
        xlabel = L"\omega \, (\mathrm{cm^{-1}})",
        ylabel = L"J(\omega)",
        xlabelsize = ls1,
        ylabelsize = ls1,
        xgridvisible = false,
        ygridvisible = false,
    )

    
    # Approximation as solid line
    lines!(ax1, w, j_approx, color = color1, linewidth = lw1, label = labels[1])
    # Reference as dashed line
    lines!(ax1, w, j_ref, color = color2, linewidth = lw2, linestyle = :dash, label = labels[2])

    xlims!(ax1, (w[1], w[end]))
    ylims!(ax1, (0, nothing))
    ax1.xticklabelsize = ls2
    ax1.yticklabelsize = ls2
    axislegend(ax1, position = :rt, labelsize = 14)

    # Lower panel: error
    ax2 = Axis(fig[2, 1],
        xlabel = L"\omega \, (\mathrm{cm^{-1}})",
        ylabel = L"J_\mathrm{approx} - J_\mathrm{ref}",
        xlabelsize = ls1,
        ylabelsize = ls1,
        xgridvisible = false,
        ygridvisible = false,
    )

    lines!(ax2, w, error, color = color3, linewidth = lw1)
    hlines!(ax2, [0.0], color = :gray, linewidth = 1, linestyle = :dash)

    xlims!(ax2, (w[1], w[end]))
    ax2.xticklabelsize = ls2
    ax2.yticklabelsize = ls2

    # Link x-axes
    linkxaxes!(ax1, ax2)

    save(filename, fig)
    return fig
end

data = readdlm("fmo_smoothed.txt")
col1 = data[:, 1]
col2 = data[:, 2] 
ω = Float64.(col1) 
J = Float64.(col2)
J[J .< 0.0] .= 0.0

r = open("jrat_fmo.bin", "r") do io
    deserialize(io)
end

plot_aaa_comparison(ω, r.(ω), J, "figure/aaa_fmo.png"; labels=("AAA", "Reference"))