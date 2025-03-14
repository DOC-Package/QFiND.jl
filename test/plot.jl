using QFiND
using CairoMakie
using LaTeXStrings 

function plot_qnsd(sbeta::Function, Ω_min::Real, Ω_max::Real, N_freq::Int)
    w = range(Ω_min, Ω_max, length=N_freq) |> collect
    j = sbeta.(w)
    # Plot the spectral density and save the plot
    fig = Figure()
    ax = Axis(fig[1, 1],
    xlabel = L"\text{Frequency} (\mathrm{cm^{-1}})",
    ylabel = L"\text{Spectral Density}",
    xlabelsize=20,
    ylabelsize=20,
    xgridvisible = false,
    ygridvisible = false,
    )
    lines!(ax, w, j, color = :red, linewidth = 2)
    xlims!(ax, (Ω_min, Ω_max))
    ylims!(ax, (0, nothing))
    #axislegend(ax)
    save("qnsd.png", fig)
end

function plot_bcf(
    wk::AbstractVector{Float64}, 
    gk::AbstractVector{Float64}, 
    bcf::Function, 
    Tc::Real, 
    N_t::Int)

    t = range(0, Tc, length=N_t)
    reference = bcf.(t)
    approx = sumexp.(t, Ref(wk*icm2ifs), Ref(gk*icm2ifs^2.0))
    error = (approx - reference) ./ abs(bcf(0.0))

    # Style parameters (matching Python style)
    ls1    = 20      # label font size
    ls2    = 15      # tick label font size
    lw1    = 2.5     # line width for main lines
    lw2    = 2.5     # line width for reference lines
    color1 = :orangered
    color2 = :royalblue
    color3 = :black

    # Create a figure with 2 rows and 1 column.
    fig = Figure(size = (800, 700))

    # Top axis: Bath correlation function (BCF) plot.
    ax1 = Axis(fig[1, 1],
        xlabel = "",
        ylabel = L"C(t)",
        xlabelsize = ls2,
        ylabelsize = ls2
    )

    # Bottom axis: Error plot.
    ax2 = Axis(fig[2, 1],
        xlabel = L"t (\mathrm{fs})",
        ylabel = L"\delta C(t)",
        xlabelsize = ls2,
        ylabelsize = ls2
    )

    # Plot on the top axis: Approximate BCF (real and imaginary parts)
    lines!(ax1, t, real.(approx),
        label = L"\mathrm{Re}\,C(t)",
        color = color1,
        linewidth = lw1
    )
    lines!(ax1, t, imag.(approx),
        label = L"\mathrm{Im}\,C(t)",
        color = color2,
        linewidth = lw1
    )

    # Plot on the top axis: Reference (reference) BCF as dashed lines.
    lines!(ax1, t, real.(reference),
        label = L"\text{Reference}",
        color = color3,
        linestyle = :dash,
        linewidth = lw2
    )
    lines!(ax1, t, imag.(reference),
        color = color3,
        linestyle = :dash,
        linewidth = lw2
    )

    # Plot on the bottom axis: Error (real and imaginary parts)
    lines!(ax2, t, real.(error),
        color = color1,
        linewidth = lw1
    )
    lines!(ax2, t, imag.(error),
        color = color2,
        linewidth = lw1
    )

    # Add a legend to the top axis at the top-right position.
    axislegend(ax1, position = :rt, labelsize = ls2)

    # Save the figure as a PNG file.
    save("bcf.png", fig)

    return fig
end
