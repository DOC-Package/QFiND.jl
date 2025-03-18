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
    #
    save("qnsd.png", fig)
end

function plot_bcf(
    wk::AbstractVector{Float64}, 
    gk::AbstractVector{<:Number}, 
    bcf::Function, 
    Tc::Real, 
    N_t::Int,
    filename::String)

    t = range(0, Tc, length=N_t)
    reference = bcf.(t)
    approx = sumexp.(t, Ref(wk*icm2ifs), Ref(gk*icm2ifs^2.0))
    error = (approx - reference) ./ abs(bcf(0.0))

    ls1    = 20      # label font size
    ls2    = 15      # tick label font size
    lw1    = 2.5     # line width for main lines
    lw2    = 2.5     # line width for reference lines
    color1 = :orangered
    color2 = :royalblue
    color3 = :black

    fig = Figure(size = (800, 700))

    ax1 = Axis(fig[1, 1],
        xlabel = "",
        ylabel = L"C(t)",
        xlabelsize = ls2,
        ylabelsize = ls2
    )
    ax2 = Axis(fig[2, 1],
        xlabel = L"t (\mathrm{fs})",
        ylabel = L"\delta C(t)",
        xlabelsize = ls2,
        ylabelsize = ls2
    )
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
    lines!(ax2, t, real.(error),
        color = color1,
        linewidth = lw1
    )
    lines!(ax2, t, imag.(error),
        color = color2,
        linewidth = lw1
    )
    axislegend(ax1, position = :rt, labelsize = ls2)
    save(filename, fig)
    return fig
end

function plot_bcf(
    U::AbstractMatrix{ComplexF64}, 
    zk::AbstractVector{ComplexF64}, 
    reference::AbstractVector{ComplexF64},
    t::AbstractVector{<:Real},
    filename::String)

    cmax = maximum(abs.(reference))
    approx = U * zk
    error = (approx - reference) ./ cmax

    ls1    = 20      # label font size
    ls2    = 15      # tick label font size
    lw1    = 2.5     # line width for main lines
    lw2    = 2.5     # line width for reference lines
    color1 = :orangered
    color2 = :royalblue
    color3 = :black

    fig = Figure(size = (800, 700))

    ax1 = Axis(fig[1, 1],
        xlabel = "",
        ylabel = L"C(t)",
        xlabelsize = ls2,
        ylabelsize = ls2
    )
    ax2 = Axis(fig[2, 1],
        xlabel = L"t (\mathrm{fs})",
        ylabel = L"\delta C(t)",
        xlabelsize = ls2,
        ylabelsize = ls2
    )
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
    lines!(ax2, t, real.(error),
        color = color1,
        linewidth = lw1
    )
    lines!(ax2, t, imag.(error),
        color = color2,
        linewidth = lw1
    )
    axislegend(ax1, position = :rt, labelsize = ls2)
    save(filename, fig)
    return fig
end

function plot_basis_time(U::AbstractMatrix{ComplexF64}, t::AbstractVector{<:Real})
    ls1    = 20      # label font size
    ls2    = 15      # tick label font size
    lw1    = 2.5     # line width for main lines
    lw2    = 2.5     # line width for reference lines

    fig = Figure(size = (900, 700))

    ax1 = Axis(fig[1, 1],
    xlabel = L"t (\mathrm{fs})",
    ylabel = L"\mathrm{Re}\;U(t)",
    xlabelsize = ls2,
    ylabelsize = ls2
    )
    ax2 = Axis(fig[2, 1],
    xlabel = L"ω (\mathrm{cm^{-1}})",
    ylabel = L"\mathrm{Im}\;U(t)",
    xlabelsize = ls2,
    ylabelsize = ls2
    )

    num_lines = size(U, 2)
    colors_top    = [cgrad(:OrRd)[i] for i in LinRange(0, 1, num_lines)]
    colors_bottom = [cgrad(:Blues)[i] for i in LinRange(0, 1, num_lines)]

    for i in 1:num_lines
    lines!(ax1, t, real.(U[:, i]),
    label = string(i),
    color = colors_top[i],
    linewidth = lw1
    )
    end
    for i in 1:num_lines
    lines!(ax2, t, imag.(U[:, i]),
    label = string(i),
    color = colors_bottom[i],
    linewidth = lw1
    )
    end
    if num_lines <= 12
        leg1 = Legend(fig, ax1, labelsize = ls2)
        leg2 = Legend(fig, ax2, labelsize = ls2)
        fig[1, 2] = leg1
        fig[2, 2] = leg2
    end
    save("basis_time.png", fig)
    return fig
end

function plot_basis_freq(V::AbstractMatrix{ComplexF64}, ω::AbstractVector{<:Real})
    ls1    = 20      # label font size
    ls2    = 15      # tick label font size
    lw1    = 2.5     # line width for main lines
    lw2    = 2.5     # line width for reference lines

    fig = Figure(size = (900, 700))

    ax1 = Axis(fig[1, 1],
    xlabel = L"ω (\mathrm{cm}^{-1})",
    ylabel = L"\mathrm{Re}\;V(ω)",
    xlabelsize = ls2,
    ylabelsize = ls2
    )
    ax2 = Axis(fig[2, 1],
    xlabel = L"ω (\mathrm{cm^{-1}})",
    ylabel = L"\mathrm{Im}\;V(ω)",
    xlabelsize = ls2,
    ylabelsize = ls2
    )

    num_lines = size(V, 2)
    colors_top    = [cgrad(:OrRd)[i] for i in LinRange(0, 1, num_lines)]
    colors_bottom = [cgrad(:Blues)[i] for i in LinRange(0, 1, num_lines)]

    for i in 1:num_lines
    lines!(ax1, ω, real.(V[:, i]),
    label = string(i),
    color = colors_top[i],
    linewidth = lw1
    )
    end
    for i in 1:num_lines
    lines!(ax2, ω, imag.(V[:, i]),
    label = string(i),
    color = colors_bottom[i],
    linewidth = lw1
    )
    end

    if num_lines <= 12
        leg1 = Legend(fig, ax1, labelsize = ls2)
        leg2 = Legend(fig, ax2, labelsize = ls2)
        fig[1, 2] = leg1
        fig[2, 2] = leg2
    end
    save("basis_freq.png", fig)
    return fig
end