using QFiND
using CairoMakie
using LaTeXStrings 

function plot_qnsd(
    ω::AbstractVector{<:Real},
    approx::AbstractVector{<:Real},
    reference::AbstractVector{<:Real},
    filename::String)

    error = (approx - reference) ./ maximum(reference)

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
        ylabel = L"J(\omega)",
        xlabelsize = ls2,
        ylabelsize = ls2
    )
    ax2 = Axis(fig[2, 1],
        xlabel = L"\omega (\mathrm{cm}^{-1})",
        ylabel = L"\delta J(\omega)",
        xlabelsize = ls2,
        ylabelsize = ls2
    )
    lines!(ax1, ω, approx,
        label = L"\text{Approximation}",
        color = color1,
        linewidth = lw1
    )
    lines!(ax1, ω, reference,
        label = L"\text{Reference}",
        color = color3,
        linestyle = :dash,
        linewidth = lw2
    )
    lines!(ax2, ω, error,
        color = color1,
        linewidth = lw1
    )
    axislegend(ax1, position = :rt, labelsize = ls2)
    save(filename, fig)
    return fig
end

function plot_qnsd(w::AbstractVector{<:Real}, j::AbstractVector{<:Real}, filename::String)
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
    xlims!(ax, (w[1], w[end]))
    ylims!(ax, (0, nothing))
    #
    save(filename, fig)
end

function plot_qnsd(sbeta::Function, Ω_min::Real, Ω_max::Real, N_freq::Int, filename::String)
    w = range(Ω_min, Ω_max, length=N_freq) |> collect
    j = sbeta.(w)
    return plot_qnsd(w, j, filename)
end

function plot_bcf(
    t::AbstractVector{<:Real},
    approx::AbstractVector{ComplexF64},
    reference::AbstractVector{ComplexF64},
    filename::String)

    error = (approx - reference) ./ abs(reference[1])

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
    ω::AbstractVector{Float64}, 
    g::AbstractVector{<:Number}, 
    bcf::Function, 
    Tc::Real, 
    N_t::Int,
    filename::String)

    t = range(0, Tc, length=N_t)
    plot_bcf(ω, g, bcf, t, filename)
end

function plot_bcf(
    ω::AbstractVector{Float64}, 
    g::AbstractVector{<:Number}, 
    reference::AbstractVector{<:Number},   
    t::AbstractVector{<:Real},
    filename::String)

    approx = bcf_approx.(t, Ref(ω), Ref(g))
    return plot_bcf(t, approx, reference, filename)
end

function plot_bcf(
    ω::AbstractVector{Float64}, 
    g::AbstractVector{<:Number}, 
    bcf::Function, 
    t::AbstractVector{<:Real},
    filename::String)

    approx = bcf_approx.(t, Ref(ω), Ref(g))
    return plot_bcf(t, approx, bcf.(t), filename)
end

function plot_bcf(
    a::AbstractVector{ComplexF64}, 
    c::AbstractVector{ComplexF64}, 
    bcf::Function, 
    Tc::Real, 
    N_t::Int,
    filename::String)

    t = range(0, Tc, length=N_t)
    reference = bcf.(t)
    approx = bcf_approx.(t, a, c)
    return plot_bcf(t, approx, reference, filename)
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
    return plot_bcf(t, approx, reference, error, filename)
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