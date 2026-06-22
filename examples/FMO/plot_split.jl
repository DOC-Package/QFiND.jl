include("../plot.jl")
using RationalFunctionApproximation
using Serialization
using DelimitedFiles
using CairoMakie
using LaTeXStrings

function plot_split_sd(
    w::AbstractVector{<:Real}, 
    j_total::AbstractVector{<:Real},
    j_narrow::AbstractVector{<:Real}, 
    j_broad::AbstractVector{<:Real}, 
    j_ref::AbstractVector{<:Real}, 
    filename::String;
    threshold::Real=NaN)

    ls1 = 20      # label font size
    ls2 = 15      # tick label font size
    lw1 = 2.5     # line width for main lines
    lw2 = 2.0     # line width for reference (dashed)

    fig = Figure(size = (800, 900))

    # Panel 1: Total comparison
    ax1 = Axis(fig[1, 1],
        ylabel = L"J(\omega)",
        xlabelsize = ls1,
        ylabelsize = ls1,
        xgridvisible = false,
        ygridvisible = false,
        title = "Total"
    )
    lines!(ax1, w, j_total, color = :red, linewidth = lw1, label = "AAA Total")
    lines!(ax1, w, j_ref, color = :black, linewidth = lw2, linestyle = :dash, label = "Reference")
    xlims!(ax1, (w[1], w[end]))
    ylims!(ax1, (0, nothing))
    ax1.xticklabelsize = ls2
    ax1.yticklabelsize = ls2
    axislegend(ax1, position = :rt, labelsize = 12)

    # Panel 2: Narrow linewidth contribution
    threshold_str = isnan(threshold) ? "" : " (|Im| < $(threshold))"
    ax2 = Axis(fig[2, 1],
        ylabel = L"J(\omega)",
        xlabelsize = ls1,
        ylabelsize = ls1,
        xgridvisible = false,
        ygridvisible = false,
        title = "Narrow linewidth" * threshold_str
    )
    lines!(ax2, w, j_narrow, color = :blue, linewidth = lw1, label = "Narrow")
    lines!(ax2, w, j_ref, color = :black, linewidth = lw2, linestyle = :dash, label = "Reference")
    xlims!(ax2, (w[1], w[end]))
    ax2.xticklabelsize = ls2
    ax2.yticklabelsize = ls2
    axislegend(ax2, position = :rt, labelsize = 12)

    # Panel 3: Broad linewidth contribution
    threshold_str = isnan(threshold) ? "" : " (|Im| ≥ $(threshold))"
    ax3 = Axis(fig[3, 1],
        xlabel = L"\omega \, (\mathrm{cm^{-1}})",
        ylabel = L"J(\omega)",
        xlabelsize = ls1,
        ylabelsize = ls1,
        xgridvisible = false,
        ygridvisible = false,
        title = "Broad linewidth" * threshold_str
    )
    lines!(ax3, w, j_broad, color = :green, linewidth = lw1, label = "Broad")
    xlims!(ax3, (w[1], w[end]))
    ax3.xticklabelsize = ls2
    ax3.yticklabelsize = ls2
    axislegend(ax3, position = :rt, labelsize = 12)

    # Link x-axes
    linkxaxes!(ax1, ax2, ax3)

    save(filename, fig)
    return fig
end

# Evaluate partial fraction expansion
function eval_partial_fraction(w::AbstractVector{<:Real}, poles::Vector{ComplexF64}, 
                                residues::Vector{ComplexF64}, r_inf::Real=0.0; clamp_negative::Bool=true)
    j = zeros(length(w))
    for (i, ω) in enumerate(w)
        val = r_inf
        for (p, r) in zip(poles, residues)
            val += real(r / (ω - p))
        end
        j[i] = clamp_negative ? max(val, 0.0) : val
    end
    return j
end

# Load data
data = readdlm("fmo_smoothed.txt")
ω = Float64.(data[:, 1])
J = Float64.(data[:, 2])
J[J .< 0.0] .= 0.0

# Load AAA approximation (Barycentric form)
r = open("r_fmo.bin", "r") do io
    deserialize(io)
end

# Get poles, residues, and r_inf
pols, res = QFiND.AAAresidues(r)
r_inf = QFiND.AAArinf(r)

# Split by linewidth threshold
threshold = 10.0  # Adjust this value as needed
narrow, broad = QFiND.split_poles_by_linewidth(pols, res; threshold=threshold)

println("Total poles: ", length(pols))
println("Narrow linewidth poles (|Im| < $threshold): ", length(narrow.poles))
println("Broad linewidth poles (|Im| >= $threshold): ", length(broad.poles))

# Evaluate each contribution (without clamping for accuracy check)
j_total = r.(ω)
j_narrow_raw = eval_partial_fraction(ω, narrow.poles, narrow.residues, r_inf; clamp_negative=false)
j_broad_raw = eval_partial_fraction(ω, broad.poles, broad.residues, 0.0; clamp_negative=false)
j_combined_raw = j_narrow_raw .+ j_broad_raw 

# Clamped versions for plotting
j_narrow = eval_partial_fraction(ω, narrow.poles, narrow.residues, r_inf; clamp_negative=true)
j_broad = eval_partial_fraction(ω, broad.poles, broad.residues, 0.0; clamp_negative=true)

# Compute errors
using LinearAlgebra
err_total = norm(j_total .- J)
err_narrow = norm(j_narrow .- J)
err_broad = norm(j_broad .- J)
err_combined_raw = norm(j_combined_raw .- J)

println("\n--- Error Analysis ---")
println("Error (AAA total):         ", err_total)
println("Error (Narrow only):       ", err_narrow)
println("Error (Broad only):        ", err_broad)
println("Error (Narrow+Broad raw):  ", err_combined_raw)
println("Max |combined - total|:    ", maximum(abs.(j_combined_raw .- j_total)))

# Plot
plot_split_sd(ω, j_total, j_narrow, j_broad, J, "figure/split_fmo.png"; threshold=threshold)

println("\nPlot saved to figure/split_fmo.png")
