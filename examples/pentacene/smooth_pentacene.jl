using DelimitedFiles
using LinearAlgebra
using Printf
using Statistics

function tricube_weight(u::Float64)
    t = abs(u)
    return t >= 1.0 ? 0.0 : (1.0 - t^3)^3
end

function local_poly_value(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    idx::UnitRange{Int},
    x0::Float64;
    degree::Int=2,
)
    xwin = Float64.(view(x, idx)) .- x0
    ywin = Float64.(view(y, idx))
    radius = maximum(abs, xwin)

    if radius == 0.0
        return ywin[1]
    end

    weights = tricube_weight.(xwin ./ radius)
    if all(iszero, weights)
        return y[Int(clamp(round(x0 - first(x) + 1), 1, length(y)))]
    end

    cols = ntuple(p -> xwin .^ (p - 1), degree + 1)
    design = hcat(cols...)
    sqrtw = sqrt.(weights)
    weighted_design = design .* sqrtw
    weighted_y = ywin .* sqrtw

    coeffs = try
        weighted_design \ weighted_y
    catch
        return dot(weights, ywin) / sum(weights)
    end

    return coeffs[1]
end

function smoothing_window(n::Int, span::Float64, degree::Int)
    window = max(degree + 3, round(Int, span * n))
    window = min(window, n)
    return isodd(window) ? window : max(degree + 3, window + 1 > n ? window - 1 : window + 1)
end

function centered_window(i::Int, n::Int, window::Int)
    half = window ÷ 2
    lo = max(1, i - half)
    hi = min(n, i + half)

    if hi - lo + 1 < window
        if lo == 1
            hi = min(n, window)
        else
            lo = max(1, n - window + 1)
        end
    end

    return lo:hi
end

function loess_smooth(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real};
    span::Float64=0.01,
    degree::Int=2,
    passes::Int=2,
)
    n = length(x)
    @assert n == length(y)
    @assert 0.0 < span <= 1.0
    @assert degree >= 0
    @assert passes >= 1

    window = smoothing_window(n, span, degree)
    smoothed = Float64.(y)

    for _ in 1:passes
        current = similar(smoothed)
        for i in eachindex(smoothed)
            idx = centered_window(i, n, window)
            current[i] = local_poly_value(x, smoothed, idx, Float64(x[i]); degree=degree)
        end
        smoothed .= current
    end

    return smoothed
end

function mean_conversion_ratio(col2::AbstractVector{<:Real}, col3::AbstractVector{<:Real})
    mask = (abs.(col2) .> 0.0) .& (abs.(col3) .> 0.0)
    return mean(Float64.(col2[mask] ./ col3[mask]))
end

function smooth_file(
    input_path::AbstractString,
    output_path::AbstractString;
    span::Float64=0.01,
    degree::Int=2,
    passes::Int=2,
)
    data = readdlm(input_path)
    size(data, 2) == 3 || throw(ArgumentError("expected a 3-column data file: $input_path"))

    omega = Float64.(data[:, 1])
    column2 = Float64.(data[:, 2])
    column3 = Float64.(data[:, 3])

    smoothed2 = loess_smooth(omega, column2; span=span, degree=degree, passes=passes)
    smoothed2 .= max.(smoothed2, 0.0)

    ratio = mean_conversion_ratio(column2, column3)
    smoothed3 = smoothed2 ./ ratio

    output = hcat(omega, smoothed2, smoothed3)
    open(output_path, "w") do io
        for row in eachrow(output)
            @printf(io, "%12.4f %16.6e %16.6e\n", row[1], row[2], row[3])
        end
    end

    rms = sqrt(sum((smoothed2 .- column2) .^ 2) / length(column2))
    println("wrote $(output_path)")
    println("span=$(span), degree=$(degree), passes=$(passes)")
    println("RMS difference (column 2): $(rms)")
end

input_path = length(ARGS) >= 1 ? ARGS[1] : "pentacene_h.dat"
output_path = length(ARGS) >= 2 ? ARGS[2] : "pentacene_h_smoothed.dat"
span = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.01
degree = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 2
passes = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 2

smooth_file(input_path, output_path; span=span, degree=degree, passes=passes)