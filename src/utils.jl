# Description: Utility functions for the project

function equispaced_grid(
    Ω_min::Real, 
    Ω_max::Real, 
    tc::Real; 
    n_freq::Integer=1000, 
    n_time::Integer=200, 
    scale::Real=1.0)

    t = collect(range(0, tc, length=n_time))
    ω = collect(range(Ω_min, Ω_max, length=n_freq))
    ω = ω .* scale
    return t, ω
end

function equispaced_grid(T_c::Real; n_time::Integer=200)
    t = collect(range(0, T_c, length=n_time))
    return t
end

function equispaced_grid(Ω_min::Real, Ω_max::Real; n_freq::Integer=1000, scale::Real=1.0)
    ω = collect(range(Ω_min, Ω_max, length=n_freq))
    ω = ω .* scale
    return ω
end

function create_integrand(sbeta::Vector{<:Real}, t::Vector{<:Real}, ω::Vector{<:Real}) :: Matrix{Float64}
    N_t = length(t)
    N_ω = length(ω)
    f = zeros(Float64, 2*N_t, N_ω)
    # Fill first N_t rows (real part)
    for i in 1:N_t
        for j in 1:N_ω
            f[i, j] = sbeta[j] * cos(ω[j] * t[i])
        end
    end
    # Fill next N_t rows (imaginary part)
    for i in (N_t+1):(2*N_t)
        for j in 1:N_ω
            f[i, j] = -sbeta[j] * sin(ω[j] * t[i - N_t])
        end
    end
    return f
end

function create_integrand_c(sbeta::Vector{<:Real}, t::Vector{<:Real}, ω::Vector{<:Real}) :: Matrix{ComplexF64}
    N_t = length(t)
    N_ω = length(ω)
    f = zeros(ComplexF64, N_t, N_ω)
    for i in 1:N_t
        for j in 1:N_ω
            f[i, j] = sbeta[j] * exp(-1.0im * ω[j] * t[i])
        end
    end

    return f
end

function linsys_weight(cc::AbstractVector{<:Complex{<:Real}}, idx::Vector{Int}, D::Matrix{ComplexF64})
    frank = size(D, 2)
    c = cc[idx[1:frank]]
    # Solve the linear system: Kr * g = c
    g = D \ c
    err = norm(D * g - c)

    return g, err
end


# Calculate the sum of the exponential function.
function sumexp(t::Float64, a::AbstractVector{Float64}, c::AbstractVector{Float64}) :: ComplexF64
    return sum(c .* exp.(-im .* a .* t))
end

function sumexp(t::Float64, a::AbstractVector{ComplexF64}, c::AbstractVector{ComplexF64}) :: ComplexF64
    return sum(c .* exp.(-a .* t))
end

function sumexp(t::Float64, a::AbstractVector{Float64}, c::AbstractVector{ComplexF64}) :: ComplexF64
    return sum(c .* exp.(-im .* a .* t))
end

function ir_time(t::Float64, U::AbstractVector{ComplexF64}, z::AbstractVector{ComplexF64}) :: ComplexF64
    return sum(c .* exp.(-im .* a .* t))
end
