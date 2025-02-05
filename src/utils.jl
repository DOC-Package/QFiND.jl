# Description: Utility functions for the project

"""
    equispaced_mesh(N_t, N_ω, tc, Ω_min, Ω_max)

Generate equispaced time and frequency grids.
"""
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

# Calculate the sum of the exponential function.
function sumexp(t::Float64, a::AbstractVector{Float64}, c::AbstractVector{Float64}) :: ComplexF64
    return sum(c .* exp.(-im .* a .* t))
end

# Calculate the sum of the exponential function.
function sumexp(t::Float64, a::AbstractVector{ComplexF64}, c::AbstractVector{ComplexF64}) :: ComplexF64
    return sum(c .* exp.(-a .* t))
end