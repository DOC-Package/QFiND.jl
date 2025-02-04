# Description: Utility functions for the project

"""
    equispaced_mesh(N_t, N_ω, tc, Ω_min, Ω_max)

Generate equispaced time and frequency grids.
"""
function equispaced_grid(N_t::Integer, N_ω::Integer, tc::Real, Ω_min::Real, Ω_max::Real)
    t = collect(range(0, tc, length=N_t))
    ω = collect(range(Ω_min, Ω_max, length=N_ω))
    ω = ω .* icm2ifs
    return t, ω
end

function equispaced_grid(N_t::Integer, tc::Real)
    t = collect(range(0, tc, length=N_t))
    return t
end

function equispaced_grid(N_ω::Integer, Ω_min::Real, Ω_max::Real)
    ω = collect(range(Ω_min, Ω_max, length=N_ω))
    ω = ω .* icm2ifs
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

# Produce the discrete data of the bath correaltion function.
function discrete_bcf(bcf::Function, t_min::Real, t_max::Real, N_t::Int) :: Tuple{Vector{Float64}, Vector{ComplexF64}}
    t = collect(range(t_min, t_max, length=N_t))
    cc = bcf.(t)
    return t, cc
end