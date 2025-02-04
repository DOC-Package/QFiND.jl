# Description: Utility functions for the project

# Calculate the sum of the exponential function.
function sumexp(t::Float64, a::AbstractVector{Float64}, c::AbstractVector{Float64}) :: ComplexF64
    return sum(c .* exp.(-im .* a .* t))
end

# Calculate the sum of the exponential function.
function sumexp(t::Float64, a::AbstractVector{ComplexF64}, c::AbstractVector{ComplexF64}) :: ComplexF64
    return sum(c .* exp.(-a .* t))
end