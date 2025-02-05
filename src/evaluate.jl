"""
Calculate the error between the approximated and exact correlation functions.
"""
function evaluate_error(
    ωk::AbstractVector{Float64}, 
    gk::AbstractVector{Float64}, 
    bcf::Function, 
    Tc::Real, 
    N_t::Int)

    t = range(0, Tc, length=N_t)
    reference = bcf.(t)
    approx = sumexp.(t, Ref(ωk*icm2ifs), Ref(gk*icm2ifs^2.0))
    error = approx - reference
    
    # Compute the normalized errors.
    norm_val = abs(bcf(0.0))
    normalized_max_error = maximum(abs.(error)) / norm_val
    normalized_avg_error = sum(abs.(error)) / (norm_val * N_t)
    
    println("Normalized maximum error: ", normalized_max_error)
    println("Normalized mean absolute error: ", normalized_avg_error)
end

function evaluate_error(
    ωk::AbstractVector{Float64}, 
    gk::AbstractVector{Float64}, 
    bcf::Function, 
    t::AbstractVector{<:Real})

    reference = bcf.(t)
    approx = sumexp.(t, Ref(ωk*icm2ifs), Ref(gk*icm2ifs^2.0))
    error = approx - reference
    
    # Compute the normalized errors.
    norm_val = abs(bcf(0.0))
    normalized_max_error = maximum(abs.(error)) / norm_val
    normalized_avg_error = sum(abs.(error)) / (norm_val * N_t)
    
    println("Normalized maximum error: ", normalized_max_error)
    println("Normalized mean absolute error: ", normalized_avg_error)
end

function evaluate_error(
    ωk::AbstractVector{Float64}, 
    gk::AbstractVector{Float64}, 
    reference::AbstractVector{Complex{<:Real}},
    t::AbstractVector{<:Real})

    reference = bcf.(t)
    approx = sumexp.(t, Ref(ωk*icm2ifs), Ref(gk*icm2ifs^2.0))
    error = approx - reference
    
    # Compute the normalized errors.
    norm_val = abs(bcf(0.0))
    normalized_max_error = maximum(abs.(error)) / norm_val
    normalized_avg_error = sum(abs.(error)) / (norm_val * N_t)
    
    println("Normalized maximum error: ", normalized_max_error)
    println("Normalized mean absolute error: ", normalized_avg_error)
end
