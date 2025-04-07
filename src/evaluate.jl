"""
Calculate the error between the approximated and exact correlation functions.
"""
function evaluate_error(
    t::AbstractVector{<:Real},
    approx::AbstractVector{<:Number},
    reference::AbstractVector{<:Number} 
    )

    error = approx - reference    
    # Compute the normalized errors.
    N_t = length(t)
    norm_val = abs(reference[1])
    normalized_max_error = maximum(abs.(error)) / norm_val
    normalized_avg_error = sum(abs.(error)) / (norm_val * N_t)
    
    println("Normalized maximum error: ", normalized_max_error)
    println("Normalized mean absolute error: ", normalized_avg_error)
end

function evaluate_error(
    ω::AbstractVector{Float64}, 
    g::AbstractVector{<:Number}, 
    bcf::Function, 
    T_max::Real, 
    N_t::Int)

    t = range(0, T_max, length=N_t)
    evaluate_error(ω, g, t, bcf)
end

function evaluate_error(
    ω::AbstractVector{Float64}, 
    g::AbstractVector{<:Number}, 
    t::AbstractVector{<:Real},
    bcf::Function
    )

    evaluate_error(ω, g, t, bcf.(t))
end

function evaluate_error(
    ω::AbstractVector{Float64}, 
    g::AbstractVector{<:Number}, 
    t::AbstractVector{<:Real},
    reference::AbstractVector{<:Number}
    )

    approx = bcf_approx.(t, Ref(ω), Ref(g))
    evaluate_error(t, approx, reference)
end

function evaluate_error(
    U::AbstractMatrix{<:Number}, 
    zk::AbstractVector{<:Number}, 
    reference::AbstractVector{ComplexF64})

    approx = U * zk
    error = approx - reference
    
    # Compute the normalized errors.
    norm_val = maximum(abs.(reference))
    normalized_max_error = maximum(abs.(error)) / norm_val
    normalized_avg_error = sum(abs.(error)) / (norm_val * size(reference,1))
    
    println("Normalized maximum error: ", normalized_max_error)
    println("Normalized mean absolute error: ", normalized_avg_error)
end
