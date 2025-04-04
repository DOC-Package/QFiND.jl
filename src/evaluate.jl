"""
Calculate the error between the approximated and exact correlation functions.
"""
function evaluate_error(
    approx::AbstractVector{<:Number},
    reference::AbstractVector{<:Number}, 
    t::AbstractVector{<:Real})

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
    ωk::AbstractVector{Float64}, 
    gk::AbstractVector{<:Number}, 
    bcf::Function, 
    Tc::Real, 
    N_t::Int)

    t = range(0, Tc, length=N_t)
    evaluate_error(ωk, gk, bcf, t)
end

function evaluate_error(
    ωk::AbstractVector{Float64}, 
    gk::AbstractVector{<:Number}, 
    bcf::Function, 
    t::AbstractVector{<:Real})

    evaluate_error(ωk, gk, bcf.(t), t)
end

function evaluate_error(
    ωk::AbstractVector{Float64}, 
    gk::AbstractVector{<:Number}, 
    reference::AbstractVector{<:Number}, 
    t::AbstractVector{<:Real})

    approx = bcf_approx.(t, Ref(ωk*icm2ifs), Ref(gk*icm2ifs^2.0))
    evaluate_error(approx, reference, t)
end

function evaluate_error(
    ak::AbstractVector{ComplexF64}, 
    ck::AbstractVector{ComplexF64}, 
    bcf::Function, 
    Tc::Real, 
    N_t::Int)

    t = range(0, Tc, length=N_t)
    reference = bcf.(t)
    approx = bcf_approx.(t, Ref(ak*icm2ifs), Ref(ck*icm2ifs^2.0))
    error = approx - reference
    
    # Compute the normalized errors.
    norm_val = abs(bcf(0.0))
    normalized_max_error = maximum(abs.(error)) / norm_val
    normalized_avg_error = sum(abs.(error)) / (norm_val * N_t)
    
    println("Normalized maximum error: ", normalized_max_error)
    println("Normalized mean absolute error: ", normalized_avg_error)
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
