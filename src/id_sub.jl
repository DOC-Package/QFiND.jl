"""
    id_freq_eps(f::AbstractMatrix{Float64}, eps::Real, rnd::Bool)

Perform interpolative decomposition on matrix `f` with error tolerance `eps`.

# Arguments
- `f::AbstractMatrix`: Input matrix (of size 2M×N).
- `eps::Real`: Error tolerance for the decomposition.
- `rnd::Bool`: If `true`, use the randomized variant.

# Returns
A tuple `(frank, idx, B, err)` where
- `frank::Int` is the computed rank of the approximation.
- `idx::Vector{Int}` contains the indices of the selected columns.
- `B::AbstractMatrix` is the skeleton matrix (the selected columns of `f`).
- `err::Real` is the maximum absolute error between `f` and its reconstruction.
"""
function id_freq(f::AbstractMatrix{Float64}, eps::Real, rnd::Bool)
    
    # Compute the interpolative decomposition using the tolerance mode.
    opts = LRAOptions(rtol=eps, sketch=:none)
    id_obj = idfact(f, opts)
    idx = id_obj[:sk]
    frank = length(idx)
    
    # The skeleton matrix B is simply the columns of f0 indexed by idx.
    B = f[:, idx]

    # The reconstructed approximation.
    f1 = ID(f, id_obj)

    # Compute the error.
    err = norm(Matrix(f1) - f)
    
    return frank, idx, B, err
end

"""
    id_freq_rank(f::AbstractMatrix{Float64}, frank::Int, rnd::Bool)

Perform interpolative decomposition on matrix `f` with a specified rank `frank`.

# Arguments
- `f::AbstractMatrix`: Input matrix (of size 2M×N).
- `frank::Integer`: Desired rank for the approximation.
- `rnd::Bool`: If `true`, use the randomized variant.

# Returns
A tuple `(idx, B, err)` where
- `idx::Vector{Int}` contains the indices of the selected columns.
- `B::AbstractMatrix` is the skeleton matrix (the selected columns of `f`).
- `err::Real` is the maximum absolute error between `f` and its reconstruction.
"""
function id_freq(f::AbstractMatrix{Float64}, frank::Int, rnd::Bool)

    # Compute the interpolative decomposition using the tolerance mode.
    opts = LRAOptions(rank=frank, sketch=:none)
    id_obj = idfact(f, opts)
    idx = id_obj[:sk]
    frank = length(idx)
    
    # The skeleton matrix B is simply the columns of f0 indexed by idx.
    B = f[:, idx]

    # The reconstructed approximation.
    f1 = ID(f, id_obj)

    # Compute the error.
    err = norm(Matrix(f1) - f)
    
    return idx, B, err
end
