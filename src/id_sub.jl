"""
    id_freq_eps(f::AbstractMatrix{Float64}, eps::Real, rnd::Bool)

Perform interpolative decomposition on matrix `f` with error tolerance `eps`.
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

"""
    id_freq_eps(f::AbstractMatrix{Float64}, eps::Real, rnd::Bool)

Perform interpolative decomposition on matrix `f` with error tolerance `eps`.
"""
function id_freq(f::AbstractMatrix{ComplexF64}, eps::Real, rnd::Bool)
    
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
    id_freq_rank(f::AbstractMatrix{ComplexF64}, frank::Int, rnd::Bool)

Perform interpolative decomposition on matrix `f` with a specified rank `frank`.
"""
function id_freq(f::AbstractMatrix{ComplexF64}, frank::Int, rnd::Bool)

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


"""
    id_time(f::AbstractMatrix{ComplexF64}, frank::Int, rnd::Bool)

Perform interpolative decomposition on matrix `f` with a specified rank `frank`.
"""
function id_time(B::AbstractMatrix{Float64}, brank::Int, rnd::Bool)

    # Compute the interpolative decomposition using the tolerance mode.
    opts = LRAOptions(rank=brank, sketch=:none)
    BT = transpose(B)
    id_obj = idfact(BT, opts)
    idx = id_obj[:sk]
    
    # The skeleton matrix B is simply the columns of f0 indexed by idx.
    D = Matrix(transpose(BT[:, idx]))

    # The reconstructed approximation.
    BT1 = ID(BT, id_obj)

    # Compute the error.
    err = norm(Matrix(BT1) - BT)
    
    return idx, D, err
end

function id_time(B::AbstractMatrix{ComplexF64}, brank::Int, rnd::Bool)

    # Compute the interpolative decomposition using the tolerance mode.
    opts = LRAOptions(rank=brank, sketch=:none)
    BT = transpose(B)
    id_obj = idfact(BT, opts)
    idx = id_obj[:sk]
    
    # The skeleton matrix B is simply the columns of f0 indexed by idx.
    D = Matrix(transpose(BT[:, idx]))

    # The reconstructed approximation.
    BT1 = ID(BT, id_obj)

    # Compute the error.
    err = norm(Matrix(BT1) - BT)
    
    return idx, D, err
end
