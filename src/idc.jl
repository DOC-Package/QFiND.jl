function sort_and_rescale(wk::AbstractVector{<:Real}, gk::AbstractVector{<:Number})
    # Sort frequencies and corresponding coefficients in ascending order
    perm = sortperm(wk)
    wk = wk[perm]
    gk = gk[perm]
    Nsp = length(wk)
    # Rescale the frequencies and coefficients
    wk = wk ./ icm2ifs
    gk = gk ./ icm2ifs^2.0
    
    return Nsp, wk, gk
end

function id_discr_c_sub(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    ω::AbstractVector{<:Real}, idx::Vector{Int}, B::AbstractMatrix{<:Number}, rand::Bool=false)
    # Estimated frequencies: use the first frank_final indices.
    frank = size(B, 2)
    wk = ω[idx[1:frank]]
    
    # Estimate weights via NNLS
    idx1, D, err = id_time(B, frank, rand)
    zk, err = linsys_weight(cc, idx1, D)
    gk = zk .* sbeta[idx[1:frank]]
    Nsp, wk, gk = sort_and_rescale(wk, gk)
    
    println("Error in LS: ", err)
    println("Number of sample points: ", frank)
    
    return (nsp=Nsp, freq=wk, coef=gk)
end

function id_discr_c(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    t::AbstractVector{<:Real}, ω::AbstractVector{<:Real}, eps::Real; rand::Bool=false)
    # Create the core matrix f.
    fmat = create_integrand_c(sbeta, t, ω)

    # Perform the Interpolative Decomposition (ID).
    frank, idx, B, err1 = id_freq(fmat, eps, rand)
    println("Error in ID: ", err1)
    
    return id_discr_c_sub(sbeta, cc, ω, idx, B, rand)
end

function id_discr_c(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    t::AbstractVector{<:Real}, ω::AbstractVector{<:Real}, frank::Int; rand::Bool=false)
    # Create the core matrix f.
    fmat = create_integrand_c(sbeta, t, ω)
    
    # Perform the Interpolative Decomposition (ID).
    idx, B, err1 = id_freq(fmat, frank, rand)
    
    return id_discr_c_sub(sbeta, cc, ω, idx, B, rand)
end

function id_discr_c(qnsd::Function, bcf::Function, N_t::Integer, N_ω::Integer, 
    tc::Real, Ω_min::Real, Ω_max::Real, eps::Real; rand::Bool=false)
    t, ω = equispaced_grid(Ω_min, Ω_max, tc; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return id_discr_c(sbeta, cc, t, ω, eps; rand=rand)
end

function id_discr_c(qnsd::Function, bcf::Function, N_t::Integer, N_ω::Integer, 
    tc::Real, Ω_min::Real, Ω_max::Real, frank::Int; rand::Bool=false)
    t, ω = equispaced_grid(Ω_min, Ω_max, tc; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return id_discr_c(sbeta, cc, t, ω, frank; rand=rand)
end

function id_discr_c(dataset::InitialDataSetID, eps::Real; rand::Bool=false)
    return id_discr_c(dataset.qnsd, dataset.bcf, dataset.time, dataset.freq, eps; rand=rand)
end

function id_discr_c(dataset::InitialDataSetID, frank::Int; rand::Bool=false)
    return id_discr_c(dataset.qnsd, dataset.bcf, dataset.time, dataset.freq, frank; rand=rand)
end