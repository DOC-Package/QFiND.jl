function sort_and_rescale(wk::AbstractVector{<:Real}, gk::AbstractVector{<:Number})
    # Sort frequencies and corresponding coeffficients in ascending order
    perm = sortperm(wk)
    wk = wk[perm]
    gk = gk[perm]
    Nsp = length(wk)
    # Rescale the frequencies and coeffficients
    wk = wk ./ icm2ifs
    gk = gk ./ icm2ifs
    
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
    gk = sqrt(zk .* sbeta[idx[1:frank]])
    Nsp, wk, gk = sort_and_rescale(wk, gk)
    
    println("Error in LS: ", err)
    println("Number of sample points: ", frank)
    
    return (nsp=Nsp, freq=wk, coeff=gk)
end

function id_discr_c(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    ω::AbstractVector{<:Real}, t::AbstractVector{<:Real}, eps::Real; rand::Bool=false)
    # Create the core matrix f.
    fmat = create_integrand_c(sbeta, t, ω)

    # Perform the Interpolative Decomposition (ID).
    frank, idx, B, err1 = id_freq(fmat, eps, rand)
    println("Error in ID: ", err1)
    
    return id_discr_c_sub(sbeta, cc, ω, idx, B, rand)
end

function id_discr_c(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    ω::AbstractVector{<:Real}, t::AbstractVector{<:Real}, frank::Int; rand::Bool=false)
    # Create the core matrix f.
    fmat = create_integrand_c(sbeta, t, ω)
    
    # Perform the Interpolative Decomposition (ID).
    idx, B, err1 = id_freq(fmat, frank, rand)
    
    return id_discr_c_sub(sbeta, cc, ω, idx, B, rand)
end

function id_discr_c(qnsd::Function, bcf::Function, Ω_min::Real, Ω_max::Real, T_max::Real, 
    N_ω::Integer, N_t::Integer, eps::Real; rand::Bool=false)
    t, ω = equispaced_grid(Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return id_discr_c(sbeta, cc, ω, t, eps; rand=rand)
end

function id_discr_c(qnsd::Function, bcf::Function, Ω_min::Real, Ω_max::Real, T_max::Real, 
    N_ω::Integer, N_t::Integer, frank::Int; rand::Bool=false)
    t, ω = equispaced_grid(Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return id_discr_c(sbeta, cc, ω, t, frank; rand=rand)
end

function id_discr_c(dataset::InitialDataSetID, eps::Real; rand::Bool=false)
    return id_discr_c(dataset.qnsd, dataset.bcf, dataset.freq, dataset.time, eps; rand=rand)
end

function id_discr_c(dataset::InitialDataSetID, frank::Int; rand::Bool=false)
    return id_discr_c(dataset.qnsd, dataset.bcf, dataset.freq, dataset.time, frank; rand=rand)
end