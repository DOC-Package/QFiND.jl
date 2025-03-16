function svd_basis_decomp_sub(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    ω::AbstractVector{<:Real}, U::AbstractMatrix{<:Number}, rand::Bool=false)
    # Estimate weights via NNLS
    frank = size(U, 2)
    idx, D, err = id_time(U, frank, rand)
    zk, err = linsys_weight(cc, idx, D)
    
    println("Number of basis functions: ", frank)
    println("Error in LS2: ", err)
    
    return (nsp=frank, basis=U, coef=zk)
end

function svd_basis_decomp(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    t::AbstractVector{<:Real}, ω::AbstractVector{<:Real}, eps::Real; rand::Bool=false)
    # Create the core matrix f.
    fmat = create_integrand_c(sbeta, t, ω)
    # Perform the SVD
    U, sv, V = svd(fmat)
    frank = count(>(eps * sv[1]), sv)
    U = U[:, 1:frank]
    V = V[:, 1:frank]
    
    return svd_basis_decomp_sub(sbeta, cc, ω, U, rand)
end

function svd_basis_decomp(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    t::AbstractVector{<:Real}, ω::AbstractVector{<:Real}, frank::Int; rand::Bool=false)
    # Create the core matrix f.
    fmat = create_integrand_c(sbeta, t, ω)
    # Perform the SVD
    U, sv, V = svd(fmat)
    U = U[:, 1:frank]
    V = V[:, 1:frank]
    
    return svd_basis_decomp_sub(sbeta, cc, ω, U, rand)
end

function svd_basis_decomp(qnsd::Function, bcf::Function, N_t::Integer, N_ω::Integer, 
    tc::Real, Ω_min::Real, Ω_max::Real, eps::Real; rand::Bool=false)
    t, ω = equispaced_grid(Ω_min, Ω_max, tc; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return svd_basis_decomp(sbeta, cc, t, ω, eps; rand=rand)
end

function svd_basis_decomp(qnsd::Function, bcf::Function, N_t::Integer, N_ω::Integer, 
    tc::Real, Ω_min::Real, Ω_max::Real, frank::Int; rand::Bool=false)
    t, ω = equispaced_grid(Ω_min, Ω_max, tc; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return svd_basis_decomp(sbeta, cc, t, ω, frank; rand=rand)
end

function svd_basis_decomp(dataset::InitialDataSetID, eps::Real; rand::Bool=false)
    return svd_basis_decomp(dataset.qnsd, dataset.bcf, dataset.time, dataset.freq, eps; rand=rand)
end

function svd_basis_decomp(dataset::InitialDataSetID, frank::Int; rand::Bool=false)
    return svd_basis_decomp(dataset.qnsd, dataset.bcf, dataset.time, dataset.freq, frank; rand=rand)
end
