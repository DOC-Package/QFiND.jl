function svd_intermed_decomp_sub(cc::AbstractVector{<:Complex{<:Real}}, 
    U::AbstractMatrix{<:Number}, V::Matrix{ComplexF64}, Dt::AbstractMatrix{<:Number}, rand::Bool=false)
    # Estimate weights via NNLS
    frank = size(U, 2)
    idx, B, err = id_time(U, frank, rand)
    zk, err = linsys_weight(cc, idx, B)
    
    println("Number of basis functions: ", frank)
    println("Error in LS2: ", err)
    return (nsp=frank, basis_time=U, basis_freq=V, coef=zk, Dt=Dt)
end

function derivative_matrix(sv::Vector{Float64}, V::Matrix{ComplexF64}, ω::Vector{Float64})
    N = length(ω)
    dω = abs(ω[2] - ω[1])
    W = Diagonal(ones(ComplexF64, N)) .* dω
    Ω = Diagonal(ω)
    Σ = Diagonal(sv)
    invΣ = Diagonal(1 ./ sv)
    A = V' * (Ω * W) * V
    return Σ * A * invΣ .* (-1im)
end

function derivative_matrix2(sv::Vector{Float64}, V::Matrix{ComplexF64}, ω::Vector{Float64})
    N = length(ω)
    dω = abs(ω[2] - ω[1])
    W = dω .* Diagonal(ones(ComplexF64, N))
    Ω = Diagonal(ω)
    Σ = Diagonal(sv)
    invΣ = Diagonal(1 ./ sv)
    A = V' * Ω * W * V
    return -im * Σ * A * invΣ
end

function svd_intermed_decomp(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    t::AbstractVector{<:Real}, ω::AbstractVector{<:Real}, eps::Real; rand::Bool=false)
    # Create the core matrix f
    fmat = create_integrand_c(sbeta, t, ω)
    # Perform the SVD
    U, sv, V = svd(fmat)
    frank = count(>(eps * sv[1]), sv)
    sv = sv[1:frank]
    U = U[:, 1:frank]
    V = V[:, 1:frank]
    # fit V using AAA algorithm
    fitV = []
    #for i in 1:frank
    #    res = aaa(ω, V[:, i]; tol=1e-10)
    #    fitV[i] = res
    #end
    Dt = derivative_matrix(sv, V, ω)
    
    return svd_intermed_decomp_sub(cc, U, V, Dt, rand)
end

function svd_intermed_decomp(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    t::AbstractVector{<:Real}, ω::AbstractVector{<:Real}, frank::Int; rand::Bool=false)
    # Create the core matrix f
    fmat = create_integrand_c(sbeta, t, ω)
    # Perform the SVD
    U, sv, V = svd(fmat)
    U = U[:, 1:frank]
    V = V[:, 1:frank]
    Dt = derivative_matrix(sv, V, ω)
    
    return svd_intermed_decomp_sub(cc, U, V, Dt, rand)
end

function svd_intermed_decomp(qnsd::Function, bcf::Function, N_t::Integer, N_ω::Integer, 
    tc::Real, Ω_min::Real, Ω_max::Real, eps::Real; rand::Bool=false)
    t, ω = equispaced_grid(Ω_min, Ω_max, tc; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return svd_intermed_decomp(sbeta, cc, t, ω, eps; rand=rand)
end

function svd_intermed_decomp(qnsd::Function, bcf::Function, N_t::Integer, N_ω::Integer, 
    tc::Real, Ω_min::Real, Ω_max::Real, frank::Int; rand::Bool=false)
    t, ω = equispaced_grid(Ω_min, Ω_max, tc; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return svd_intermed_decomp(sbeta, cc, t, ω, frank; rand=rand)
end

function svd_intermed_decomp(dataset::InitialDataSetSVD, eps::Real; rand::Bool=false)
    return svd_intermed_decomp(dataset.qnsd, dataset.bcf, dataset.time, dataset.freq, eps; rand=rand)
end

function svd_intermed_decomp(dataset::InitialDataSetSVD, frank::Int; rand::Bool=false)
    return svd_intermed_decomp(dataset.qnsd, dataset.bcf, dataset.time, dataset.freq, frank; rand=rand)
end
