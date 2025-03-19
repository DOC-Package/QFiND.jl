function calc_coeff(cc::AbstractVector{<:Complex{<:Real}}, U::AbstractMatrix{<:Number})
    frank = size(U, 2)
    idx, B, err = id_time(U, frank, false)
    z, err = linsys_weight(cc, idx, B)
    return z
end

function calc_coeff(ω::AbstractVector{<:Real}, sv::Vector{Float64}, V::Matrix{ComplexF64})
    dω = abs(ω[2] - ω[1])
    W = ones(ComplexF64, length(ω)) .* dω
    z = V' * W
    z = z .* sv
    return z
end

function calc_derivative_matrix(sv::Vector{Float64}, V::Matrix{ComplexF64}, ω::Vector{Float64})
    Ω = Diagonal(ω)
    Σ = Diagonal(sv)
    invΣ = Diagonal(1 ./ sv)
    A = V' * (Ω * V)
    return Σ * A * invΣ .* (-1im)
end

function svd_intermed_decomp(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    t::AbstractVector{<:Real}, ω::AbstractVector{<:Real}, eps::Real)
    fmat = create_integrand_c(sbeta, t, ω)
    # Perform the SVD
    U, sv, V = svd(fmat)
    frank = count(>(eps * sv[1]), sv)
    sv = sv[1:frank]
    U = U[:, 1:frank]
    V = V[:, 1:frank]
    z = calc_coeff(cc, U)
    Dt = calc_derivative_matrix(sv, V, ω)
    println("Number of basis functions: ", frank)
    return (nsp=frank, basis_time=U, basis_freq=V, coef=z, Dt=Dt)
end

function svd_intermed_decomp(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    t::AbstractVector{<:Real}, ω::AbstractVector{<:Real}, frank::Int)
    fmat = create_integrand_c(sbeta, t, ω)
    # Perform the SVD
    U, sv, V = svd(fmat)
    U = U[:, 1:frank]
    V = V[:, 1:frank]
    z = calc_coeff(cc, U)
    Dt = derivative_matrix(sv, V, ω)
    return (nsp=frank, basis_time=U, basis_freq=V, coef=z, Dt=Dt)
end

function svd_intermed_decomp(qnsd::Function, bcf::Function, N_t::Integer, N_ω::Integer, 
    tc::Real, Ω_min::Real, Ω_max::Real, eps::Real)
    t, ω = equispaced_grid(Ω_min, Ω_max, tc; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return svd_intermed_decomp(sbeta, cc, t, ω, eps)
end

function svd_intermed_decomp(qnsd::Function, bcf::Function, N_t::Integer, N_ω::Integer, 
    tc::Real, Ω_min::Real, Ω_max::Real, frank::Int)
    t, ω = equispaced_grid(Ω_min, Ω_max, tc; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return svd_intermed_decomp(sbeta, cc, t, ω, frank)
end

function svd_intermed_decomp(dataset::InitialDataSetSVD, eps::Real)
    return svd_intermed_decomp(dataset.qnsd, dataset.bcf, dataset.time, dataset.freq, eps)
end

function svd_intermed_decomp(dataset::InitialDataSetSVD, frank::Int)
    return svd_intermed_decomp(dataset.qnsd, dataset.bcf, dataset.time, dataset.freq, frank)
end
