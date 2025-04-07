
"""
    nnls_weight(t, B)

Estimate the coefficients (amplitudes) by solving a nonnegative least‐squares (NNLS)
problem. 
"""
function nnls_weight(cc::AbstractVector{<:Complex{<:Real}}, t::Vector{<:Real}, B::AbstractMatrix{<:Real}) 
    N_t = length(t)
    # check if cc has the correct length
    if length(cc) != N_t
        throw(ArgumentError("The length of the BCF vector must be equal to the length of the time vector."))
    end
    # Construct the vector c = [real(cc); imag(cc)].
    c = zeros(Float64, 2*N_t)
    c[1:N_t] .= real.(cc)
    c[N_t+1:2*N_t] .= imag.(cc)
    # Solve the NNLS problem: minimize ||B*g - c|| subject to g ≥ 0.
    z = nonneg_lsq(B, c; alg=:nnls) |> vec
    err = norm(B * z - c)
    return z, err
end

function nnls_weight_time(cc::AbstractVector{<:Complex{<:Real}}, t::Vector{<:Real}, B::AbstractMatrix{<:Real})
    N_t = length(t)
    # check if cc has the correct length
    if length(cc) != N_t
        throw(ArgumentError("The length of the BCF vector must be equal to the length of the time vector."))
    end
    frank = size(B, 2)
    idx, D, err = id_time(B, frank, false)
    # Construct the vector c = [real(cc); imag(cc)].
    c = vcat(real.(cc), imag.(cc))
    c = c[idx[1:frank]]
    # Solve the NNLS problem: minimize ||B*g - c|| subject to g ≥ 0.
    z = nonneg_lsq(D, c; alg=:nnls) |> vec
    err = norm(D * z - c)
    return z, err
end

function sort_and_rescale(wk::AbstractVector{<:Real}, zk::AbstractVector{<:Real}, gk::AbstractVector{<:Real})
    # Sort frequencies and corresponding coefficients in ascending order
    perm = sortperm(wk)
    wk = wk[perm]
    zk = zk[perm]
    gk = gk[perm]
    
    Nsp = length(wk)
    # Remove any zero coefficients.
    if minimum(gk) == 0.0
        keep = findall(x -> x > 0.0, gk)
        wk = wk[keep]
        zk = zk[keep]
        gk = gk[keep]
        Nsp = length(wk)
    end

    # Rescale the frequencies and coefficients
    wk = wk ./ icm2ifs
    zk = zk ./ icm2ifs
    gk = gk ./ icm2ifs
    
    return Nsp, wk, zk, gk
end

function id_discr_sub(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}},
    ω::AbstractVector{<:Real}, t::AbstractVector{<:Real}, idx::Vector{Int}, B::AbstractMatrix{<:Real})
    # Estimated frequencies: use the first frank_final indices.
    frank = size(B, 2)
    wk = ω[idx[1:frank]]
    
    # Estimate weights via NNLS
    zk, err2 = nnls_weight(cc, t, B)
    gk = sqrt.(zk .* sbeta[idx[1:frank]] .* (2.0 / π))
    Nsp, wk, zk, gk = sort_and_rescale(wk, zk, gk)
    
    println("Error in NNLS: ", err2)
    println("Number of sample points: ", Nsp)
    
    return (nsp=Nsp, freq=wk, coeff=gk, weight=zk, frank=frank)
end

function id_discr(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    ω::AbstractVector{<:Real}, t::AbstractVector{<:Real}, eps::Real; rand::Bool=false)
    # Create the core matrix f.
    fmat = create_integrand(sbeta, t, ω)

    # Perform the Interpolative Decomposition (ID).
    frank, idx, B, err1 = id_freq(fmat, eps, rand)
    println("Error in ID: ", err1)
    
    return id_discr_sub(sbeta, cc, ω, t, idx, B)
end

function id_discr(sbeta::AbstractVector{<:Real}, cc::AbstractVector{<:Complex{<:Real}}, 
    ω::AbstractVector{<:Real}, t::AbstractVector{<:Real}, frank::Int; rand::Bool=false)
    # Create the core matrix f.
    fmat = create_integrand(sbeta, t, ω)
    
    # Perform the Interpolative Decomposition (ID).
    idx, B, err1 = id_freq(fmat, frank, rand)
    
    return id_discr_sub(sbeta, cc, ω, t, idx, B)
end

function id_discr(qnsd::Function, bcf::Function, Ω_min::Real, Ω_max::Real, T_max::Real, 
    N_ω::Integer, N_t::Integer, eps::Real; rand::Bool=false)
    t, ω = equispaced_grid(Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return id_discr(sbeta, cc, ω, t, eps; rand=rand)
end

function id_discr(qnsd::Function, bcf::Function, Ω_min::Real, Ω_max::Real, T_max::Real, 
    N_ω::Integer, N_t::Integer, frank::Int; rand::Bool=false)
    t, ω = equispaced_grid(Ω_min, Ω_max, T_max; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    sbeta = qnsd.(ω; scale=icm2ifs)
    cc = bcf.(t)
    return id_discr(sbeta, cc, ω, t, frank; rand=rand)
end

function id_discr(dataset::InitialDataSetID, eps::Real; rand::Bool=false)
    return id_discr(dataset.qnsd, dataset.bcf, dataset.freq, dataset.time, eps; rand=rand)
end

function id_discr(dataset::InitialDataSetID, frank::Int; rand::Bool=false)
    return id_discr(dataset.qnsd, dataset.bcf, dataset.freq, dataset.time, frank; rand=rand)
end