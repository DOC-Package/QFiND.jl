"""
    create_integrand(t, w)

Construct the matrix `f` to be decomposed. The matrix has size `(2*N_t, N_ω)`.

# Arguments
- `sbeta::Function`: Quantum noise spectral density.
- `t::Vector{<:Real}`: Time grid (length `N_t`).
- `w::Vector{<:Real}`: Frequency grid (length `N_ω`).

# Returns
A real matrix `f` of size `(2*N_t, N_ω)`.
"""
function create_integrand(sbeta::Function, t::Vector{<:Real}, ω::Vector{<:Real})
    N_t = length(t)
    N_ω = length(ω)
    f = zeros(Float64, 2*N_t, N_ω)
    # Fill first N_t rows (real part)
    for i in 1:N_t
        for j in 1:N_ω
            f[i, j] = sbeta(ω[j]; scale=icm2ifs) * cos(ω[j] * t[i])
        end
    end
    # Fill next N_t rows (imaginary part)
    for i in (N_t+1):(2*N_t)
        for j in 1:N_ω
            f[i, j] = -sbeta(ω[j]; scale=icm2ifs) * sin(ω[j] * t[i - N_t])
        end
    end

    return f
end

function create_integrand(sbeta::Vector{<:Real}, t::Vector{<:Real}, ω::Vector{<:Real}) :: Matrix{Float64}
    N_t = length(t)
    N_ω = length(ω)
    f = zeros(Float64, 2*N_t, N_ω)
    # Fill first N_t rows (real part)
    for i in 1:N_t
        for j in 1:N_ω
            f[i, j] = sbeta[j] * cos(ω[j] * t[i])
        end
    end
    # Fill next N_t rows (imaginary part)
    for i in (N_t+1):(2*N_t)
        for j in 1:N_ω
            f[i, j] = -sbeta[j] * sin(ω[j] * t[i - N_t])
        end
    end

    return f
end

"""
    nnls_weight(t, B)

Estimate the coefficients (amplitudes) by solving a nonnegative least‐squares (NNLS)
problem. 

# Arguments
- `bcf::Function`: Function that returns the BCF at a given time.
- `t::Vector{<:Real}`: Time grid (length `N_t`).
- `B::AbstractMatrix{<:Real}`: A matrix of size `(2*N_t, r)`, where `r` is the ID rank.

# Returns
- `g` is the vector of estimated coefficients (of length `r`),
- `err` is the residual norm from the NNLS solve.
"""
function nnls_weight(bcf::Function, t::Vector{<:Real}, B::AbstractMatrix{<:Real}) 
    N_t = length(t)
    cc = bcf.(t)
    # Construct the vector c = [real(cc); imag(cc)].
    c = zeros(Float64, 2*N_t)
    c[1:N_t] .= real.(cc)
    c[N_t+1:2*N_t] .= imag.(cc)
    # Solve the NNLS problem: minimize ||B*g - c|| subject to g ≥ 0.
    g = nonneg_lsq(B, c; alg=:nnls) |> vec
    err = norm(B * g - c)
    return g, err
end

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
    g = nonneg_lsq(B, c; alg=:nnls) |> vec
    err = norm(B * g - c)
    return g, err
end

function sort_and_rescale(wk::AbstractVector{<:Real}, gk::AbstractVector{<:Real})
    # Sort frequencies and corresponding coefficients in ascending order
    perm = sortperm(wk)
    wk = wk[perm]
    gk = gk[perm]
    
    Nsp = length(wk)
    # Remove any zero coefficients.
    if minimum(gk) == 0.0
        keep = findall(x -> x > 0.0, gk)
        wk = wk[keep]
        gk = gk[keep]
        Nsp = length(wk)
    end

    # Rescale the frequencies and coefficients
    wk = wk ./ icm2ifs
    gk = gk ./ icm2ifs^2.0
    
    return Nsp, wk, gk
end

"""
    id_discr(N_t, N_ω, tc, Ω_min, Ω_max, eps, frank; rand=false)

Estimate the frequencies and coefficients using interpolative decomposition (ID) and NNLS.

# Arguments
- `sbeta::Function`: Quantum noise spectral density.
- `N_t::Integer`: Number of time points.
- `N_ω::Integer`: Number of frequency points.
- `tc::Real`: Maximum time value.
- `Ω_min::Real`: Minimum frequency value.
- `Ω_max::Real`: Maximum frequency value.
- `eps::Real`: Error tolerance for ID (used when `frank < 1`).
- `frank`: Desired rank for the ID. If `frank < 1`, the algorithm uses `eps` to choose the rank.
- `rand::Bool` (optional): Whether to use a randomized ID algorithm. Default is `false`.

# Returns
A tuple `(Nsp, wk, zk, frank)` where:
- `Nsp` is the number of estimated frequencies,
- `wk` is the vector of estimated frequencies,
- `zk` is the vector of estimated coefficients,
- `frank` is the rank actually used.
"""
function id_discr(dataset::DiscreteDataSetID, eps::Real; rand::Bool=false)
   
    # Create the core matrix f.
    fmat = create_integrand(dataset.qnsd, dataset.time, dataset.freq)
    
    # Perform the Interpolative Decomposition (ID).
    println("Starting ID")
    frank, idx, B, err1 = id_freq(fmat, eps, rand)
    println("Rank of f: ", frank)
    
    # Estimated frequencies: use the first frank_final indices.
    wk = dataset.freq[idx[1:frank]]
    
    # Estimate weights via NNLS
    zk, err2 = nnls_weight(dataset.bcf, dataset.time, B)
    gk = zk .* dataset.qnsd[idx[1:frank]]
    println(typeof(zk))
    println(typeof(gk))

    Nsp, wk, gk = sort_and_rescale(wk, gk)
    
    println("Number of sample points: ", Nsp)
    println("Error in ID: ", err1)
    println("Error in NNLS: ", err2)
    
    return (nsp=Nsp, freq=wk, coef=gk, frank=frank)
end

function id_discr(sbeta::Function, bcf::Function, t::AbstractVector{<:Real}, ω::AbstractVector{<:Real}, eps::Real; rand::Bool=false)
    
    # Create the core matrix f.
    fmat = create_integrand(sbeta, t, ω)
    
    # Perform the Interpolative Decomposition (ID).
    frank, idx, B, err1 = id_freq(fmat, eps, rand)
    println("Rank of f: ", frank)
    
    # Estimated frequencies: use the first frank_final indices.
    wk = ω[idx[1:frank]]
    
    # Estimate weights via NNLS
    zk, err2 = nnls_weight(bcf, t, B)
    gk = zk .* sbeta.(wk; scale=icm2ifs)
    
    Nsp, wk, gk = sort_and_rescale(wk, gk)
    
    println("Number of sample points: ", Nsp)
    println("Error in ID: ", err1)
    println("Error in NNLS: ", err2)
    
    return (nsp=Nsp, freq=wk, coef=gk)
end

function id_discr(sbeta::Function, bcf::Function, t::AbstractVector{<:Real}, ω::AbstractVector{<:Real}, frank::Int; rand::Bool=false)
    
    # Create the core matrix f.
    fmat = create_integrand(sbeta, t, ω)
    
    # Perform the Interpolative Decomposition (ID).
    idx, B, err1 = id_freq(fmat, frank, rand)
    println("Rank of f: ", frank)
    
    # Estimated frequencies: use the first frank_final indices.
    wk = ω[idx[1:frank]]
    
    # Estimate weights via NNLS
    zk, err2 = nnls_weight(bcf, t, B)
    gk = zk .* sbeta.(wk; scale=icm2ifs)
    
    Nsp, wk, gk = sort_and_rescale(wk, gk)
    
    println("Number of sample points: ", Nsp)
    println("Error in ID: ", err1)
    println("Error in NNLS: ", err2)
    
    return (nsp=Nsp, freq=wk, coef=gk)
end

function id_discr(sbeta::Function, bcf::Function, N_t::Integer, N_ω::Integer, tc::Real, Ω_min::Real, Ω_max::Real, eps::Real; rand::Bool=false)
    # Generate equispaced time and frequency grids.
    t, ω = equispaced_grid(Ω_min, Ω_max, tc; n_freq=N_ω, n_time=N_t, scale=icm2ifs)
    
    return id_discr(sbeta, bcf, t, ω, eps; rand=rand)
end

function id_discr(sbeta::Function, bcf::Function, N_t::Integer, N_ω::Integer, tc::Real, Ω_min::Real, Ω_max::Real, frank::Int; rand::Bool=false)
    # Generate equispaced time and frequency grids.
    t, ω = equispaced_grid(Ω_min, Ω_max, tc; n_freq=N_ω, n_time=N_t, scale=icm2ifs)

    return id_discr(sbeta, bcf, t, ω, frank; rand=rand)
end
