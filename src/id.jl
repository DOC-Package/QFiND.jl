"""
    equispaced_mesh(N_t, N_w, tc, omega_min, omega_max)

Generate equispaced time and frequency grids.

# Arguments
- `N_t::Integer`: Number of time points.
- `N_w::Integer`: Number of frequency points.
- `tc::Real`: Maximum time value.
- `omega_min::Real`: Minimum frequency value.
- `omega_max::Real`: Maximum frequency value.

# Returns
A tuple `(t, w)` where
- `t` is a vector of time points (length `N_t`),
- `w` is a vector of frequency points (length `N_w`), scaled by `icm2ifs`.
"""
function equispaced_mesh(N_t::Integer, N_w::Integer, tc::Real, omega_min::Real, omega_max::Real)
    t = collect(range(0, tc, length=N_t))
    w = collect(range(omega_min, omega_max, length=N_w))
    w = [ω * icm2ifs for ω in w]
    return t, w
end

"""
    create_integrand(t, w)

Construct the matrix `f` to be decomposed. The matrix has size `(2*N_t, N_w)`.

# Arguments
- `sbeta::Function`: Quantum noise spectral density.
- `t::Vector{<:Real}`: Time grid (length `N_t`).
- `w::Vector{<:Real}`: Frequency grid (length `N_w`).

# Returns
A real matrix `f` of size `(2*N_t, N_w)`.
"""
function create_integrand(sbeta::Function, t::Vector{<:Real}, w::Vector{<:Real})
    N_t = length(t)
    N_w = length(w)
    f = zeros(Float64, 2*N_t, N_w)
    # Fill first N_t rows (real part)
    for i in 1:N_t
        for j in 1:N_w
            f[i, j] = sbeta(w[j]; scale=icm2ifs) * cos(w[j] * t[i])
        end
    end
    # Fill next N_t rows (imaginary part)
    for i in (N_t+1):(2*N_t)
        for j in 1:N_w
            f[i, j] = -sbeta(w[j]; scale=icm2ifs) * sin(w[j] * t[i - N_t])
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
    g = nonneg_lsq(B, c; alg=:nnls)
    err = norm(B * g - c)
    return g, err
end

"""
    id_discr(N_t, N_w, tc, omega_min, omega_max, eps, frank; rand=false)

Estimate the frequencies and coefficients using interpolative decomposition (ID) and NNLS.

# Arguments
- `sbeta::Function`: Quantum noise spectral density.
- `N_t::Integer`: Number of time points.
- `N_w::Integer`: Number of frequency points.
- `tc::Real`: Maximum time value.
- `omega_min::Real`: Minimum frequency value.
- `omega_max::Real`: Maximum frequency value.
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
function id_discr(sbeta::Function, bcf::Function, N_t::Integer, N_w::Integer, tc::Real, omega_min::Real, omega_max::Real, eps::Real; rand::Bool=false)
    # Generate equispaced time and frequency meshes.
    t, w = equispaced_mesh(N_t, N_w, tc, omega_min, omega_max)

    # Create the core matrix f.
    fmat = create_integrand(sbeta, t, w)
    
    # Perform the Interpolative Decomposition (ID).
    println("Starting ID")
    frank, idx, B, err1 = id_freq(fmat, eps, rand)
    println("Rank of f: ", frank)
    
    # Estimated frequencies: use the first frank_final indices.
    wk = w[idx[1:frank]]
    
    # Estimate weights via NNLS
    zk, err2 = nnls_weight(bcf, t, B)
    
    # Sort frequencies and corresponding coefficients in ascending order
    perm = sortperm(wk)
    wk = wk[perm]
    zk = zk[perm]
    
    Nsp = frank
    # Remove any zero coefficients.
    if minimum(zk) == 0.0
        keep = findall(x -> x > 0.0, zk)
        wk = wk[keep]
        zk = zk[keep]
        Nsp = length(wk)
    end

    # Coefficients
    gk = zk .* sbeta.(wk; scale=icm2ifs)

    # Rescale the frequencies and coefficients
    wk = wk ./ icm2ifs
    gk = gk ./ icm2ifs^2.0
    
    println("Number of sample points: ", Nsp)
    println("Error in ID: ", err1)
    println("Error in NNLS: ", err2)
    
    return (nsp=Nsp, freq=wk, coef=gk, frank=frank)
end

function id_discr(sbeta::Function, bcf::Function, N_t::Integer, N_w::Integer, tc::Real, omega_min::Real, omega_max::Real, frank::Int; rand::Bool=false)
    # Generate equispaced time and frequency meshes.
    t, w = equispaced_mesh(N_t, N_w, tc, omega_min, omega_max)
    
    # Create the core matrix f.
    fmat = create_integrand(sbeta, t, w)
    
    # Perform the Interpolative Decomposition (ID).
    idx, B, err1 = id_freq(fmat, frank, rand)
    println("Rank of f: ", frank)
    
    # Estimated frequencies: use the first frank_final indices.
    wk = w[idx[1:frank]]
    
    # Estimate weights via NNLS
    zk, err2 = nnls_weight(bcf, t, B)
    
    # Sort frequencies and corresponding coefficients in ascending order
    perm = sortperm(wk)
    wk = wk[perm]
    zk = zk[perm]
    
    Nsp = frank
    # Remove any zero coefficients.
    if minimum(zk) == 0.0
        keep = findall(x -> x > 0.0, zk)
        wk = wk[keep]
        zk = zk[keep]
        Nsp = length(wk)
    end

    # Coefficients
    gk = zk .* sbeta.(wk; scale=icm2ifs)

    # Rescale the frequencies and coefficients
    wk = wk ./ icm2ifs
    gk = gk ./ icm2ifs^2.0
    
    println("Number of sample points: ", Nsp)
    println("Error in ID: ", err1)
    println("Error in NNLS: ", err2)
    
    return (nsp=Nsp, freq=wk, coef=gk)
end
