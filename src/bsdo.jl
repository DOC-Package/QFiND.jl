"""
    orthpoly_discretization(wj, N)

Compute the discretization of orthogonal polynomials.
"""
function orthpoly_discretization(wj::Array{Float64,2}, N::Int)
    # Compute recurrence coefficients from Lanczos
    ab = lanczos(wj, N) 

    # Construct tridiagonal matrix M
    diag_main = diagm(0 => ab[:, 1])
    off_diag_elements = sqrt.(ab[2:end, 2])
    diag_upper = diagm(1 => off_diag_elements)
    diag_lower = diagm(-1 => off_diag_elements)
    M = diag_main + diag_upper + diag_lower

    # Compute the eigenvalues
    eig = eigen(Symmetric(M))
    wd = eig.values
    zd = eig.vectors[1, :].^2

    return wd, zd
end

"""
    bsdo(sbeta, Ω_min::Real, Ω_max::Real, N_ω::Int, M_sp::Int)

Compute the discretization of the QNSD using the BSDO method.
"""
function bsdo_discr(sbeta::Function, Ω_min::Real, Ω_max::Real, N_ω::Int, M_sp::Int)

    # Generate the frequencies
    ω = range(Ω_min, Ω_max, length=N_ω) |> collect

    # Compute the quantum noise spectral density
    S = sbeta.(w)
    S = max.(S, 0.0)

    # Combine w and j into a 2D array
    ωS = hcat(w, S)

    # Discretize
    freq, coef = orthpoly_discretization(ωS, M_sp)

    # Normalize zk
    if Omega_min < 0
        norm1, err1 = quadgk(x -> sbeta(x), Ω_min, 0)
        norm2, err2 = quadgk(x -> sbeta(x), 0, Ω_max)
        norm = norm1 + norm2
    else
        norm, err = quadgk(x -> sbeta(x), Ω_min, Ω_max)
    end
    coef .*= norm

    return (freq = freq, coef = coef)
end
