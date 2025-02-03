"""
    orthpoly_discretization(wj, N)

Compute the discretization of orthogonal polynomials.

Parameters
----------
- wj::Matrix{Float64}
    An n×2 array where:
      - wj[:,1] is the set of support points of the spectral density (frequencies).
      - wj[:,2] is the spectral density values.
- N::Int
    The number of discretization points.

Returns
-------
- wd::Vector{Float64}
    The discretized frequencies.
- zd::Vector{Float64}
    The coefficients.
"""
function orthpoly_discretization(wj::Array{Float64,2}, N::Int)
    # Compute recurrence coefficients (Alpha,Beta) from Lanczos
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
    bsdo(N_w, Omega_min, Omega_max, Msp)

Compute the discretization of the spectral density using the BSDO method.

Parameters
----------
- N_w::Int
    The number of discretization points.
- Omega_min::Float64
    The minimum frequency.
- Omega_max::Float64
    The maximum frequency.
- Msp::Int
    Number of polynomials or states to consider in orthpoly discretization.

Returns
-------
A tuple `(wk, zk)` where:
- wk::Vector{Float64}
    The discretized frequencies.
- gk::Vector{Float64}
    The scaled/summed weights for those frequencies.
"""
function bsdo_discr(sbeta::Function, Omega_min::Real, Omega_max::Real, N_w::Int, Msp::Int)

    # Generate the frequencies
    w = range(Omega_min, Omega_max, length=N_w) |> collect

    # Compute the quantum noise spectral density
    j = sbeta.(w)
    j = max.(j, 0.0)

    # Combine w and j into a 2D array
    wj = hcat(w, j)

    # Discretize
    wk, gk = orthpoly_discretization(wj, Msp)

    # Normalize zk
    if Omega_min < 0
        norm1, err1 = quadgk(x -> sbeta(x), Omega_min, 0)
        norm2, err2 = quadgk(x -> sbeta(x), 0, Omega_max)
        norm = norm1 + norm2
    else
        norm, err = quadgk(x -> sbeta(x), Omega_min, Omega_max)
    end
    gk .*= norm

    return (freq = wk, coef = gk)
end
