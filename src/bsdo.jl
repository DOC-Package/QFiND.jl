"""
    ab = lanczos(xw, N)

Lanczos algorithm for the computation of the first N recurrence

Given the discrete inner product whose nodes are contained 
in the first column, and whose weights are contained in the 
second column, of the (Ncap)x2 array xw, the call

generates the first N recurrence coefficients ab of the 
corresponding discrete orthogonal polynomials. The N alpha-
coefficients are stored in the first column, the N beta-
coefficients in the second column, of the Nx2 array ab.

The code is adapted from the routine RKPW in
W.B. Gragg and W.J. Harrod, "The numerically stable 
reconstruction of Jacobi matrices from spectral data", 
Numer. Math. 44 (1984), 317-335.
"""

function lanczos(xw::AbstractArray{Float64,2}, N::Int)
    Ncap = size(xw, 1)

    # Check range of N
    if N <= 0 || N > Ncap
        error("N out of range")
    end

    # Create p0, p1 arrays (size: Ncap+1)
    p0 = zeros(Float64, Ncap + 1)
    p1 = zeros(Float64, Ncap + 1)

    p0[2:end] .= xw[:, 1]
    p1[1] = xw[1, 2]

    # Main loop
    for n in 1:(Ncap - 1)
        pn   = xw[n + 1, 2]
        gam  = 1.0
        sig  = 0.0
        t    = 0.0
        xlam = xw[n + 1, 1]
        for k in 1:(n + 1)
            rho = p1[k] + pn
            tmp = gam * rho
            tsig = sig

            if rho <= 0
                gam = 1.0
                sig = 0.0
            else
                gam = p1[k] / rho
                sig = pn    / rho
            end

            tk = sig * (p0[k] - xlam) - gam * t
            p0[k] = p0[k] - (tk - t)
            t = tk

            if sig <= 0
                pn = tsig * p1[k]
            else
                pn = (t^2) / sig
            end
            tsig = sig
            p1[k] = tmp
        end
    end

    ab = hcat(p0[1:N], p1[1:N])

    return ab
end

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
    bsdo(sbeta, ω::AbstractVector{Float64}, M_sp::Int)

Compute the discretization of the QNSD using the BSDO method.
"""
function bsdo_discr(sbeta::Function, ω::AbstractVector{Float64}, M_sp::Int)

    # Compute the quantum noise spectral density
    S = sbeta.(ω)
    S = max.(S, 0.0)

    # Combine w and j into a 2D array
    ωS = hcat(ω, S)

    # Discretize
    ωk, gk = orthpoly_discretization(ωS, M_sp)

    # Normalize zk
    if Ω_min < 0
        norm1, err1 = quadgk(x -> sbeta(x), Ω_min, 0)
        norm2, err2 = quadgk(x -> sbeta(x), 0, Ω_max)
        norm = norm1 + norm2
    else
        norm, err = quadgk(x -> sbeta(x), Ω_min, Ω_max)
    end
    gk .*= norm

    return (freq = ωk, coef = gk)
end


function bsdo_discr(sbeta::Function, Ω_min::Real, Ω_max::Real, M_sp::Int; n_lanczos::Int=1000)

    # Generate the frequencies
    ω = range(Ω_min, Ω_max, length=n_lanczos) |> collect

    return bsdo_discr(sbeta, ω, M_sp)
end
