function spectral_decomposition(sd::SpectralDensity, Temp::Real; npade::Int=1)
    β = ħ * 1e15 / (kb * Temp)
    exponents = ComplexF64[]
    coeffs    = ComplexF64[]

    # Contributions from the spectral density
    poles_J  = sd_poles(sd; scale=icm2ifs)
    res_J    = sd_residues(sd; scale=icm2ifs)
    for (ω_p, r_p) in zip(poles_J, res_J)
        if imag(ω_p) < 0
            A = r_p * (coth((β*ω_p)/2) + 1)
            push!(exponents, 1.0im * ω_p)
            push!(coeffs,    (-1.0im) * A)
        end
    end

    # Contributions from the Pade approximation
    if npade > 0
        ξ, η = padeN_Nm1(npade)
        poles_pade = ξ ./ β
        for (ω_p, η_j) in zip(poles_pade, η)
            A = (2im * η_j / β) * sd(1im * ω_p; scale=icm2ifs)
            push!(exponents, ω_p)
            push!(coeffs,    A)
        end
    end

    return exponents, coeffs
end

function padeN_Nm1(nlt::Int)
    δ(x,y) = ==(x,y)
    bn = n -> 2*n + 1

    LAM = zeros(Float64, 2*nlt, 2*nlt)
    for n in 1:(2*nlt)
        for m in 1:(2*nlt)
            LAM[m, n] = (δ(m,n+1) + δ(m,n-1)) / sqrt(bn(m) * bn(n))
        end
    end
    evals = eigvals(LAM)
    evals_sorted = sort(evals)
    xi = [2.0 / real(evals_sorted[2*nlt - n + 1]) for n in 1:nlt]

    LAMt = zeros(Float64, 2*nlt-1, 2*nlt-1)
    for n in 1:(2*nlt-1)
        for m in 1:(2*nlt-1)
            LAMt[m, n] = (δ(m,n+1) + δ(m,n-1)) / sqrt(bn(m+1) * bn(n+1))
        end
    end
    evalt = eigvals(LAMt)
    evalt_sorted = sort(evalt)
    zeta_vec = [2.0 / real(evalt_sorted[2*nlt - n]) for n in 1:(nlt-1)]

    eta = zeros(Float64, nlt)
    for n in 1:nlt
        nume = 1.0
        deno = 1.0
        for j in 1:(nlt-1)
            nume *= (zeta_vec[j]^2 - xi[n]^2)
        end
        for j in 1:nlt
            if j != n
                deno *= (xi[j]^2 - xi[n]^2)
            end
        end
        eta[n] = 0.5 * nlt * bn(nlt+1) * nume / deno
    end

    return xi, eta
end

function padeN_N(nlt::Int)
    δ(x,y) = ==(x,y)
    bn = n -> 2*n + 1
    RN = 1.0 / (4.0 * (nlt + 1) * bn[nlt+1])
    
    LAM = zeros(2*nlt+1, 2*nlt+1)
    for n in 1:(2*nlt+1)
        for m in 1:(2*nlt+1)
            LAM[m, n] = (δ(m,n+1) + δ(m,n-1)) / sqrt(bn(m) * bn(n))
        end
    end
    evals = eigvals(LAM)
    evals_sorted = sort(evals) 
    xi = [2.0 / real(evals_sorted[2*nlt+2 - n]) for n in 1:nlt ]
    
    LAMt = zeros(2*nlt-1, 2*nlt-1)
    for n in 1:(2*nlt-1)
        for m in 1:(2*nlt-1)
            LAMt[m, n] = (δ(m,n+1) + δ(m,n-1)) / sqrt(bn(m+1) * bn(n+1))
        end
    end
    evalt = eigen(LAMt).values
    evalt_sorted = sort(evalt)
    zeta = [ 2.0 / evalt_sorted[2*nlt - n] for n in 1:(nlt-1) ]
    
    eta = zeros(Float64, nlt)
    for n in 1:nlt
        nume = prod( (zeta[j]^2 - xi[n]^2) for j in 1:(nlt-1) )
        deno = prod( (xi[j]^2 - xi[n]^2) for j in 1:nlt if j != n )
        eta[n] = 0.5 * RN * nume / deno
    end

    return xi, eta, RN
end