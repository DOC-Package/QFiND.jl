function psd(sd::SpectralDensity, Temp::Real, npade::Int; pade_type::Symbol=:Nm1)
    β = ħ * 1e15 / (kb * Temp)
    expon = ComplexF64[]
    coeff = ComplexF64[]

    # Contributions from the spectral density
    poles_J = sd_poles(sd; scale=icm2ifs)
    res_J = sd_residues(sd; scale=icm2ifs)
    for (ω_p, r_p) in zip(poles_J, res_J)
        if imag(ω_p) < 0
            A = r_p * (coth((β*ω_p)/2) + 1)
            push!(expon, 1.0im * ω_p)
            push!(coeff,    (-1.0im) * A)
        end
    end

    # Contributions from the Pade approximation
    if npade > 0
        if pade_type == :Nm1
            ξ, η = padeN_Nm1(npade)
        elseif pade_type == :N
            ξ, η, _ = padeN_N(npade)
        else
            throw(ArgumentError("pade_type must be :N or :Nm1"))
        end
        
        poles_pade = ξ ./ β
        for (ω_p, η_j) in zip(poles_pade, η)
            A = (2im * η_j / β) * sd(1im * ω_p; scale=icm2ifs)
            push!(expon, ω_p)
            push!(coeff, A)
        end
    end

    return expon, coeff
end

function tpsd(sd::SpectralDensity, Temp::Real, npade::Int, tol::Real; pade_type::Symbol=:Nm1)
    β = ħ * 1e15 / (kb * Temp)
    expon = ComplexF64[]
    coeff = ComplexF64[]

    # Contributions from the spectral density
    poles_J = sd_poles(sd; scale=icm2ifs)
    res_J = sd_residues(sd; scale=icm2ifs)
    for (ω_p, r_p) in zip(poles_J, res_J)
        if imag(ω_p) < 0
            A = r_p * (coth((β*ω_p)/2) + 1)
            push!(expon, 1.0im * ω_p)
            push!(coeff, (-1.0im) * A)
        end
    end

    # Contributions from the Pade approximation
    expon_psd = ComplexF64[]
    coeff_psd = ComplexF64[]
    if npade > 0
        if pade_type == :Nm1
            ξ, η = padeN_Nm1(npade)
        elseif pade_type == :N
            ξ, η, _ = padeN_N(npade)
        else
            throw(ArgumentError("pade_type must be :N or :Nm1"))
        end
        
        poles_pade = ξ ./ β
        for (ω_p, η_j) in zip(poles_pade, η)
            A = (2im * η_j / β) * sd(1im * ω_p; scale=icm2ifs)
            push!(expon_psd, ω_p)
            push!(coeff_psd, A)
        end
    end
    expon_tpsd, coeff_tpsd = QFiND.balanced_truncation(expon_psd, coeff_psd, tol)
    push!(expon, expon_tpsd...)
    push!(coeff, coeff_tpsd...)

    return expon, coeff
end

function tpsd(sd::SpectralDensity, Temp::Real, npade::Int, ntrun::Int; pade_type::Symbol=:Nm1)
    # ntpsd must be smaller than npade
    if ntrun > npade
        throw(ArgumentError("ntpsd must be smaller than or equal to npade."))
    end
    β = ħ * 1e15 / (kb * Temp)
    expon = ComplexF64[]
    coeff = ComplexF64[]

    # Contributions from the spectral density
    poles_J = sd_poles(sd; scale=icm2ifs)
    res_J = sd_residues(sd; scale=icm2ifs)
    for (ω_p, r_p) in zip(poles_J, res_J)
        if imag(ω_p) < 0
            A = r_p * (coth((β*ω_p)/2) + 1)
            push!(expon, 1.0im * ω_p)
            push!(coeff, (-1.0im) * A)
        end
    end

    # Contributions from the Pade approximation
    expon_psd = ComplexF64[]
    coeff_psd = ComplexF64[]
    if npade > 0
        if pade_type == :Nm1
            ξ, η = padeN_Nm1(npade)
        elseif pade_type == :N
            ξ, η, _ = padeN_N(npade)
        else
            throw(ArgumentError("pade_type must be :N or :Nm1"))
        end

        poles_pade = ξ ./ β
        for (ω_p, η_j) in zip(poles_pade, η)
            A = (2im * η_j / β) * sd(1im * ω_p; scale=icm2ifs)
            push!(expon_psd, ω_p)
            push!(coeff_psd, A)
        end
    end
    expon_tpsd, coeff_tpsd = QFiND.balanced_truncation(expon_psd, coeff_psd, ntrun)
    push!(expon, expon_tpsd...)
    push!(coeff, coeff_tpsd...)

    return expon, coeff
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
    zeta = [2.0 / real(evalt_sorted[2*nlt - n]) for n in 1:(nlt-1)]

    eta = zeros(Float64, nlt)
    for n in 1:nlt
        nume = 1.0
        deno = 1.0
        for j in 1:(nlt-1)
            nume *= (zeta[j]^2 - xi[n]^2)
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
    RN = 1.0 / (4.0 * (nlt + 1) * bn(nlt+1))
    
    LAM = zeros(2*nlt+1, 2*nlt+1)
    for n in 1:(2*nlt+1)
        for m in 1:(2*nlt+1)
            LAM[m, n] = (δ(m,n+1) + δ(m,n-1)) / sqrt(bn(m) * bn(n))
        end
    end
    evals = eigvals(LAM)
    evals_sorted = sort(evals) 
    xi = [2.0 / real(evals_sorted[2*nlt+2 - n]) for n in 1:nlt ]
    
    LAMt = zeros(2*nlt, 2*nlt)
    for n in 1:2*nlt
        for m in 1:2*nlt
            LAMt[m, n] = (δ(m,n+1) + δ(m,n-1)) / sqrt(bn(m+1) * bn(n+1))
        end
    end
    evalt = eigen(LAMt).values
    evalt_sorted = sort(evalt)
    zeta = [ 2.0 / evalt_sorted[2*nlt - n + 1] for n in 1:nlt ]
    
    eta = zeros(Float64, nlt)
    for n in 1:nlt
        nume = 1.0
        deno = 1.0
        for j in 1:nlt
            nume *= (zeta[j]^2 - xi[n]^2)
        end
        for j in 1:nlt
            if j != n
                deno *= (xi[j]^2 - xi[n]^2)
            end
        end
        eta[n] = 0.5 * RN * nume / deno
        
    end

    return xi, eta, RN
end