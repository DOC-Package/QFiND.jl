function padeN_Nm1(nlt::Int)
    δ(x,y) = ==(x,y)
    bn = n -> 2*n + 1

    # Build LAM (size: (2*nlt)×(2*nlt))
    LAM = zeros(Float64, 2*nlt, 2*nlt)
    for n in 1:(2*nlt)
        for m in 1:(2*nlt)
            LAM[m, n] = (δ(m,n+1) + δ(m,n-1)) / sqrt(bn(m) * bn(n))
        end
    end

    # Compute eigenvalues of LAM.
    evals = eigvals(LAM)
    evals_sorted = sort(evals)
    xik = [2.0 / real(evals_sorted[2*nlt - n + 1]) for n in 1:nlt]

    # Build LAMt (size: (2*nlt-1)×(2*nlt-1))
    LAMt = zeros(Float64, 2*nlt-1, 2*nlt-1)
    for n in 1:(2*nlt-1)
        for m in 1:(2*nlt-1)
            LAMt[m, n] = (δ(m,n+1) + δ(m,n-1)) / sqrt(bn(m+1) * bn(n+1))
        end
    end

    evalt = eigvals(LAMt)
    evalt_sorted = sort(evalt)
    zeta_vec = [2.0 / real(evalt_sorted[2*nlt - n]) for n in 1:(nlt-1)]

    # Compute etak:
    etak = zeros(Float64, nlt)
    for n in 1:nlt
        nume = 1.0
        deno = 1.0
        for j in 1:(nlt-1)
            nume *= (zeta_vec[j]^2 - xik[n]^2)
        end
        for j in 1:nlt
            if j != n
                deno *= (xik[j]^2 - xik[n]^2)
            end
        end
        etak[n] = 0.5 * nlt * bn(nlt+1) * nume / deno
    end

    return xik, etak
end

function padeN_N(nlt::Int)
    δ(x,y) = ==(x,y)
    bn = n -> 2*n + 1
    RN = 1.0 / (4.0 * (nlt + 1) * bn[nlt+1])
    
    # Allocate LAM of size (2*nlt+1)×(2*nlt+1)
    LAM = zeros(2*nlt+1, 2*nlt+1)
    for n in 1:(2*nlt+1)
        for m in 1:(2*nlt+1)
            LAM[m, n] = (δ(m,n+1) + δ(m,n-1)) / sqrt(bn(m) * bn(n))
        end
    end
    evals = eigvals(LAM)
    evals_sorted = sort(evals)  # ascending order
    xik = [2.0 / real(evals_sorted[2*nlt+2 - n]) for n in 1:nlt ]
    
    # Allocate LAMt of size (2*nlt-1)×(2*nlt-1)
    LAMt = zeros(2*nlt-1, 2*nlt-1)
    for n in 1:(2*nlt-1)
        for m in 1:(2*nlt-1)
            LAMt[m, n] = (δ(m,n+1) + δ(m,n-1)) / sqrt(bn(m+1) * bn(n+1))
        end
    end
    evalt = eigen(LAMt).values
    evalt_sorted = sort(evalt)
    zeta = [ 2.0 / evalt_sorted[2*nlt - n] for n in 1:(nlt-1) ]
    
    # Compute etak: For each n=1:nlt, do a product over j.
    etak = zeros(Float64, nlt)
    for n in 1:nlt
        nume = prod( (zeta[j]^2 - xik[n]^2) for j in 1:(nlt-1) )
        deno = prod( (xik[j]^2 - xik[n]^2) for j in 1:nlt if j != n )
        etak[n] = 0.5 * RN * nume / deno
    end

    return xik, etak, RN
end
