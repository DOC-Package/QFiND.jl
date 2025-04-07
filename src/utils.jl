# Utility functions for the project

function equispaced_grid(
    Ω_min::Real, 
    Ω_max::Real, 
    tc::Real; 
    n_freq::Integer=1000, 
    n_time::Integer=200, 
    scale::Real=1.0)

    t = collect(range(0, tc, length=n_time))
    ω = collect(range(Ω_min, Ω_max, length=n_freq))
    ω = ω .* scale
    return t, ω
end

function equispaced_grid(T_c::Real; n_time::Integer=200)
    t = collect(range(0, T_c, length=n_time))
    return t
end

function equispaced_grid(Ω_min::Real, Ω_max::Real; n_freq::Integer=1000, scale::Real=1.0)
    ω = collect(range(Ω_min, Ω_max, length=n_freq))
    ω = ω .* scale
    return ω
end

function equispaced_discr(sbeta::Function, Ω_min, Ω_max, M_sp)
    ωk = collect(range(Ω_min, Ω_max, length=M_sp))
    # if ωk has 0 as an element, remove it
    ωk = filter(x -> x != 0, ωk)
    gk = sbeta.(ωk)
    return (freq = ωk, coef = gk)
end

function gausslegendre_discr(sbeta::Function, Ω_min, Ω_max, M_sp)
    roots, weights = gausslegendre(M_sp)
    factor = (Ω_max - Ω_min) / 2
    ωk = factor .* roots .+ (Ω_max + Ω_min) / 2
    gk = weights .* factor .* sbeta.(ωk)
    return (freq = ωk, coef = gk)
end

function create_integrand(sbeta::Vector{<:Real}, t::Vector{<:Real}, ω::Vector{<:Real}) :: Matrix{Float64}
    N_t = length(t)
    N_ω = length(ω)
    f = zeros(Float64, 2*N_t, N_ω)
    # Fill first N_t rows (real part)
    for i in 1:N_t
        for j in 1:N_ω
            f[i, j] = sbeta[j] * cos(ω[j] * t[i]) / π
        end
    end
    # Fill next N_t rows (imaginary part)
    for i in (N_t+1):(2*N_t)
        for j in 1:N_ω
            f[i, j] = -sbeta[j] * sin(ω[j] * t[i - N_t]) / π
        end
    end
    return f
end

function create_integrand_c(sbeta::Vector{<:Real}, t::Vector{<:Real}, ω::Vector{<:Real}) :: Matrix{ComplexF64}
    N_t = length(t)
    N_ω = length(ω)
    f = zeros(ComplexF64, N_t, N_ω)
    for i in 1:N_t
        for j in 1:N_ω
            f[i, j] = sbeta[j] * exp(-1.0im * ω[j] * t[i]) / π
        end
    end
    return f
end

function linsys_weight(cc::AbstractVector{<:Complex{<:Real}}, idx::Vector{Int}, D::Matrix{ComplexF64})
    frank = size(D, 2)
    c = cc[idx[1:frank]]
    g = D \ c
    err = norm(D * g - c)
    return g, err
end

# Calculate the sum of the exponential function.
function bcf_approx(t::Float64, ω::AbstractVector{Float64}, g::AbstractVector{<:Number}) :: ComplexF64
    ω = ω .* icm2ifs
    g = g .* icm2ifs
    return sum(g.^2.0 .* exp.(-im .* ω .* t) ./ 2.0)
end

function save_freq_coeff(freq::Vector{Float64}, coeff::Vector{Float64}, filename::String) 
    open(filename, "w") do io
        write(io, "-"^41 * "\n") 
        write(io, " "^4 * "Frequencies" * " "^10 * "Coefficients\n")
        write(io, "-"^41 * "\n")  
        
        for (ω, g) in zip(freq, coeff)
            write(io, @sprintf("%0.12e    %0.12e\n", 
                               ω, g))
        end
    end
end

function lawson(x::Vector{<:Real}, f::Function, r::Barycentric, nsteps::Integer)
    x1 = setdiff(x, r.nodes)
    ⍺, β = lawson(x1, f.(x1), r.nodes, r.values, r.weights, nsteps)
    return Barycentric(r.nodes, ⍺ ./ β, β, ⍺)
end

function lawson(x::Vector{<:Real}, f::Vector{<:Real}, r::Barycentric, nsteps::Integer)
    idx = findall(xi -> !(xi in r.nodes), x)
    x1 = x[idx]
    f1 = f[idx]
    ⍺, β = lawson(x1, f1, r.nodes, r.values, r.weights, nsteps)
    return Barycentric(r.nodes, ⍺ ./ β, β, ⍺)
end