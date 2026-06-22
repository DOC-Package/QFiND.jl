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
function bcf_approx(t::Float64, ω::AbstractVector{Float64}, g::AbstractVector{Float64}) :: ComplexF64
    ω = ω .* icm2ifs
    g = g .* icm2ifs
    return sum(g.^2.0 .* exp.(-im .* ω .* t) ./ 2.0)
end

function bcf_approx(t::Float64, a::AbstractVector{<:Number}, c::AbstractVector{<:Number}) :: ComplexF64
    return sum(c .* exp.(-a .* t))
end

function bcf_approx(t::AbstractVector{<:Real}, a::AbstractVector{<:Number}, c::AbstractVector{<:Number}) :: Vector{ComplexF64}
    return [sum(c .* exp.(-a .* t_i)) for t_i in t]
end


function save_freq_coeff(freq::Vector{Float64}, coeff::Vector{Float64}, filename::String) 
    open(filename, "w") do io
        write(io, "="^46 * "\n") 
        write(io, " "^1 * "Frequencies [cm^-1]" * " "^5 * "Coefficients [cm^-1]\n")
        write(io, "="^46 * "\n")  
        
        for (ω, g) in zip(freq, coeff)
            write(io, @sprintf("%0.12e      %0.12e\n", 
                               ω, g))
        end
    end
end

function save_expon_coeff(expon::Vector{ComplexF64}, coeff::Vector{ComplexF64}, filename::String)
    open(filename, "w") do io
        write(io, "# coeff (real)    coeff (imag)    expon (real)    expon (imag)\n")
        for (c1, c2) in zip(coeff, expon)
            write(io, @sprintf("%0.12e  %0.12e  %0.12e  %0.12e\n", 
                               real(c1), imag(c1), real(c2), imag(c2)))
        end
    end
end

function save_expon_coeff_union(expon::Vector{ComplexF64}, coeff::Vector{ComplexF64}, filename::String)
    open(filename, "w") do io
        for (c1, c2) in zip(coeff, expon)
            write(io, @sprintf("%0.12e  %0.12e  %0.12e  %0.12e  %0.12e  %0.12e\n", 
                               real(c1), imag(c1), 0.0, 0.0, real(c2), imag(c2)))
            write(io, @sprintf("%0.12e  %0.12e  %0.12e  %0.12e  %0.12e  %0.12e\n", 
                               0.0, 0.0, real(c1), -imag(c1), real(c2), -imag(c2)))
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

function AAApoles(r::Barycentric{T,S}) where {T,S}
    w = weights(r)
    nonzero = @. !iszero(w)
    z, w = nodes(r)[nonzero], w[nonzero]
    m = length(w)
    B = diagm( [zero(T); ones(T, m)] )
    E = [zero(T) transpose(w); ones(T, m) diagm(z) ];
    pol = []  # put it into scope
    try
        pol = filter( isfinite, eigvals(E, B) )
    catch
        # generalized eigen not available in extended precision, so:
        λ = filter( z->abs(z)>1e-13, eigvals(E\B) )
        pol = 1 ./ λ
        (; α, β) = schur(complex(E), complex(B))
        pol = filter( isfinite, α ./ β )
    end
    return pol
end

function AAAresidues(r::Barycentric)
    ζ = AAApoles(r)
    res = similar( complex(ζ) )
    z, y, w = nodes(r), values(r), weights(r)
    for (i, t) in pairs(ζ)
        numer = sum( w*y / (t-z) for (z, y, w) in zip(z, y, w))
        denomdiff = -sum( w / (t-z)^2 for (z, w) in zip(z, w))
        res[i] = numer / denomdiff
    end
    return ζ, res
end

function AAArinf(r::Barycentric)
    r_inf = sum(r.weights .* r.values) / sum(r.weights)
    return r_inf
end

function AAAPartialFraction(r::Barycentric)
    pols, res = AAAresidues(r)
    rinf = AAArinf(r)
    return pols, res, rinf
end

"""
    split_poles_by_linewidth(poles, residues; threshold)

Split poles and residues into two groups based on the imaginary part (linewidth) of the poles.

# Arguments
- `poles::Vector{ComplexF64}`: Vector of poles
- `residues::Vector{ComplexF64}`: Vector of residues corresponding to poles
- `threshold::Real`: Threshold for the absolute value of the imaginary part

# Returns
- `narrow`: Named tuple (poles, residues) for poles with |Im(p)| < threshold (narrow linewidth)
- `broad`: Named tuple (poles, residues) for poles with |Im(p)| >= threshold (broad linewidth)
"""
function split_poles_by_linewidth(poles::Vector{ComplexF64}, residues::Vector{ComplexF64}; threshold::Real)
    narrow_idx = findall(p -> abs(imag(p)) < threshold, poles)
    broad_idx = findall(p -> abs(imag(p)) >= threshold, poles)
    
    narrow = (poles = poles[narrow_idx], residues = residues[narrow_idx])
    broad = (poles = poles[broad_idx], residues = residues[broad_idx])
    
    return narrow, broad
end

"""
    split_poles_by_linewidth(r::Barycentric; threshold)

Split poles and residues from a Barycentric rational approximation based on linewidth.
"""
function split_poles_by_linewidth(r::Barycentric; threshold::Real)
    pols, res = AAAresidues(r)
    return split_poles_by_linewidth(pols, res; threshold=threshold)
end

"""
    filter_poles_by_linewidth(poles, residues; min_width=0.0, max_width=Inf)

Filter poles and residues based on the imaginary part (linewidth) range.

# Arguments
- `poles::Vector{ComplexF64}`: Vector of poles
- `residues::Vector{ComplexF64}`: Vector of residues
- `min_width::Real`: Minimum linewidth (default: 0.0)
- `max_width::Real`: Maximum linewidth (default: Inf)

# Returns
- Named tuple (poles, residues) for poles with min_width <= |Im(p)| < max_width
"""
function filter_poles_by_linewidth(poles::Vector{ComplexF64}, residues::Vector{ComplexF64}; 
                                    min_width::Real=0.0, max_width::Real=Inf)
    idx = findall(p -> min_width <= abs(imag(p)) < max_width, poles)
    return (poles = poles[idx], residues = residues[idx])
end

"""
    sort_poles_by_linewidth(poles, residues; rev=false)

Sort poles and residues by the absolute value of the imaginary part (linewidth).

# Arguments
- `poles::Vector{ComplexF64}`: Vector of poles
- `residues::Vector{ComplexF64}`: Vector of residues
- `rev::Bool`: If true, sort in descending order (default: false, ascending)

# Returns
- Named tuple (poles, residues) sorted by linewidth
"""
function sort_poles_by_linewidth(poles::Vector{ComplexF64}, residues::Vector{ComplexF64}; rev::Bool=false)
    idx = sortperm(abs.(imag.(poles)); rev=rev)
    return (poles = poles[idx], residues = residues[idx])
end

"""
    poles_to_exponents(poles, residues)

Convert poles and residues of Sβ(ω) to exponential sum coefficients for correlation function.

Given the partial fraction decomposition:
    S̃_β(ω) = Σ_k η_k / (ω - z_k)

The correlation function is:
    C̃(t) = Σ_k c_k exp(-a_k t)

where z'_k are poles in the lower half plane (Im(z_k) < 0),
c_k = -i η_k, and a_k = i z'_k.

# Arguments
- `poles::Vector`: All poles z_k (can be real or complex)
- `residues::Vector`: Corresponding residues η_k

# Returns
- Named tuple (c, a) where C(t) = Σ_k c_k exp(-a_k t)
  - `c`: Exponential coefficients c_k = -i η_k
  - `a`: Exponents a_k = i z'_k (should have Re(a_k) > 0 for decay)
"""
function poles_to_exponents(poles::Vector, residues::Vector)
    poles_c = complex.(poles)
    residues_c = complex.(residues)
    
    # Select poles in the lower half plane (Im(z) < 0)
    lower_idx = findall(z -> imag(z) < 0, poles_c)
    z_lower = poles_c[lower_idx]
    η_lower = residues_c[lower_idx]
    
    # c_k = -i η_k, a_k = i z'_k
    c = -im .* η_lower
    a = im .* z_lower
    
    return (c = c, a = a)
end

"""
    poles_to_exponents(r::Barycentric)

Convert Barycentric rational approximation to exponential sum coefficients.
"""
function poles_to_exponents(r::Barycentric)
    pols, res = AAAresidues(r)
    return poles_to_exponents(pols, res)
end