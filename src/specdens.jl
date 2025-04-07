# Define Abstract Type for Spectral Density Models
abstract type SpectralDensity <: Function end

# Power-law with exponential cutoff
struct PowerLawExpSD <: SpectralDensity
    s :: Float64   # exponent
    γ :: Float64   # cutoff frequency
    α :: Float64   # coefficient 
    reorgene :: Float64   # reorganization energy
end

function PowerLawExpSD(s::Float64, γ::Float64; alpha::Union{Float64,Nothing}=nothing, reorgene::Union{Float64,Nothing}=nothing)
    if (alpha === nothing) == (reorgene === nothing)
        throw(ArgumentError("Please provide either α or reorgene, but not both."))
    end

    if alpha !== nothing
        reorgene = gamma(s)
        return PowerLawExpSD(s, γ, alpha, reorgene)
    else
        alpha = reorgene / gamma(s)
        return PowerLawExpSD(s, γ, alpha, reorgene)
    end
end

# Tannor-Meyer
struct TannorMeyerSD <: SpectralDensity
    Ω :: Vector{Float64}  # vector of central frequencies
    Γ :: Vector{Float64}  # vector of widths
    λ :: Vector{Float64}  # vector of intensities
end
TannorMeyerSD(Ω::Real, Γ::Real, λ::Real) = TannorMeyerSD([Float64(Ω)], [Float64(Γ)], [Float64(λ)])

# Brownian oscillator
struct BrownianSD <: SpectralDensity
    Ω :: Vector{Float64}  # vector of local frequencies
    Γ :: Vector{Float64}  # vector of damping coefficients
    λ :: Vector{Float64}  # vector of coupling strengths
end
BrownianSD(Ω::Real, Γ::Real, λ::Real) = BrownianSD([Float64(Ω)], [Float64(Γ)], [Float64(λ)])

struct DrudeSD <: SpectralDensity
    γ :: Float64
    λ :: Float64
end

struct AAAfittedSD <: SpectralDensity
    bary::Barycentric
    reorgene::Float64
    scale_reorgene::Float64
end
AAAfittedSD(bary::Barycentric; lb::Real=0.0, ub::Real=Inf) = AAAfittedSD(bary, reorganization_energy(bary; lb=lb,ub=ub), 1.0)
function AAAfittedSD(bary::Barycentric, reorgene::Float64; lb::Real=0.0, ub::Real=Inf)
    rawene = reorganization_energy(bary; lb=lb,ub=ub)
    scale_reorgene = reorgene / rawene
    return AAAfittedSD(bary, reorgene, scale_reorgene)
end

struct RationalSD <: SpectralDensity
    poles::Vector{ComplexF64}
    residues::Vector{ComplexF64}
end

struct WideBandSD <: SpectralDensity
    Γ :: Float64
end

function reorganization_energy(sd::SpectralDensity; lb::Real=0.0, ub::Real=Inf) :: Float64
    integrand(ω) = sd(ω) / ω / 2.0
    E_r, err = quadgk(ω -> integrand(ω), lb, ub; atol=1e-12)
    return E_r
end

function reorganization_energy(bary::Barycentric; lb::Real=0.0,       ub::Real=Inf)
    integrand(ω) = (evaluate(bary, ω)) / ω / 2.0
    E_r, err = quadgk(ω -> integrand(ω), lb, ub; atol=1e-12)
    return E_r
end

function reorganization_energy(ω::Vector{Float64}, g::Vector{Float64})
    return sum(g.^2.0 ./ abs.(ω)) / 2.0
end

"""
    (specdens::PowerLawExpSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Power-law with exponential cutoff model.
"""
function (specdens::PowerLawExpSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    s = specdens.s
    γ = specdens.γ * scale
    α = specdens.α * scale
    return sgn * 2.0 * α * γ^(-s) * ω^s * exp(-ω / γ)
end


"""
    (specdens::TannorMeyerSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Tannor-Meyer model.
"""
function (specdens::TannorMeyerSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    Ω = specdens.Ω .* scale
    Γ = specdens.Γ .* scale
    λ = specdens.λ .* scale
    res = 0.0
    for i in eachindex(specdens.Ω)
        p    = 8.0 / π * specdens.Γ[i] * specdens.λ[i] * (specdens.Ω[i]^2 + specdens.Γ[i]^2)
        deno = ((ω + specdens.Ω[i])^2 + specdens.Γ[i]^2) * ((ω - specdens.Ω[i])^2 + specdens.Γ[i]^2)
        res += p * ω / deno
    end
    return sgn * res
end

function sd_poles(sd::TannorMeyerSD)
    poles = ComplexF64[]
    for (Ω, Γ, _) in zip(sd.Ω, sd.Γ, sd.λ)
        push!(poles, Ω + im*(Γ/2))
        push!(poles, Ω - im*(Γ/2))
    end
    return poles
end

function sd_residues(sd::TannorMeyerSD)
    res = ComplexF64[]
    for (Ω, Γ, λ) in zip(sd.Ω, sd.Γ, sd.λ)
        res1 = λ * (Γ/2)^2 / (im*Γ)
        res2 = -λ * (Γ/2)^2 / (im*Γ)
        push!(res, res1)
        push!(res, res2)
    end
    return res
end

"""
    (specdens::BrownianSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Brownian oscillator model.
"""
function (specdens::BrownianSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    Ω = specdens.Ω .* scale
    Γ = specdens.Γ .* scale
    λ = specdens.λ .* scale
    res = 0.0
    for i in eachindex(Ω)
        p    = 4.0 / π * Γ[i] * λ[i] * Ω[i]^2
        deno = (ω^2 - Ω[i]^2)^2 + (Γ[i]^2 * ω^2)
        res += p * ω / deno
    end
    return sgn * res
end



function sd_nodes(sd::BrownianSD)
    poles = ComplexF64[]
    for (Ω, Γ, _) in zip(sd.Ω, sd.Γ, sd.λ)
        delta = sqrt(Ω^2 - (Γ/2)^2)
        push!(poles, delta - im*(Γ/2))
        push!(poles, -delta - im*(Γ/2))
    end
    return poles
end

function sd_residues(sd::BrownianSD)
    res = ComplexF64[]
    for (Ω, Γ, λ) in zip(sd.Ω, sd.Γ, sd.λ)
        delta = sqrt(Ω^2 - (Γ/2)^2)
        res1 = - (λ * Ω^2) / (4im*delta)
        res2 =   (λ * Ω^2) / (4im*delta)
        push!(res, res1)
        push!(res, res2)
    end
    return res
end

"""
    (specdens::AAAfittedSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Brownian oscillator model.
"""
function (specdens::AAAfittedSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω) / scale
    res = evaluate(specdens.bary, ω)
    return sgn * scale * res * specdens.scale_reorgene
end


function (specdens::RationalSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    poles = specdens.poles .* scale
    residues = specdens.residues .* scale
    res = 0.0
    for i in eachindex(specdens.poles)
        res += residues[i] / (ω - poles[i])
    end
    res = real(res)
    if res < 0.0
        res = 0.0
    end
    return sgn * scale * res
end


"""
    (specdens::DrudeSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Drude model.
"""
function (specdens::DrudeSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    γ = specdens.γ .* scale
    λ = specdens.λ .* scale
    res = 2.0 * λ * γ * ω / (ω^2 + γ^2)
    return sgn * res
end


# Define Abstract Type for the QNSD
abstract type QuantumNoiseSpectralDensity <: Function end

# Bosonic QNSD
struct BosonicQNSD <: QuantumNoiseSpectralDensity
    specdens :: SpectralDensity  
    Temp :: Float64          
end

struct BosonicQNSD_HighT <: QuantumNoiseSpectralDensity
    specdens :: SpectralDensity  
    Temp :: Float64          
end


function f_BE(ω::Float64, Temp::Float64)
    β = ħ * 1e15 / (kb * Temp)
    return 1.0 / (exp(β * icm2ifs * ω) - 1.0)
end


function (b::BosonicQNSD)(ω::Float64; scale::Float64=1.0) :: Float64
    if b.Temp == 0.0
        return (b.specdens)(ω; scale=scale)
    else
        β = ħ * 1e15 / (kb * b.Temp)
        factor = (1.0 / tanh(0.5 * β * ω * icm2ifs / scale) + 1.0) / 2.0
        return (b.specdens)(ω; scale=scale) * factor
    end
end

function (b::BosonicQNSD_HighT)(ω::Float64; scale::Float64=1.0) :: Float64
    if b.Temp == 0.0
        return (b.specdens)(ω; scale=scale) 
    else
        β = ħ * 1e15 / (kb * b.Temp)
        factor = 2.0 / (β * ω * icm2ifs / scale) / 2.0
        return (b.specdens)(ω; scale=scale) * factor
    end
end

# Fermionic QNSD
abstract type FermionicQNSD <: QuantumNoiseSpectralDensity end
struct FermionicQNSD_Plus <: FermionicQNSD 
    specdens :: SpectralDensity
    Temp :: Float64      
    ChemPot :: Float64
end
struct FermionicQNSD_Minus <: FermionicQNSD
    specdens :: SpectralDensity
    Temp :: Float64    
    ChemPot :: Float64  
end

function (f::FermionicQNSD_Plus)(ω::Float64; scale::Float64=1.0) :: Float64
    if b.Temp == 0.0
        return (b.specdens)(ω; scale=scale) 
    else
        β = ħ * 1e15 / (kb * b.Temp)
        factor = (1.0 / tanh(0.5 * β * icm2ifs / scale * ω) + 1.0) / 2.0
        return (b.specdens)(ω; scale=scale) * factor
    end
end





