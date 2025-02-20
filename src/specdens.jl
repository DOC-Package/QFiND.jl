
# Define Abstract Type for Spectral Density Models
abstract type SpectralDensity <: Function end

# Power-law with exponential cutoff
struct PowerLawExpSD <: SpectralDensity
    s    :: Real   # exponent
    α    :: Real   # coefficient 
    γ :: Real   # cutoff frequency
end

# Tannor-Meyer
struct TannorMeyerSD <: SpectralDensity
    Ω :: Vector{Float64}  # vector of central frequencies
    Γ :: Vector{Float64}  # vector of widths
    λ :: Vector{Float64}  # vector of intensities
end

# Brownian oscillator
struct BrownianSD <: SpectralDensity
    Ω :: Vector{Float64}  # vector of local frequencies
    Γ :: Vector{Float64}  # vector of damping coefficients
    λ :: Vector{Float64}  # vector of coupling strengths
end

struct AAAfittedSD <: SpectralDensity
    bary::Barycentric
end

"""
    (specdens::PowerLawExpSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Power-law with exponential cutoff model.
"""
function (specdens::PowerLawExpSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    s = specdens.s
    α = specdens.α
    γ = specdens.γ * scale
    return sgn * π * α * γ^(1.0 - s) * ω^s * exp(-ω / γ)
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
        p    = 4.0 * specdens.Γ[i] * specdens.λ[i] * (specdens.Ω[i]^2 + specdens.Γ[i]^2)
        deno = ((ω + specdens.Ω[i])^2 + specdens.Γ[i]^2) * ((ω - specdens.Ω[i])^2 + specdens.Γ[i]^2)
        res += p * ω / deno
    end
    return sgn*res
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
    for i in eachindex(specdens.Ω)
        p    = 2.0 * specdens.Γ[i] * specdens.λ[i] * specdens.Ω[i]^2
        deno = (ω^2 - specdens.Ω[i]^2)^2 + (specdens.Γ[i]^2 * ω^2)
        res += p * ω / deno
    end
    return sgn*res
end

"""
    (specdens::AAAfittedSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Brownian oscillator model.
"""
function (specdens::AAAfittedSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω) / scale
    res = evaluate(specdens.bary, ω)
    return sgn * scale * res
end


# Define Abstract Type for the QNSD
abstract type QuantumNoiseSpectralDensity <: Function end

# Bosonic QNSD
struct BosonicQNSD <: QuantumNoiseSpectralDensity
    specdens :: SpectralDensity  
    Temp :: Float64          
end

"""
    

"""

function (b::BosonicQNSD)(ω::Float64; scale::Float64=1.0) :: Float64
    β = ħ * 1e15 / (kb * b.Temp)
    if b.Temp == 0.0
        return (b.specdens)(ω; scale=scale) / π
    else
        factor = (1.0 / tanh(0.5 * β * icm2ifs / scale * ω) + 1.0) / (2π)
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
    β = ħ * 1e15 / (kb * f.Temp)
    if b.Temp == 0.0
        return (b.specdens)(ω; scale=scale) / π
    else
        factor = (1.0 / tanh(0.5 * β * icm2ifs / scale * ω) + 1.0) / (2π)
        return (b.specdens)(ω; scale=scale) * factor
    end
end





