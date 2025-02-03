
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


"""
    (sd::PowerLawExpSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Power-law with exponential cutoff model.
"""
function (sd::PowerLawExpSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    s = sd.s
    α = sd.α
    γ = sd.γ * scale
    return sgn * π * α * γ^(1.0 - s) * ω^s * exp(-ω / γ)
end


"""
    (sd::TannorMeyerSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Tannor-Meyer model.
"""
function (sd::TannorMeyerSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    Ω = sd.Ω .* scale
    Γ = sd.Γ .* scale
    λ = sd.λ .* scale
    res = 0.0
    for i in eachindex(sd.Ω)
        p    = 4.0 * sd.Γ[i] * sd.λ[i] * (sd.Ω[i]^2 + sd.Γ[i]^2)
        deno = ((ω + sd.Ω[i])^2 + sd.Γ[i]^2) * ((ω - sd.Ω[i])^2 + sd.Γ[i]^2)
        res += p * ω / deno
    end
    return sgn*res
end


"""
    (sd::BrownianSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Brownian oscillator model.
"""
function (sd::BrownianSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    Ω = sd.Ω .* scale
    Γ = sd.Γ .* scale
    λ = sd.λ .* scale
    res = 0.0
    for i in eachindex(sd.Ω)
        p    = 2.0 * sd.Γ[i] * sd.λ[i] * sd.Ω[i]^2
        deno = (ω^2 - sd.Ω[i]^2)^2 + (sd.Γ[i]^2 * ω^2)
        res += p * ω / deno
    end
    return sgn*res
end


# Define Abstract Type for the QNSD
abstract type QuantumNoiseSpectralDensity <: Function end

# Bosonic QNSD
struct BosonicQNSD <: QuantumNoiseSpectralDensity
    sd :: SpectralDensity  
    Temp :: Float64          
end

"""
    

"""

function (b::BosonicQNSD)(ω::Float64; scale::Float64=1.0)
    β = ħ * 1e15 / (kb * b.Temp)
    if b.Temp == 0.0
        return (b.sd)(ω; scale=scale) / π
    else
        factor = (1.0 / tanh(0.5 * β * icm2ifs / scale * ω) + 1.0) / (2π)
        return (b.sd)(ω; scale=scale) * factor
    end
end



# Fermionic QNSD
abstract type FermionicQNSD <: QuantumNoiseSpectralDensity end
struct FermionicQNSD_Plus <: FermionicQNSD 
    sd :: SpectralDensity
    Temp :: Float64      
end
struct FermionicQNSD_Minus <: FermionicQNSD
    sd :: SpectralDensity
    Temp :: Float64      
end






