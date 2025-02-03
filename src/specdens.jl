export PowerLawExpSD, TannorMeyerSD, BrownianSD, make_sdens, make_sbeta

# Define Abstract Type for Spectral Density Models
abstract type SpectralDensity end

# Power-law with exponential cutoff
struct PowerLawExpSD <: SpectralDensity
    s    :: Real   # exponent
    α    :: Real   # coefficient 
    gamc :: Real   # cutoff frequency
end

# Tannor-Meyer
struct TannorMeyerSD <: SpectralDensity
    Omg :: Vector{Float64}  # vector of central frequencies
    Gam :: Vector{Float64}  # vector of widths
    Lam :: Vector{Float64}  # vector of intensities
end

# Brownian oscillator
struct BrownianSD <: SpectralDensity
    Omg :: Vector{Float64}  # vector of local frequencies
    Gam :: Vector{Float64}  # vector of damping coefficients
    Lam :: Vector{Float64}  # vector of coupling strengths
end


"""
    evaluate(sd::PowerLawExpSD, omega::Float64) -> Float64

Compute the spectral density for the Power-law with exponential cutoff model.
"""
function compute_spectral_density(sd::PowerLawExpSD, omega::Float64; scale::Float64=1.0)::Float64
    s = sd.s
    α = sd.α
    gamc = sd.gamc * scale
    return π * α * gamc^(1.0 - s) * omega^s * exp(-omega / gamc)
end

"""
    compute_density(sd::TannorMeyerSD, omega::Float64) -> Float64

Compute the spectral density for the Tannor-Meyer model.
"""
function compute_spectral_density(sd::TannorMeyerSD, omega::Float64; scale::Float64=1.0)::Float64
    Omg = sd.Omg .* scale
    Gam = sd.Gam .* scale
    Lam = sd.Lam .* scale
    res = 0.0
    for i in eachindex(sd.Omg)
        p    = 4.0 * sd.Gam[i] * sd.Lam[i] * (sd.Omg[i]^2 + sd.Gam[i]^2)
        deno = ((omega + sd.Omg[i])^2 + sd.Gam[i]^2) * ((omega - sd.Omg[i])^2 + sd.Gam[i]^2)
        res += p * omega / deno
    end
    return res
end

"""
    compute_density(sd::BrownianSD, omega::Float64) -> Float64

Compute the spectral density for the Brownian oscillator model.
"""
function compute_spectral_density(sd::BrownianSD, omega::Float64; scale::Float64=1.0)::Float64
    Omg = sd.Omg .* scale
    Gam = sd.Gam .* scale
    Lam = sd.Lam .* scale
    res = 0.0
    for i in eachindex(sd.Omg)
        p    = 2.0 * sd.Gam[i] * sd.Lam[i] * sd.Omg[i]^2
        deno = (omega^2 - sd.Omg[i]^2)^2 + (sd.Gam[i]^2 * omega^2)
        res += p * omega / deno
    end
    return res
end

"""
    sdens(sd::SpectralDensity, omega::Float64; scale::Float64=1.0) -> Float64

"""
sdens(sd::SpectralDensity, omega::Float64; scale::Float64=1.0) = compute_spectral_density(sd, omega; scale = scale)

"""
    make_sdens(sd::SpectralDensity; scale::Float64=1.0) -> Function
"""
function make_sdens(sd::SpectralDensity; scale::Float64=1.0)
    return (omega::Float64) -> sign(omega) * sdens(sd, abs(omega); scale=scale)
end

"""
    make_sbeta(sd::SpectralDensity, Temp::Float64, icm2ifs::Real, hbar::Real, kb::Real; scale::Float64=1.0)
        -> Function
"""
function make_sbeta(sd::SpectralDensity, Temp::Float64; scale::Float64=1.0)
    beta = hbar * 1e15 / (kb * Temp)
    return (omega::Float64) -> begin
        if Temp == 0.0
            return sign(omega) * sdens(sd, abs(omega); scale=scale) / π
        else
            factor = (1.0 / tanh(0.5 * beta * icm2ifs / scale * omega) + 1.0) / (2π)
            return sign(omega) * sdens(sd, abs(omega); scale=scale) * factor
        end
    end
end
