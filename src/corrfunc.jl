abstract type BathCorrelationFunction <: Function end
abstract type BosonicBathCorrelationFunction <: BathCorrelationFunction end
struct BosonicBCF <: BosonicBathCorrelationFunction 
    specdens::Function
    Temp::Real
    Ω_c::Real
end
struct BosonicBCF_Real <: BosonicBathCorrelationFunction 
    specdens::Function
    Temp::Real
    Ω_c::Real
end
struct BosonicBCF_Imag <: BosonicBathCorrelationFunction 
    specdens::Function
    Temp::Real
    Ω_c::Real
end


function make_bcf(specdens::Function, Temp::Real, Ω_c::Real)
    return t -> begin
        return bcf_re(specdens, Temp, Ω_c, t) + 1.0im * bcf_im(specdens, Ω_c, t)
    end
end

# Real part of a BCF for a time t.
function bcf_re(specdens::Function, Temp::Real, Ω_c::Real, t::Real)
    if Temp == 0.0
        integrand_0K(ω) = specdens(ω; scale=icm2ifs) * cos(ω * t) / π
        (res, err) = quadgk(ω -> integrand_0K(ω), 0.0, Ω_c * icm2ifs; atol=1e-12)
    else
        β = ħ * 1e15 / (kb * Temp)
        integrand(ω) = (specdens(ω; scale=icm2ifs) / tanh(0.5 * β * ω)) * cos(ω*t) / π
        (res, err) = quadgk(ω -> integrand(ω), 0.0, Ω_c * icm2ifs; atol=1e-12)
    end
    return res
end

# Imaginary part of a BCF for a time t.
function bcf_im(specdens::Function, Ω_c::Real, t::Real)
    integrand(ω) = - specdens(ω; scale=icm2ifs) * sin(ω * t) / π
    (res, err) = quadgk(ω -> integrand(ω), 0.0, Ω_c * icm2ifs; atol=1e-12)
    return res
end
