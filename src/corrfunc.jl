abstract type BathCorrelationFunction <: Function end

# 
abstract type BosonicBathCorrelationFunction <: BathCorrelationFunction end

# Bosonic BCF with real and imaginary parts
struct BosonicBCF <: BosonicBathCorrelationFunction 
    specdens::SpectralDensity
    Temp::Real
    lb::Real
    ub::Real
    rtol::Union{Real,Nothing}
    atol::Union{Real,Nothing}
end
BosonicBCF(specdens::SpectralDensity, Temp::Real; ub::Real=Inf, rtol::Union{Real,Nothing}=nothing, atol::Union{Real,Nothing}=nothing) = BosonicBCF(specdens, Temp, 0.0, ub, rtol, atol)


# Bosonic BCF with real part only
struct BosonicBCF_Real <: BosonicBathCorrelationFunction 
    specdens::SpectralDensity
    Temp::Real
    lb::Real
    ub::Real
    rtol::Union{Real,Nothing} 
    atol::Union{Real,Nothing}
end
BosonicBCF_Real(specdens::SpectralDensity, Temp::Real; ub::Real=Inf, rtol::Union{Real,Nothing}=nothing, atol::Union{Real,Nothing}=nothing) = BosonicBCF_Real(specdens, Temp, 0.0, ub, rtol, atol)

# Bosonic BCF with imaginary part only
struct BosonicBCF_Imag <: BosonicBathCorrelationFunction 
    specdens::SpectralDensity
    Temp::Real
    lb::Real
    ub::Real
    rtol::Union{Real,Nothing}
    atol::Union{Real,Nothing}
end
BosonicBCF_Imag(specdens::SpectralDensity, Temp::Real; ub::Real=Inf, rtol::Union{Real,Nothing}=nothing, atol::Union{Real,Nothing}=nothing) = BosonicBCF_Imag(specdens, Temp, 0.0, ub, rtol, atol)

struct BosonicBCF_dt <: BosonicBathCorrelationFunction 
    specdens::SpectralDensity
    Temp::Real
    lb::Real
    ub::Real
    rtol::Union{Real,Nothing}
    atol::Union{Real,Nothing}
end
BosonicBCF_dt(specdens::SpectralDensity, Temp::Real; ub::Real=Inf, rtol::Union{Real,Nothing}=nothing, atol::Union{Real,Nothing}=nothing) = BosonicBCF_dt(specdens, Temp, 0.0, ub, rtol, atol)


# Bosonic BCF
function _compute_bcf(b::BosonicBCF, t::Real)
    # Real part
    integrand_re = (ω::Float64) -> begin
        if b.Temp == 0.0
            b.specdens(ω; scale=icm2ifs) * cos(ω * t) / π
        else
            β = ħ * 1e15 / (kb * b.Temp)
            (b.specdens(ω; scale=icm2ifs) / tanh(0.5 * β * ω)) * cos(ω * t) / π
        end
    end
    S, err1 = quadgk(integrand_re, b.lb * icm2ifs, b.ub * icm2ifs; rtol=b.rtol, atol=b.atol)

    # Imag part
    integrand_im = (ω::Float64) -> -b.specdens(ω; scale=icm2ifs) * sin(ω * t) / π
    A, err2 = quadgk(integrand_im, b.lb * icm2ifs, b.ub * icm2ifs; rtol=b.rtol, atol=b.atol)
    return S + 1.0im * A, max(err1, err2)
end

function (b::BosonicBCF)(t::Real; ret_err::Bool=false)
    val, err = _compute_bcf(b, t)
    ret_err ? (val, err) : val
end

# Real part of a BCF
function _compute_bcf_real(b::BosonicBCF_Real, t::Real)
    integrand_re = (ω::Float64) -> begin
        if b.Temp == 0.0
            b.specdens(ω; scale=icm2ifs) * cos(ω * t) / π
        else
            β = ħ * 1e15 / (kb * b.Temp)
            (b.specdens(ω; scale=icm2ifs) / tanh(0.5 * β * ω)) * cos(ω * t) / π
        end
    end
    S, err = quadgk(integrand_re, b.lb * icm2ifs, b.ub * icm2ifs; rtol=b.rtol, atol=b.atol)
    return S, err
end

function (b::BosonicBCF_Real)(t::Real; ret_err::Bool=false)
    S, err = _compute_bcf_real(b, t)
    ret_err ? (S, err) : S
end

# Imaginary part of a BCF for a time t.
function _compute_bcf_imag(b::BosonicBCF_Imag, t::Real)
    integrand_im = (ω::Float64) -> -b.specdens(ω; scale=icm2ifs) * sin(ω * t) / π
    A, err = quadgk(integrand_im, b.lb * icm2ifs, b.ub * icm2ifs; rtol=b.rtol, atol=b.atol)
    return A, err
end

function (b::BosonicBCF_Imag)(t::Real; ret_err::Bool=false)
    A, err = _compute_bcf_imag(b, t)
    ret_err ? (A, err) : A
end

# Time derivative of Bosonic BCF
function _compute_bcf_dt(b::BosonicBCF_dt, t::Real)
    # Real part 
    integrand_re = (ω::Float64) -> begin
        if b.Temp == 0.0
            - ω * b.specdens(ω; scale=icm2ifs) * sin(ω * t) / π
        else
            β = ħ * 1e15 / (kb * b.Temp)
            - ω * (b.specdens(ω; scale=icm2ifs) / tanh(0.5 * β * ω)) * sin(ω * t) / π
        end
    end
    S, err1 = quadgk(integrand_re, b.lb * icm2ifs, b.ub * icm2ifs; rtol=b.rtol, atol=b.atol)

    # Imaginary part 
    integrand_im(ω) = - ω * b.specdens(ω; scale=icm2ifs) * cos(ω * t) / π
    A, err2 = quadgk(integrand_im, b.lb * icm2ifs, b.ub * icm2ifs; rtol=b.rtol, atol=b.atol)
    return S + 1.0im * A, max(err1, err2)
end

function (b::BosonicBCF_dt)(t::Real; ret_err::Bool=false)
    val, err = _compute_bcf_dt(b, t)
    ret_err ? (val, err) : val
end

# Fermionic BCF
abstract type FermionicBathCorrelationFunction <: BathCorrelationFunction end

#
struct FermionicBCF_Minus <: FermionicBathCorrelationFunction 
    specdens::Function
    Temp::Real
    ub::Real
end

# 
struct FermionicBCF_Plus <: FermionicBathCorrelationFunction 
    specdens::Function
    Temp::Real
    ub::Real
end

# Fermionic BCF Minus
function (b::FermionicBCF_Minus)(t::Real) :: ComplexF64
    # Real part of the BCF
    μ = b.ChemPot
    integrand_re = (ε::Float64) -> begin
        if b.Temp == 0.0
            b.specdens(ω; scale=icm2ifs) * cos(ε * t) / 2π
        else
            β = ħ * 1e15 / (kb * b.Temp)
            (b.specdens(ω; scale=icm2ifs) / tanh(0.5 * β * ω)) * cos(ω * t) / π
        end
    end
    Cr, err1 = quadgk(ω -> integrand_re(ω), b.lb * icm2ifs, b.ub * icm2ifs; atol=1e-12)

    # Imaginary part of the BCF
    integrand_im(ω) = -b.specdens(ω; scale=icm2ifs) * sin(ω * t) / π
    Ci, err2 = quadgk(ω -> integrand_im(ω), b.lb * icm2ifs, b.ub * icm2ifs; atol=1e-12)

    return Cr + 1.0im * Ci, max(err1, err2)
end

# Fermionic BCF
function (b::FermionicBCF_Plus)(t::Real) :: ComplexF64
    # Real part of the BCF
    integrand_re = (ω::Float64) -> begin
        if b.Temp == 0.0
            b.specdens(ω; scale=icm2ifs) * cos(ω * t) / π
        else
            β = ħ * 1e15 / (kb * b.Temp)
            (b.specdens(ω; scale=icm2ifs) / tanh(0.5 * β * ω)) * cos(ω * t) / π
        end
    end
    Cr, err1 = quadgk(ω -> integrand_re(ω), b.lb * icm2ifs, b.ub * icm2ifs; atol=1e-12)

    # Imaginary part of the BCF
    integrand_im(ω) = -b.specdens(ω; scale=icm2ifs) * sin(ω * t) / π
    Ci, err2 = quadgk(ω -> integrand_im(ω), b.lb * icm2ifs, b.ub * icm2ifs; atol=1e-12)

    return Cr + 1.0im * Ci, max(err1, err2)
end