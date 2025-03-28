abstract type BathCorrelationFunction <: Function end

# 
abstract type BosonicBathCorrelationFunction <: BathCorrelationFunction end

# Bosonic BCF with real and imaginary parts
struct BosonicBCF <: BosonicBathCorrelationFunction 
    specdens::SpectralDensity
    Temp::Real
    uplim::Real
end
BosonicBCF(specdens::SpectralDensity, Temp::Real) = BosonicBCF(specdens, Temp, Inf)

# Bosonic BCF with real part only
struct BosonicBCF_Real <: BosonicBathCorrelationFunction 
    specdens::SpectralDensity
    Temp::Real
    uplim::Real
end
BosonicBCF_Real(specdens::SpectralDensity, Temp::Real) = BosonicBCF_Real(specdens, Temp, Inf)

# Bosonic BCF with imaginary part only
struct BosonicBCF_Imag <: BosonicBathCorrelationFunction 
    specdens::SpectralDensity
    Temp::Real
    uplim::Real
end
BosonicBCF_Imag(specdens::SpectralDensity, Temp::Real) = BosonicBCF_Imag(specdens, Temp, Inf)

struct BosonicBCF_dt <: BosonicBathCorrelationFunction 
    specdens::SpectralDensity
    Temp::Real
    uplim::Real
end
BosonicBCF_dt(specdens::SpectralDensity, Temp::Real) = BosonicBCF_dt(specdens, Temp, Inf)


# Bosonic BCF
function (b::BosonicBCF)(t::Real) :: ComplexF64
    # Real part of the BCF
    integrand_re = (ω::Float64) -> begin
        if b.Temp == 0.0
            b.specdens(ω; scale=icm2ifs) * cos(ω * t) / π
        else
            β = ħ * 1e15 / (kb * b.Temp)
            (b.specdens(ω; scale=icm2ifs) / tanh(0.5 * β * ω)) * cos(ω * t) / π
        end
    end
    S, err1 = quadgk(ω -> integrand_re(ω), 0.0, b.uplim * icm2ifs; atol=1e-12)

    # Imaginary part of the BCF
    integrand_im(ω) = -b.specdens(ω; scale=icm2ifs) * sin(ω * t) / π
    A, err2 = quadgk(ω -> integrand_im(ω), 0.0, b.uplim * icm2ifs; atol=1e-12)

    return S + 1.0im * A
end

# Real part of a BCF for a time t.
function (b::BosonicBCF_Real)(t::Real) :: Float64
    integrand_re = (ω::Float64) -> begin
        if b.Temp == 0.0
            b.specdens(ω; scale=icm2ifs) * cos(ω * t) / π
        else
            β = ħ * 1e15 / (kb * b.Temp)
            (b.specdens(ω; scale=icm2ifs) / tanh(0.5 * β * ω)) * cos(ω * t) / π
        end
    end
    S, err = quadgk(ω -> integrand_re(ω), 0.0, b.uplim * icm2ifs; atol=1e-12)
    return S
end

# Imaginary part of a BCF for a time t.
function (b::BosonicBCF_Imag)(t::Real) :: Float64
    integrand_im(ω) = -b.specdens(ω; scale=icm2ifs) * sin(ω * t) / π
    A, err = quadgk(ω -> integrand_im(ω), 0.0, v.uplim * icm2ifs; atol=1e-12)
    return A
end

function reorganization_energy(sd::SpectralDensity; uplim::Real=Inf) :: Float64
    integrand(ω) = sd(ω) / ω / π
    E_r, err = quadgk(ω -> integrand(ω), 0.0, uplim; atol=1e-12)
    return E_r
end

# Bosonic BCF
function (b::BosonicBCF_dt)(t::Real) :: ComplexF64
    # Real part of the BCF
    integrand_re = (ω::Float64) -> begin
        if b.Temp == 0.0
            - ω * b.specdens(ω; scale=icm2ifs) * sin(ω * t) / π
        else
            β = ħ * 1e15 / (kb * b.Temp)
            - ω * (b.specdens(ω; scale=icm2ifs) / tanh(0.5 * β * ω)) * sin(ω * t) / π
        end
    end
    S, err1 = quadgk(ω -> integrand_re(ω), 0.0, b.uplim * icm2ifs; atol=1e-12)

    # Imaginary part of the BCF
    integrand_im(ω) = - ω * b.specdens(ω; scale=icm2ifs) * cos(ω * t) / π
    A, err2 = quadgk(ω -> integrand_im(ω), 0.0, b.uplim * icm2ifs; atol=1e-12)

    return S + 1.0im * A
end

# Fermionic BCF
abstract type FermionicBathCorrelationFunction <: BathCorrelationFunction end

#
struct FermionicBCF_Minus <: FermionicBathCorrelationFunction 
    specdens::Function
    Temp::Real
    uplim::Real
end

# 
struct FermionicBCF_Plus <: FermionicBathCorrelationFunction 
    specdens::Function
    Temp::Real
    uplim::Real
end

# Fermionic BCF Minus
function (b::FermionicBCF_Minus)(t::Real) :: ComplexF64
    # Real part of the BCF
    μ = b.ChemPot
    integrand_re = (ε::Float64) -> begin
        if b.Temp == 0.0
            b.specdens(ω; scale=icm2ifs) * () * cos(ε * t) / 2π
        else
            β = ħ * 1e15 / (kb * b.Temp)
            (b.specdens(ω; scale=icm2ifs) / tanh(0.5 * β * ω)) * cos(ω * t) / π
        end
    end
    Cr, err1 = quadgk(ω -> integrand_re(ω), 0.0, b.uplim * icm2ifs; atol=1e-12)

    # Imaginary part of the BCF
    integrand_im(ω) = -b.specdens(ω; scale=icm2ifs) * sin(ω * t) / π
    Ci, err2 = quadgk(ω -> integrand_im(ω), 0.0, b.uplim * icm2ifs; atol=1e-12)

    return Cr + 1.0im * Ci
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
    Cr, err1 = quadgk(ω -> integrand_re(ω), 0.0, b.uplim * icm2ifs; atol=1e-12)

    # Imaginary part of the BCF
    integrand_im(ω) = -b.specdens(ω; scale=icm2ifs) * sin(ω * t) / π
    Ci, err2 = quadgk(ω -> integrand_im(ω), 0.0, b.uplim * icm2ifs; atol=1e-12)

    return Cr + 1.0im * Ci
end