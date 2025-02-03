export make_bcf

function make_bcf(specdens::Function, Temp::Real, Ω_c::Real)
    return t -> begin
        return bcf_re(specdens, Temp, Ω_c, t) + 1.0im * bcf_im(specdens, Ω_c, t)
    end
end

# Real part of a BCF for a time t.
function bcf_re(specdens::Function, Temp::Real, Ω_c::Real, t::Real)
    if Temp == 0.0
        integrand_0K(w) = specdens(w) * cos(w*t) / π
        (res, err) = quadgk(w -> integrand_0K(w), 0.0, Ω_c * icm2ifs; atol=1e-12)
    else
        beta = hbar * 1e15 / (kb * Temp)
        integrand(w) = (specdens(w) / tanh(0.5*beta*w)) * cos(w*t) / π
        (res, err) = quadgk(w -> integrand(w), 0.0, Ω_c * icm2ifs; atol=1e-12)
    end
    return res
end

# Imaginary part of a BCF for a time t.
function bcf_im(specdens::Function, Ω_c::Real, t::Real)
    integrand(w) = - specdens(w) * sin(w*t) / π
    (res, err) = quadgk(w -> integrand(w), 0.0, Ω_c * icm2ifs; atol=1e-12)
    return res
end
