"""
function esprit_decomp(bcf::Function, N_t::Integer, tc::Real, eps::Real) 
    ef = esprit(bcf, 0.0, tc, N_t, eps)
    return (expon = ef.expon / icm2ifs, coeff = ef.coeff / icm2ifs^2.0)
end 

function DrudeHighT(specdens::DrudeSD, Temp::Float64; scale::Float64=1.0)
    β = ħ * 1e15 / (kb * Temp)
    γ = specdens.γ .* scale
    λ = specdens.λ .* scale
    c = 2.0 * λ / β - 1.0im * λ * γ
    return (expon = γ, coeff = c)
end

function BrownianHighT(specdens::BrownianSD, Temp::Float64; scale::Float64=1.0)
    β = ħ * 1e15 / (kb * Temp)
    Ω = specdens.Ω .* scale
    Γ = specdens.Γ .* scale
    λ = specdens.λ .* scale
    a = zeros(ComplexF64, 2 * length(Ω))
    c = zeros(ComplexF64, 2 * length(Ω))
    for i in eachindex(Ω)
        omega1 = sqrt(Ω[i]^2 - Γ[i]^2 / 4.0)
        gp = Γ[i] / 2.0 - im * omega1 
        gm = Γ[i] / 2.0 + im * omega1
        a[2*i-1] = gp
        a[2*i] = gm
        c[2*i-1] = - im * λ[i] * Ω[i]^2 / (β * omega1 * gp) - λ[i] * Ω[i]^2 / (2.0 * omega1)
        c[2*i] = im * λ[i] * Ω[i]^2 / (β * omega1 * gm) + λ[i] * Ω[i]^2 / (2.0 * omega1)
    end
    return (expon = a, coeff = c)
end
"""
