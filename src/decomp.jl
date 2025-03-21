function esprit_decomp(bcf::Function, N_t::Integer, tc::Real, eps::Real) 
    ef = esprit(bcf, 0.0, tc, N_t, eps)
    return (expon = ef.expon / icm2ifs, coeff = ef.coeff / icm2ifs^2.0)
end 

function Drude_HighT(specdens::DrudeSD, Temp::Float64; scale::Float64=1.0)
    β = ħ * 1e15 / (kb * Temp)
    γ = specdens.γ .* scale
    λ = specdens.λ .* scale
    a = zeros(ComplexF64, 2 * length(γ))
    c = zeros(ComplexF64, 2 * length(γ))
    
end

function Brownian_HighT(specdens::BrownianSD; scale::Float64=1.0)
    Ω = specdens.Ω .* scale
    Γ = specdens.Γ .* scale
    λ = specdens.λ .* scale
    a = zeros(ComplexF64, 2 * length(Ω))
    c = zeros(ComplexF64, 2 * length(Ω))
    for i in eachindex(Ω)
        a[2*i-1] = Ω[i] + im * Γ[i]
        a[2*i] = Ω[i] - im * Γ[i]
        c[2*i-1] = λ[i] * Γ[i]^2 / (im * Γ[i])
        c[2*i] = -λ[i] * Γ[i]^2 / (im * Γ[i])
    end
    return (expon = a, coeff = c)
end

function TannorMeyer_HighT(specdens::TannorMeyerSD; scale::Float64=1.0)
    Ω = specdens.Ω .* scale
    Γ = specdens.Γ .* scale
    λ = specdens.λ .* scale
    a = zeros(ComplexF64, 2 * length(Ω))
    c = zeros(ComplexF64, 2 * length(Ω))
    for i in eachindex(Ω)
        a[2*i-1] = Ω[i] + im * (Γ[i] / 2)
        a[2*i] = Ω[i] - im * (Γ[i] / 2)
        c[2*i-1] = λ[i] * (Γ[i] / 2)^2 / (im * Γ[i])
        c[2*i] = -λ[i] * (Γ[i] / 2)^2 / (im * Γ[i])
    end
    return (expon = a, coeff = c)
end
