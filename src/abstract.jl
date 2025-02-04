abstract type DiscretizationMethod <: Function end
struct ID <: DiscretizationMethod end
struct BSDO <: DiscretizationMethod end

abstract type DiscreteData end
struct DiscreteDataID <: DiscreteData
    freq :: Vector{Float64}
    coef :: Vector{ComplexF64}
end

struct DiscreteDataBSDO <: DiscreteData
    freq :: Vector{Float64}
    qnsd :: Vector{ComplexF64}
end

function prepare_data(BSDO::BSDO, sbeta::Function, Omega_min::Real, Omega_max::Real, M_sp::Int; nlanczos::Int=1000) :: DiscreteDataBSDO
    # Generate the frequencies
    ω = range(Omega_min, Omega_max, length=nlanczos) |> collect

    # Compute the quantum noise spectral density
    S = sbeta.(ω)
    S = max.(S, 0.0)
end
