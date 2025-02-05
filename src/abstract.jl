abstract type Approximation <:Function end

# Discretization methods
abstract type Discretization <: Approximation end
struct DiscrID <: Discretization end
struct DiscrBSDO <: Discretization end

# Decomposition methods
abstract type Decomposition <: Approximation end
struct DecompESPRIT <: Decomposition end
struct DecompMP <: Decomposition end




abstract type DiscreteDataSet end
struct DiscreteDataSetID <: DiscreteDataSet
    time :: Vector{Float64}
    freq :: Vector{Float64}
    qnsd :: Vector{Float64}
    bcf :: Vector{ComplexF64}
end

function DiscreteData(
    method::DiscrID, 
    sbeta::Function, 
    bcf::Function, 
    Ω_min::Real, 
    Ω_max::Real, 
    T_c::Real; 
    n_freq::Int=1000, 
    n_time::Int=200) :: DiscreteDataSetID

    # Generate the frequencies
    ω = collect(range(Ω_min, Ω_max, length=n_freq))
    ω = ω .* icm2ifs
    # Compute the quantum noise spectral density
    S = sbeta.(ω; scale=icm2ifs)
    S = max.(S, 0.0)
    # Generate the time grid
    t = collect(range(0, T_c, length=n_time))
    # Compute the bath correlation function
    bcf = bcf.(t)

    return DiscreteDataSetID(t, ω, S, bcf)
end
