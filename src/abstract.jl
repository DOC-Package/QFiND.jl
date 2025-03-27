abstract type Approximation <:Function end

# Discretization methods
abstract type Discretization <: Approximation end
struct DiscrID <: Discretization end
struct DiscrBSDO <: Discretization end

abstract type Decomposition <: Approximation end
struct DecompID <: Decomposition end
struct DecompSVD <: Decomposition end
struct DecompESPRIT <: Decomposition end

# 
abstract type InitialDataSet end
struct InitialDataSetBSDO <: InitialDataSet
    freq :: Vector{Float64}
    qnsd :: Function
end
struct InitialDataSetID <: InitialDataSet
    freq :: Vector{Float64}
    time :: Vector{Float64}
    qnsd :: Vector{Float64}
    bcf :: Vector{ComplexF64}
end
struct InitialDataSetSVD <: InitialDataSet
    freq :: Vector{Float64}
    time :: Vector{Float64}
    weights :: Vector{Float64}
    qnsd :: Vector{Float64}
    bcf :: Vector{ComplexF64}
    dC :: Vector{ComplexF64}
end

function InitialData(
    method::DiscrBSDO, 
    sbeta::Function, 
    Ω_min::Real, 
    Ω_max::Real; 
    n_lanczos::Int=1000 
    ) :: InitialDataSetBSDO

    ω = collect(range(Ω_min, Ω_max, length=n_lanczos))

    return InitialDataSetBSDO(ω, sbeta)
end

function InitialData(
    method::DiscrID, 
    sbeta::Function, 
    bcf::Function, 
    Ω_min::Real, 
    Ω_max::Real, 
    T_c::Real; 
    n_freq::Int=1000, 
    n_time::Int=500) :: InitialDataSetID

    ω = collect(range(Ω_min, Ω_max, length=n_freq))
    ω = ω .* icm2ifs
    S = sbeta.(ω; scale=icm2ifs)
    S = max.(S, 0.0)
    t = collect(range(0, T_c, length=n_time))
    bcf = bcf.(t)
    return InitialDataSetID(ω, t, S, bcf)
end

function InitialData(
    method::DecompID, 
    sbeta::Function, 
    bcf::Function, 
    Ω_min::Real, 
    Ω_max::Real, 
    T_c::Real; 
    n_freq::Int=1000, 
    n_time::Int=500) :: InitialDataSetID

    ω = collect(range(Ω_min, Ω_max, length=n_freq))
    ω = ω .* icm2ifs
    S = sbeta.(ω; scale=icm2ifs)
    S = max.(S, 0.0)
    t = collect(range(0, T_c, length=n_time))
    bcf = bcf.(t)
    return InitialDataSetID(ω, t, S, bcf)
end

function InitialData(
    method::DecompSVD, 
    sbeta::Function, 
    bcf::Function, 
    dC::Function,
    Ω_min::Real, 
    Ω_max::Real, 
    T_c::Real; 
    n_freq::Int=1000, 
    n_time::Int=500) :: InitialDataSetSVD

    ω = collect(range(Ω_min, Ω_max, length=n_freq))
    ω = ω .* icm2ifs
    w = ones(length(ω))
    S = sbeta.(ω; scale=icm2ifs)
    S = max.(S, 0.0)
    t = collect(range(0, T_c, length=n_time))
    bcf = bcf.(t)
    dC = dC.(t)
    return InitialDataSetSVD(ω, t, w, S, bcf, dC)
end

