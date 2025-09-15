# Define Abstract Type for Spectral Density Models
abstract type SpectralDensity <: Function end

# SumSD type
struct SumSD <: SpectralDensity
    sd1 :: SpectralDensity
    sd2 :: SpectralDensity
end

function (s::SumSD)(ω::Float64; scale::Float64=1.0) :: Float64
    return s.sd1(ω; scale=scale) + s.sd2(ω; scale=scale)
end

function (s::SumSD)(ω::ComplexF64; scale::Float64=1.0) :: ComplexF64
    return s.sd1(ω; scale=scale) + s.sd2(ω; scale=scale)
end
# Overloading the + operator for SpectralDensity types
Base.:+(sd1::SpectralDensity, sd2::SpectralDensity) = SumSD(sd1, sd2)

function sd_poles(sd::SumSD; scale::Float64=1.0)
    poles1 = sd_poles(sd.sd1; scale=scale)
    poles2 = sd_poles(sd.sd2; scale=scale)
    return vcat(poles1, poles2)
end

function sd_residues(sd::SumSD; scale::Float64=1.0)
    residues1 = sd_residues(sd.sd1; scale=scale)
    residues2 = sd_residues(sd.sd2; scale=scale)
    return vcat(residues1, residues2)
end

# Power-law with exponential cutoff
struct PowerLawExpSD <: SpectralDensity
    s :: Float64   # exponent
    γ :: Float64   # cutoff frequency
    α :: Float64   # coefficient 
    reorgene :: Float64   # reorganization energy
end

function PowerLawExpSD(s::Float64, γ::Float64; alpha::Union{Float64,Nothing}=nothing, reorgene::Union{Float64,Nothing}=nothing)
    if (alpha === nothing) == (reorgene === nothing)
        throw(ArgumentError("Please provide either α or reorgene, but not both."))
    end

    if alpha !== nothing
        reorgene = gamma(s)
        return PowerLawExpSD(s, γ, alpha, reorgene)
    else
        alpha = reorgene / gamma(s)
        return PowerLawExpSD(s, γ, alpha, reorgene)
    end
end

# Tannor-Meyer
struct TannorMeyerSD <: SpectralDensity
    Ω :: Vector{Float64}  # vector of central frequencies
    Γ :: Vector{Float64}  # vector of widths
    λ :: Vector{Float64}  # vector of intensities
end
TannorMeyerSD(Ω::Real, Γ::Real, λ::Real) = TannorMeyerSD([Float64(Ω)], [Float64(Γ)], [Float64(λ)])

# Brownian oscillator
struct BrownianSD <: SpectralDensity
    Ω :: Vector{Float64}  # vector of local frequencies
    Γ :: Vector{Float64}  # vector of damping coefficients
    λ :: Vector{Float64}  # vector of coupling strengths
end
BrownianSD(Ω::Real, Γ::Real, λ::Real) = BrownianSD([Float64(Ω)], [Float64(Γ)], [Float64(λ)])

struct DrudeSD <: SpectralDensity
    γ :: Float64
    λ :: Float64
end

struct DiscreteSD <: SpectralDensity
    freq :: Vector{Float64}
    weight :: Vector{Float64}
end
DiscreteSD(ω::Real, g::Real) = DiscreteSD([Float64(ω)], [Float64(g)])

struct AAAfittedSD <: SpectralDensity
    bary::Barycentric
    reorgene::Float64
    scale_reorgene::Float64
end
AAAfittedSD(bary::Barycentric; lb::Real=0.0, ub::Real=Inf) = AAAfittedSD(bary, reorganization_energy(bary; lb=lb,ub=ub), 1.0)
function AAAfittedSD(bary::Barycentric, reorgene::Float64; lb::Real=0.0, ub::Real=Inf)
    rawene = reorganization_energy(bary; lb=lb,ub=ub)
    scale_reorgene = reorgene / rawene
    return AAAfittedSD(bary, reorgene, scale_reorgene)
end

struct RationalSD <: SpectralDensity
    poles::Vector{ComplexF64}
    residues::Vector{ComplexF64}
end

struct WideBandSD <: SpectralDensity
    Γ :: Float64
end

function reorganization_energy(sd::SpectralDensity; lb::Real=0.0, ub::Real=Inf) :: Float64
    integrand(ω) = sd(ω) / ω / π
    E_r, err = quadgk(ω -> integrand(ω), lb, ub; atol=1e-12)
    return E_r
end

function reorganization_energy(bary::Barycentric; lb::Real=0.0, ub::Real=Inf)
    integrand(ω) = (evaluate(bary, ω)) / ω / π
    E_r, err = quadgk(ω -> integrand(ω), lb, ub; atol=1e-12)
    return E_r
end

function reorganization_energy(ω::Vector{Float64}, g::Vector{Float64})
    return sum(g.^2.0 ./ abs.(ω)) / 2.0
end
reorganization_energy(ω::Real, g::Real) = reorganization_energy([Float64(ω)], [Float64(g)])

"""
    (specdens::PowerLawExpSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Power-law with exponential cutoff model.
"""
function (specdens::PowerLawExpSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    s = specdens.s
    γ = specdens.γ * scale
    α = specdens.α * scale
    return sgn * π * α * γ^(-s) * ω^s * exp(-ω / γ)
end


"""
    (specdens::TannorMeyerSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Tannor-Meyer model.
"""
function (specdens::TannorMeyerSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ωa = abs(ω)
    Ωt = specdens.Ω 
    Γt = specdens.Γ 
    λt = specdens.λ 
    res = 0.0
    @inbounds for i in eachindex(Ωt)
        Ωs = Ωt[i] * scale
        Γs = Γt[i] * scale
        λs = λt[i] * scale
        p    = 4.0 * Γs * λs * (Ωs^2 + Γs^2)
        deno = ((ωa + Ωs)^2 + Γs^2) * ((ωa - Ωs)^2 + Γs^2)
        res += p * ωa / deno
    end
    return sgn * res
end

function sd_poles(sd::TannorMeyerSD)
    poles = ComplexF64[]
    for (Ω, Γ, _) in zip(sd.Ω, sd.Γ, sd.λ)
        push!(poles, Ω + im*(Γ/2))
        push!(poles, Ω - im*(Γ/2))
    end
    return poles
end

function sd_residues(sd::TannorMeyerSD)
    res = ComplexF64[]
    for (Ω, Γ, λ) in zip(sd.Ω, sd.Γ, sd.λ)
        res1 = λ * (Γ/2)^2 / (im*Γ)
        res2 = -λ * (Γ/2)^2 / (im*Γ)
        push!(res, res1)
        push!(res, res2)
    end
    return res
end

"""
    (specdens::BrownianSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Brownian oscillator model.
"""
function (specdens::BrownianSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ωa = abs(ω)
    Ωt = specdens.Ω 
    Γt = specdens.Γ 
    λt = specdens.λ 
    res = 0.0
    @inbounds for i in eachindex(Ωt)
        Ωs = Ωt[i] * scale
        Γs = Γt[i] * scale
        λs = λt[i] * scale
        p    = 2.0 * Γs * λs * (Ωs^2)
        deno = (ωa^2 - Ωs^2)^2 + (Γs^2 * ωa^2)
        res += p * ωa / deno
    end
    return sgn * res
end

function (specdens::BrownianSD)(ω::ComplexF64; scale::Float64=1.0) :: ComplexF64
    Ωt = specdens.Ω 
    Γt = specdens.Γ 
    λt = specdens.λ 
    res = 0.0
    @inbounds for i in eachindex(Ωt)
        Ωs = Ωt[i] * scale
        Γs = Γt[i] * scale
        λs = λt[i] * scale
        p    = 2.0 * Γs * λs * (Ωs^2)
        deno = (ω^2 - Ωs^2)^2 + (Γs^2 * ω^2)
        res += p * ω / deno
    end
    return res
end

function sd_poles(sd::BrownianSD; scale::Float64=1.0)
    poles = ComplexF64[]
    Ωt = sd.Ω 
    Γt = sd.Γ 
    for (Ω, Γ) in zip(Ωt, Γt)
        Ω = Ω * scale
        Γ = Γ * scale
        delta = sqrt(Ω^2 - (Γ/2)^2)
        push!(poles, delta - im*(Γ/2))
        push!(poles, -delta - im*(Γ/2))
    end
    return poles
end

function sd_residues(sd::BrownianSD; scale::Float64=1.0)
    res = ComplexF64[]
    Ωt = sd.Ω
    Γt = sd.Γ
    λt = sd.λ
    for (Ω, Γ, λ) in zip(Ωt, Γt, λt)
        Ω = Ω * scale
        Γ = Γ * scale
        λ = λ * scale
        delta = sqrt(Ω^2 - (Γ/2)^2)
        res1 = - (λ * Ω^2) / (2im*delta)
        res2 =   (λ * Ω^2) / (2im*delta)
        push!(res, res1)
        push!(res, res2)
    end
    return res
end

"""
    (specdens::DrudeSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Drude-Lorentz model.
"""
function (specdens::DrudeSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ωa = abs(ω)
    γs = specdens.γ * scale 
    λs = specdens.λ * scale
    res = 2 * λs * γs * ωa / (ωa^2 + γs^2)
    return sgn * res
end

function (specdens::DrudeSD)(ω::ComplexF64; scale::Float64=1.0) :: ComplexF64
    γs = specdens.γ * scale 
    λs = specdens.λ * scale
    res = 2 * λs * γs * ω / (ω^2 + γs^2)
    return res
end

function sd_poles(sd::DrudeSD; scale::Float64=1.0)
    poles = ComplexF64[]
    γs = sd.γ * scale
    push!(poles, - 1im * γs)
    return poles
end

function sd_residues(sd::DrudeSD; scale::Float64=1.0)
    res = ComplexF64[]
    γs = sd.γ * scale
    λs = sd.λ * scale
    push!(res, λs * γs)
    return res
end

"""
    (specdens::AAAfittedSD)(ω::Float64; scale::Float64=1.0) -> Float64

Compute the spectral density for the Brownian oscillator model.
"""
function (specdens::AAAfittedSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω) / scale
    res = evaluate(specdens.bary, ω)
    return sgn * scale * res * specdens.scale_reorgene
end


function (specdens::RationalSD)(ω::Float64; scale::Float64=1.0) :: Float64
    sgn = sign(ω)
    ω = abs(ω)
    poles = specdens.poles .* scale
    residues = specdens.residues .* scale
    res = 0.0
    for i in eachindex(specdens.poles)
        res += residues[i] / (ω - poles[i])
    end
    res = real(res)
    if res < 0.0
        res = 0.0
    end
    return sgn * scale * res
end

# Define Abstract Type for the QNSD
abstract type QuantumNoiseSpectralDensity <: Function end

# Bosonic QNSD
struct BosonicQNSD <: QuantumNoiseSpectralDensity
    specdens :: SpectralDensity  
    Temp :: Float64          
end

struct BosonicQNSD_HighT <: QuantumNoiseSpectralDensity
    specdens :: SpectralDensity  
    Temp :: Float64          
end


function f_BE(ω::Float64, Temp::Float64)
    β = ħ * 1e15 / (kb * Temp)
    return 1.0 / (exp(β * icm2ifs * ω) - 1.0)
end


function (b::BosonicQNSD)(ω::Float64; scale::Float64=1.0) :: Float64
    if b.Temp == 0.0
        return (b.specdens)(ω; scale=scale)
    else
        β = ħ * 1e15 / (kb * b.Temp)
        factor = (1.0 / tanh(0.5 * β * ω * icm2ifs / scale) + 1.0) / 2.0
        return (b.specdens)(ω; scale=scale) * factor
    end
end

BosonicThermalBogoliubov(ω::Float64, g::Float64, Temp::Float64; scale::Float64=1.0) = BosonicThermalBogoliubov([Float64(ω)], [Float64(g)], Temp; scale=scale)
function BosonicThermalBogoliubov(ω::AbstractVector{Float64}, g::AbstractVector{Float64}, Temp::Float64; scale::Float64=1.0)
    β = ħ * 1e15 / (kb * Temp)
    g_p = g .* (1.0 ./ tanh.(0.5 * β .* ω .* icm2ifs / scale) .+ 1.0) ./ 2.0
    g_t = -g .* (1.0 ./ tanh.(-0.5 * β .* ω .* icm2ifs / scale) .+ 1.0) ./ 2.0
    freq = vcat(-ω, ω)
    coeff = vcat(sqrt.(g_t), sqrt.(g_p))
    perm = sortperm(freq)
    return freq[perm], coeff[perm]
end

BosonicQNSD_Discrete(ω::Float64, g::Float64, Temp::Float64; scale::Float64=1.0) = BosonicQNSD_Discrete([Float64(ω)], [Float64(g)], Temp; scale=scale)
function BosonicQNSD_Discrete(ω::AbstractVector{Float64}, g::AbstractVector{Float64}, Temp::Float64; scale::Float64=1.0) 
    β = ħ * 1e15 / (kb * Temp)
    g_p = g .* (1.0 ./ tanh.(0.5 * β .* ω .* icm2ifs / scale) .+ 1.0) ./ 2.0
    g_t = -g .* (1.0 ./ tanh.(-0.5 * β .* ω .* icm2ifs / scale) .+ 1.0) ./ 2.0
    expon = vcat(ω, -ω) 
    coeff = vcat(0.5 .* g_p, 0.5 .* g_t) 
    return expon, coeff
end

function (b::BosonicQNSD_HighT)(ω::Float64; scale::Float64=1.0) :: Float64
    if b.Temp == 0.0
        return (b.specdens)(ω; scale=scale) 
    else
        β = ħ * 1e15 / (kb * b.Temp)
        factor = 2.0 / (β * ω * icm2ifs / scale) / 2.0
        return (b.specdens)(ω; scale=scale) * factor
    end
end

# Fermionic QNSD
abstract type FermionicQNSD <: QuantumNoiseSpectralDensity end
struct FermionicQNSD_Plus <: FermionicQNSD 
    specdens :: SpectralDensity
    Temp :: Float64      
    ChemPot :: Float64
end
struct FermionicQNSD_Minus <: FermionicQNSD
    specdens :: SpectralDensity
    Temp :: Float64    
    ChemPot :: Float64  
end

# Discrete Gaussian Spectral Density
struct DiscreteGaussianSD <: SpectralDensity
    frequencies :: Vector{Float64}    
    amplitudes :: Vector{Float64}     
    width :: Float64                  
    reorgene :: Float64              # 再組織化エネルギー
end

"""
    DiscreteGaussianSD(frequencies, amplitudes, width)
    DiscreteGaussianSD(frequencies, amplitudes, width, reorgene)

離散データのスペクトル密度をガウシアンの和として表現します。

# Arguments
- `frequencies`: 周波数ベクトル
- `amplitudes`: 各周波数でのスペクトル密度値
- `width`: ガウシアン幅
- `reorgene`: 再組織化エネルギー（省略時は自動計算）

# Example
```julia
frequencies = [100.0, 200.0, 300.0]
amplitudes = [0.5, 0.8, 0.3]
width = 10.0
sd = DiscreteGaussianSD(frequencies, amplitudes, width)
# または再組織化エネルギーを指定
sd = DiscreteGaussianSD(frequencies, amplitudes, width, 0.1)
```
"""
# reorgeneを自動計算するコンストラクタ
function DiscreteGaussianSD(frequencies::Vector{Float64}, amplitudes::Vector{Float64}, width::Float64)
    # 離散版の再組織化エネルギー計算
    reorgene = sum(amplitudes ./ abs.(frequencies)) / π
    return DiscreteGaussianSD(frequencies, amplitudes, width, reorgene)
end

function (sd::DiscreteGaussianSD)(ω::Float64; scale::Float64=1.0) :: Float64
    ωs = ω * scale
    ws = sd.width * scale
    result = 0.0
    for i in eachindex(sd.frequencies)
        freq_scaled = sd.frequencies[i] * scale
        amp_scaled = sd.amplitudes[i] / scale^2  
        # Gaussian: A * ω * exp(-(ω-ω₀)²/(2σ²)) - 収束性向上のためωをかける
        result += amp_scaled * ωs * exp(-(ωs - freq_scaled)^2 / (2 * ws^2))
    end
    return result
end

# Discrete Lorentzian Spectral Density
struct DiscreteLorentzianSD <: SpectralDensity
    frequencies :: Vector{Float64}    
    amplitudes :: Vector{Float64}     
    width :: Float64                  
    reorgene :: Float64              # 再組織化エネルギー
end

"""
    DiscreteLorentzianSD(frequencies, amplitudes, width)
    DiscreteLorentzianSD(frequencies, amplitudes, width, reorgene)

離散データのスペクトル密度をLorentzianの和として表現します。

# Arguments
- `frequencies`: 周波数ベクトル
- `amplitudes`: 各周波数でのスペクトル密度値
- `width`: Lorentzian幅（半値全幅FWHM）
- `reorgene`: 再組織化エネルギー（省略時は自動計算）

# Example
```julia
frequencies = [100.0, 200.0, 300.0]
amplitudes = [0.5, 0.8, 0.3]
width = 10.0
sd = DiscreteLorentzianSD(frequencies, amplitudes, width)
# または再組織化エネルギーを指定
sd = DiscreteLorentzianSD(frequencies, amplitudes, width, 0.1)
```
"""
# reorgeneを自動計算するコンストラクタ
function DiscreteLorentzianSD(frequencies::Vector{Float64}, amplitudes::Vector{Float64}, width::Float64)
    reorgene = sum(amplitudes ./ abs.(frequencies)) / π
    return DiscreteLorentzianSD(frequencies, amplitudes, width, reorgene)
end

# 実数周波数での評価
function (sd::DiscreteLorentzianSD)(ω::Float64; scale::Float64=1.0) :: Float64
    ωs = ω * scale
    γs = sd.width * scale / 2.0  # FWHM を半値幅に変換
    result = 0.0
    for i in eachindex(sd.frequencies)
        freq_scaled = sd.frequencies[i] * scale
        amp_scaled = sd.amplitudes[i] / scale^2
        # Lorentzian: A * (γ/2) * ω / ((ω-ω₀)² + (γ/2)²)
        result += amp_scaled * γs * ωs / ((ωs - freq_scaled)^2 + γs^2)
    end
    return result
end

# 極とその留数を計算する関数
function sd_poles(sd::DiscreteLorentzianSD; scale::Float64=1.0)
    poles = ComplexF64[]
    γs = sd.width * scale / 2.0
    for freq in sd.frequencies
        freq_scaled = freq * scale
        # Lorentzianの極は ω₀ ± iγ/2
        push!(poles, freq_scaled + 1im * γs)
        push!(poles, freq_scaled - 1im * γs)
    end
    return poles
end

function sd_residues(sd::DiscreteLorentzianSD; scale::Float64=1.0)
    residues = ComplexF64[]
    γs = sd.width * scale / 2.0
    for amp in sd.amplitudes
        amp_scaled = amp / scale^2
        # Lorentzianの留数
        res = amp_scaled * γs / (2im * γs)  # = amp_scaled / (2im)
        push!(residues, res)
        push!(residues, -res)  # 共役極に対応
    end
    return residues
end
