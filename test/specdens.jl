using Test
using QFiND

@testset "specdens.jl" begin

    @testset "Ohmic exponential cutoff" begin
        γ = 4.0
        α = 2.5
        ohmic = OhmicExpSD(γ; alpha=α)
        powerlaw = PowerLawExpSD(1.0, γ; alpha=α)

        @test ohmic(1.5) ≈ π * α / γ * 1.5 * exp(-1.5 / γ)
        @test ohmic(-1.5) ≈ -ohmic(1.5)
        @test ohmic(1.5) ≈ powerlaw(1.5)
        @test reorganization_energy(ohmic) ≈ α
    end

    @testset "Ohmic Gaussian cutoff" begin
        γ = 4.0
        α = 2.5
        ohmic = OhmicGaussSD(γ; alpha=α)

        @test ohmic(1.5) ≈ π * α / γ * 1.5 * exp(-(1.5 / γ)^2)
        @test ohmic(-1.5) ≈ -ohmic(1.5)
        @test reorganization_energy(ohmic) ≈ α * sqrt(π) / 2

        normalized = OhmicGaussSD(γ; reorgene=α)
        @test reorganization_energy(normalized) ≈ α
    end
    
    # Exp
    s = 0.5
    lam = 50.0
    gamc = 50.0
    # BO
    Ω = 1400.0
    Γ = 200.0
    λ = 600.0

    sdens = PowerLawExpSD(s, gamc; reorgene=lam)
    rene = reorganization_energy(sdens)
    @test rene ≈ lam
    println("PL reorganization energy: ", rene)
    sdens = TannorMeyerSD(Ω, Γ, λ)
    rene = reorganization_energy(sdens)
    @test rene ≈ λ 
    println("TM reorganization energy: ", rene)
    sdens = BrownianSD(Ω, Γ, λ)
    rene = reorganization_energy(sdens)
    @test rene ≈ λ 
    println("BO reorganization energy: ", rene)
end
