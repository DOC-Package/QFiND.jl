using Test
using QFiND

@testset "specdens.jl" begin
    
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
