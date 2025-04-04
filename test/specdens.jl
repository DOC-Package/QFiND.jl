using Test
using QFiND

@testset "specdens.jl" begin
    
    # Exp
    s = 1.0
    alpha = 1.0
    gamc = 50.0
    # BO
    Ω = 1400.0
    Γ = 200.0
    λ = 600.0

    sdens = TannorMeyerSD(Ω, Γ, λ)
    rene = reorganization_energy(sdens)
    @test rene ≈ λ 
    println("TM reorganization energy: ", rene)
    sdens = BrownianSD(Ω, Γ, λ)
    rene = reorganization_energy(sdens)
    @test rene ≈ λ 
    println("BO reorganization energy: ", rene)
end
