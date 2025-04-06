include("../plot.jl")
using QFiND
using LinearAlgebra
using RationalFunctionApproximation
  
# spectral density
s = 0.25
alpha = 50.0
gamc = 50.0
Temp = 0.0
Ω_min = 0.0
Ω_max = 600.0
N_ω = 2000

sdens = PowerLawExpSD(s, gamc; alpha=alpha)

ω = collect(range(Ω_min, Ω_max, length=N_ω)) 
J = sdens.(ω)
r = aaa(ω, J; tol=1e-5, max_degree=100, lookahead=50)

println("degree: ", length(r.nodes))
err = norm(r.(ω) - J)
println("error: ", err)

r = mylawson(ω, J, r, 10)
err = norm(r.(ω) - J)
println("error: ", err)
#