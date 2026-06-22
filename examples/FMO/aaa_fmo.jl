include("../plot.jl")
using LinearAlgebra
using RationalFunctionApproximation
using Serialization
using DelimitedFiles

data = readdlm("fmo_smoothed.txt")
col1 = data[:, 1]
col2 = data[:, 2] 
ω = Float64.(col1) 
J = Float64.(col2)
J[J .< 0] .= 0.0 

println("AAA started")
r = aaa(ω, J; tol=1e-9, max_degree=1500, lookahead=400)
println("degree: ", length(r.nodes))
err = norm(r.(ω) - J)
println("error: ", err)

open("r_fmo.bin", "w") do io
    serialize(io, r)
end

pols, res, rinf = AAAPartialFraction(r)
Jrat = RationalSD(pols, res, rinf)
err = norm(Jrat.(ω) - J)
println("error (RationalSD): ", err)

open("jrat_fmo.bin", "w") do io
    serialize(io, Jrat)
end

Jrat = RationalSD(pols, res)
err = norm(Jrat.(ω) - J)
println("error (RationalSD wo rinf): ", err)

