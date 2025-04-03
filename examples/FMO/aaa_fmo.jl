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

"AAA started"
r = aaa(ω, J; tol=1e-14, max_degree=1200, lookahead=200)
println("degree: ", length(r.nodes))
err = norm(r.(ω) - J)
println("error: ", err)

open("r_fmo.bin", "w") do io
    serialize(io, r)
end

