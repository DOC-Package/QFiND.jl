include("../plot.jl")
using LinearAlgebra
using QFiND
using RationalFunctionApproximation
using Serialization
using DelimitedFiles

data = readdlm("atcry_smoothed.txt")
col1 = data[:, 1]
col2 = data[:, 2] 
ω = Float64.(col1) 
J = Float64.(col2)
J[J .< 0.0] .= 0.0
#idx = findall(ω .< 1.0)
#J[idx] .= 0.0

println("AAA started")
r = aaa(ω, J; tol=1e-12, max_degree=1500, lookahead=200)
println("degree: ", length(r.nodes))
err = norm(r.(ω) - J)
println("error: ", err)

sdens = AAAfittedSD(r)
plot_qnsd(sdens, 0.0, 4000, 4000, "atcry_aaa.png")

open("r_atcry.bin", "w") do io
    serialize(io, r)
end

