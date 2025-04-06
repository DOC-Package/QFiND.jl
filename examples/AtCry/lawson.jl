using LinearAlgebra
using RationalFunctionApproximation
using Serialization
using DelimitedFiles
using QFiND

data = readdlm("atcry_smoothed.txt")
col1 = data[:, 1]
col2 = data[:, 2] 
ω = Vector(Float64.(col1)) 
J = Vector(Float64.(col2))
J[J .< 0.0] .= 0.0

r = open("r_atcry.bin", "r") do io
    deserialize(io)
end

err = norm(r.(ω) - J)
println("error: ", err)
#r = mylawson(ω, J, r, 10)

r = lawson(ω, J, r, 5)
err = norm(r.(ω) - J)
println("error: ", err)

open("r_atcry_lawson.bin", "w") do io
    serialize(io, r)
end
