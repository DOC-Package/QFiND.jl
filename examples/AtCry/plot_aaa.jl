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

r = open("r_atcry.bin", "r") do io
    deserialize(io)
end

plot_qnsd(ω, r.(ω), J, "./figure/aaa_atcry.png")