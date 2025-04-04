include("../plot.jl")
using RationalFunctionApproximation
using Serialization
using DelimitedFiles

data = readdlm("fmo_smoothed.txt")
col1 = data[:, 1]
col2 = data[:, 2] 
ω = Float64.(col1) 
J = Float64.(col2)
J[J .< 0.0] .= 0.0

r = open("r_fmo.bin", "r") do io
    deserialize(io)
end

plot_qnsd(ω, r.(ω), J, "aaa_fmo.png")