include("smooth.jl")

data = readdlm("fmo.txt")
col1 = data[:, 1]
col2 = data[:, 2]
ω = Float64.(col1)
J = Float64.(col2)
S = smoothing(ω, J)

