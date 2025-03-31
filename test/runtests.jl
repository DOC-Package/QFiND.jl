using QFiND
using Test

@testset "QFiND.jl" begin
    include("bo.jl")
    include("bsdo.jl")
    include("svd.jl")
    include("id.jl")
    include("idc.jl")
end
