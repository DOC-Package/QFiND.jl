using QFiND
using Test

@testset "QFiND.jl" begin
    include("bo.jl")
    include("bsdo.jl")
    include("id.jl")
end
