using QFiND
using Test

@testset "QFiND.jl" begin
    include("specdens.jl")
    include("bo.jl")
    include("bsdo.jl")
    include("id.jl")
end
