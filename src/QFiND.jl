module QFiND

include("bsdo.jl")
include("constants.jl")
include("corrfunc.jl")
include("eval.jl")
include("id.jl")
include("id_sub.jl")
include("lanczos.jl")
include("specdens.jl")

using LinearAlgebra
using LowRankApprox
using NonNegLeastSquares
using QuadGK 
using .Constants

export bsdo_discr, make_bcf, calc_error, id_discr, make_sdens, make_sbeta, icm2ifs

end
