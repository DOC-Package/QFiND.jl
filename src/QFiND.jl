module QFiND

export icm2ifs
export sumexp, equispaced_grid, calc_error
export DiscretizationMethod, DiscrID, DiscrBSDO
export DiscreteDataSet, DiscreteDataSetID, DiscreteData
export SpectralDensity, PowerLawExpSD, TannorMeyerSD, BrownianSD, AAAfittedSD
export QuantumNoiseSpectralDensity, BosonicQNSD, FermionicQNSD
export FermionicQNSD_Plus, FermionicQNSD_Minus
export BathCorrelationFunction, BosonicBathCorrelationFunction, FermionicBathCorrelationFunction
export BosonicBCF, BosonicBCF_Real, BosonicBCF_Imag
export FermionicBCF_Plus, FermionicBCF_Minus
export bsdo_discr, id_discr

import RationalFunctionApproximation: Barycentric

include("constants.jl")
include("abstract.jl")
include("utils.jl")
include("specdens.jl")
include("corrfunc.jl")
include("bsdo.jl")
include("id_sub.jl")
include("id.jl")
include("eval.jl")

using LinearAlgebra
using LowRankApprox
using NonNegLeastSquares
using QuadGK 
using RationalFunctionApproximation

end
