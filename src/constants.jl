module Constants

export icm2ifs, ħ, kb

const icm2ifs = 299792458.0 * 1e2 * 2.0 * π * 1e-15   # c * 100 * 2π * 1e-15
const ħ = (6.62607015e-34) / (2.0 * π)          # Dirac's constant
const kb      = 1.380649e-23                          # Boltzmann's constant

end 
