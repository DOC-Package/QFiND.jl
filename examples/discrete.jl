using QFiND
  
# spectral density
omega = 100.0
g = sqrt(300.0)
Temp = 100.0
freq, coeff = ThermalBogoliubov(omega, g, Temp)
save_freq_coeff(freq, coeff, "freq_coeff_discrete.txt")

