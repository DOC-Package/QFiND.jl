using QFiND
  
omega = 100.0
g = 300.0
Temp = 100.0
freq, coeff = ThermalBogoliubov(omega, g, Temp)
save_freq_coeff(freq, coeff, "freq_coeff_discrete.txt")

expon, coeff = BosonicQNSD_Discrete(omega, g, Temp)
save_expon_coeff_union(expon, coeff, "expon_coeff_discrete.txt")

omega = [100.0, 50.0]
g = [300.0, 150.0]
Temp = 100.0
freq, coeff = ThermalBogoliubov(omega, g, Temp)
save_freq_coeff(freq, coeff, "freq_coeff_discrete2.txt")
