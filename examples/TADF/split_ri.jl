using LinearAlgebra
function split_real_imag(freq, coeff)
    N = length(freq)
    freq_real = zeros(ComplexF64, 2*N)
    coeff_real = zeros(ComplexF64, 2*N)
    freq_imag = zeros(ComplexF64, 2*N)
    coeff_imag = zeros(ComplexF64, 2*N)
    for i in 1:N
        gam = freq[i]
        c = coeff[i] 
        # store results
        freq_real[2*i-1] = gam
        freq_real[2*i]   = conj(gam)
        freq_imag[2*i-1] = gam
        freq_imag[2*i]   = conj(gam)
        coeff_real[2*i-1] = c / 2.0
        coeff_real[2*i]   = conj(c) / 2.0
        coeff_imag[2*i-1] = c / 2.0 /im
        coeff_imag[2*i]   = - conj(c) / 2.0 / im
    end

    return freq_real, coeff_real, freq_imag, coeff_imag
end