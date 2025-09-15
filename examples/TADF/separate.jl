using LinearAlgebra
function separate(freq, coeff)
    N = length(freq)
    freq_real = zeros(ComplexF64, 2*N)
    freq_imag = zeros(ComplexF64, 2*N)
    coeff_real = zeros(ComplexF64, 2*N)
    coeff_imag = zeros(ComplexF64, 2*N)
    for i in 1:N
        gam = real(freq[i])
        omg = imag(freq[i])
        a = real(coeff[i])      
        b = imag(coeff[i])
        G = [-gam omg; -omg -gam]
        
        # Gの行列のフロベニウスノルム（元のノルム）を計算
        G_norm = norm(G)
        println("G行列のノルム: ", G_norm)
        
        # eigen decomposition
        result = eigen(G)
        ev = result.values
        U = result.vectors
        
        println("固有値: ", ev)
        println("固有ベクトル: ", U)
        
        # store results
        freq_real[2*i-1] = -ev[1]
        freq_real[2*i]   = -ev[2]
        freq_imag[2*i-1] = -ev[1]
        freq_imag[2*i]   = -ev[2]
        
        # 固有ベクトルにGのノルムを使って補正
        Uinv = inv(U)
        coeff_real[2*i-1] = U[1, 1] * G_norm
        coeff_real[2*i]   = U[1, 2] * G_norm
        coeff_imag[2*i-1] = U[2, 1] * G_norm
        coeff_imag[2*i]   = U[2, 2] * G_norm
    end

    return freq_real, coeff_real, freq_imag, coeff_imag
end