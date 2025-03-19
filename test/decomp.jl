```
function bath_correlation_exponentials(sd::SpectralDensity, beta::Real; npade::Int=5)
    exp_list = ComplexF64[]
    coeff_list = ComplexF64[]
    
    # 1. スペクトル密度 J(ω) の極からの寄与
    #    ここでは sd に対して nodes(sd) と residues(sd) を用いる．
    nodesJ = nodes(sd)
    resJ = residues(sd)
    @assert length(nodesJ) == length(resJ) "J(ω)の極と留数の個数が一致していません"
    for (ω_p, r_p) in zip(nodesJ, resJ)
        # 対応する温度因子は coth(βω_p/2)+1
        factor = coth((beta*ω_p)/2) + 1
        # 寄与係数
        A = r_p * factor
        # t>0 で下半平面の極のみ寄与（実際は contour 選択による）
        if imag(ω_p) < 0
            push!(exp_list, ω_p)
            push!(coeff_list, A)
        end
    end
    
    # 2. Pade展開による温度因子の寄与
    #    coth(βω/2) ≈ 1/(βω/2) + Σ_{j=1}^{npade} [2η_j (βω/2)]/((βω/2)²+ξ_j²)
    #    から、分数項は ω の単純極 ω_j = -2i ξ_j/β を持つ．
    #    その寄与は、極 ω_j における J(ω_j) によるものとし，
    #    留数計算の結果、係数は (2η_j/β) J(ω_j) となるとする．
    pade_poles, η = pade_coth_terms(beta, npade)
    for (ω_pade, η_j) in zip(pade_poles, η)
        A = (2η_j/ beta) * sd(ω_pade)
        push!(exp_list, ω_pade)
        push!(coeff_list, A)
    end

    return exp_list, coeff_list
end


function spectral_decomposition(sdens::Function, Temp::Real, 
    pole_sd, res_sd, pole_be, res_be, f_func)

    in_contour(ω) = (imag(ω) < 0) 

    coeff = ComplexF64[]
    expon = ComplexF64[]

    for (ω, r) in zip(pole_sd, res_sd)
        if in_contour(ω)
            factor = coth((beta*ω)/2) + 1
            A = r * factor
            push!(expon, ω)
            push!(coeff, 2π * im * A)
        end
    end

    pade_poles, η = pade_coth_terms(beta, npade)
    for (ω, η_j) in zip(pade_poles, η)
        if in_contour(ω)
            A = (2η_j/ beta) * sdens(ω)
            push!(expon, ω)
            push!(coeff, 2π * im * A)
        end
    end

    return (coeff=coeff, expon=expon)
end
```