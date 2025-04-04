using Printf

function save_freq_coeff(freq::Vector{Float64}, coeff::Vector{Float64}, filename::String) 
    open(filename, "w") do io
        write(io, " "^4 * "Frequencies" * " "^10 * "Coefficients\n")
        write(io, "-"^40 * "\n")  
        
        for (f, c) in zip(freq, coeff)
            write(io, @sprintf("%0.12e    %0.12e\n", 
                               f, c))
        end
    end
end