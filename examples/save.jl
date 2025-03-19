using Printf

function save_array(filename::String, array1::Vector{Complex{T}}, array2::Vector{Complex{T}}) where T
    open(filename, "w") do io
        write(io, " "^15 * "Coefficients" * " "^25 * "Exponents\n")
        write(io, " "^10 * "Real" * " "^10 * "Imag" * " "^20 * "Real" * " "^10 * "Imag\n")
        write(io, "-"^80 * "\n")  
        
        for (c1, c2) in zip(array1, array2)
            write(io, @sprintf("%0.12e  %0.12e  %0.12e  %0.12e\n", 
                               real(c1), imag(c1), real(c2), imag(c2)))
            write(io, @sprintf("%0.12e  %0.12e  %0.12e  %0.12e\n", 
                               real(c1), -imag(c1), real(c2), -imag(c2)))
        end
    end
end

function save_array_union(filename::String, array1::Vector{Complex{T}}, array2::Vector{Complex{T}}) where T
    open(filename, "w") do io
        for (c1, c2) in zip(array1, array2)
            write(io, @sprintf("%0.12e  %0.12e  %0.12e  %0.12e  %0.12e  %0.12e\n", 
                               real(c1), imag(c1), 0.0, 0.0, real(c2), imag(c2)))
            write(io, @sprintf("%0.12e  %0.12e  %0.12e  %0.12e  %0.12e  %0.12e\n", 
                               0.0, 0.0, real(c1), -imag(c1), real(c2), -imag(c2)))
        end
    end
end