# Work Flow

1. Construct a spectral density
As an examples, we consider a Ohmic form with exponential cutoff

```math
J(\omega) = \pi\alpha
```
```julia
PowerLawExpSD(s::Float64, γ::Float64; alpha::Union{Float64,Nothing}=nothing, reorgene::Union{Float64,Nothing}=nothing) 
```
We need to specify either the coupling strength `alpha` or the reorganization energy `reorgene`.  Note that in the Ohmic case ($s=1$), the reorganization energy equals to `alpha`.


```julia
    julia> using QFiND
    
    julia> sdens = PowLowExpSD(1.0, 50.0; alpha=35.0)
    (::PowerLawExpSD) (generic function with 1 method)
```

This `sdens` works as a function.

```julia
    julia> Er = reorganization_energy(sdens)
    35.0
```

2. Construct a QNSD

```julia
    julia> sbeta = BosonicQNSD(sdens, 300.0)
    (::BosonicQNSD) (generic function with 1 method)
```

3. Construct a bath correlation function

```julia
    julia> Temp = 300.0;

    julia> bcf = BosonicBCF(sbeta, Temp)
    (::BosonicBCF) (generic function with 1 method)
```

You can specify the upper bound and threshold of the integral when the next step (4. Create an initial dataset) is very slow, such that

```julia
    julia> bcf = BosonicBCF(sbeta, 300.0)
    (::BosonicBCF) (generic function with 1 method)
```

4. Create an initial dataset
```julia

```