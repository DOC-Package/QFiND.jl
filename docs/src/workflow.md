# Work Flow

1. Define a spectral density
For examples, we consider a Ohmic form with exponential cutoff.
```math
J(\omega) = \pi\alpha
```

```julia
sdens = PowLowExpSD(1.0, 50.0, alpha=35.0)
```

This `sdens` works as a function.

2. Construct a QNSD

```julia
sbeta = QuantumNoiseSD(sdens, 300.0)
```

3. Bath Correlation Function

```julia
bcf = BosonicBCF(sbeta)
```