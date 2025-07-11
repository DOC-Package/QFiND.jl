# Decomposition Methods

```math
C(t)=\sum_k c_k \mathrm{e}^{-\gamma_k t}
```

## Exponential Fitting to BCF

ExpFit.jl library, which is a part of our package.

## Truncated Padé Spectral Decomposition

Padé spectral decomposition (PSD) is often employed to decompose Bose-Einstein distribution into a sum of Lorentzian functions

```math
\coth(\frac{\beta\hbar\omega}{2}) + 1 = \frac{2}{\beta\hbar\omega} + 1 + \frac{2}{\beta\hbar}\sum_{k=1}^N \frac{2\eta_k\omega}{\omega^2+\nu_k^2}
```

