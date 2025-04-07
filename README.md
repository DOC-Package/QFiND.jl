# QFiND

<!--![Static Badge](https://img.shields.io/badge/Version-v0.1.0-brightgreen)!-->
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://doc-package.github.io/QFiND.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://doc-package.github.io/QFiND.jl/dev/)
[![Build Status](https://github.com/DOC-Package/QFiND.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/DOC-Package/QFiND.jl/actions/workflows/CI.yml?query=branch%3Amain)

QFiND (Quantum Finite-temperature Noise Discretizer and Decomposer) provides routines for constructing an effective discrete representation of a system-bath model. In other words, this codes provides an approximation of the bath correlation function $C(t)$ for a given spectral density $J(\omega)$

$$
\begin{aligned}
C(t)&=\frac{1}{2\pi} \int_{-\infty}^{\infty} \mathrm{d}\omega J(\omega)\left[\mathrm{coth}\left(\frac{\beta \omega}{2}\right)+1\right] \mathrm{e}^{-i \omega t}\\
&\approx \sum_{k=1}^M g_k^2 \mathrm{e}^{-i\omega_k t}
\end{aligned}
$$

where $`\omega_k,g_k \in ℝ \backslash \{0\}`$.
The code allows for the estimation of frequencies and coefficients in the system plus bosonic bath model using Interpolative Decomposition (ID) and Non-negative Least Squares (NNLS). 

## Documentation

## Installation

```julia
    julia> using Pkg

    julia> Pkg.add(; url="https://github.com/DOC-package/QFiND.jl")
```

## Simple usage

For other usage, see `./examples`.


## Cite `QFiND`

If you find the framework useful in your research, we would be grateful if you could cite our publications:
- H. Takahashi and R. Borrelli, J. Chem. Phys. 161, 151101 (2024). (https://doi.org/10.1063/5.0232232) 
- H. Takahashi and R. Borrelli, J. Chem. Theory Comput. 21, 2206–2218 (2025). (https://doi.org/10.1021/acs.jctc.4c01728)

Here are the bibtex entries:
```bib
@article{TakahashiBorrelli2024JCP,
  title = {Effective Modeling of Open Quantum Systems by Low-Rank Discretization of Structured Environments},
  author = {Takahashi, Hideaki and Borrelli, Raffaele},
  year = {2024},
  month = oct,
  journal = {The Journal of Chemical Physics},
  volume = {161},
  number = {15},
  pages = {151101},
  issn = {0021-9606},
  doi = {10.1063/5.0232232}
}

@article{TakahashiBorrelli2025JCTC,
  title = {Discretization of {{Structured Bosonic Environments}} at {{Finite Temperature}} by {{Interpolative Decomposition}}: {{Theory}} and {{Application}}},
  shorttitle = {Discretization of {{Structured Bosonic Environments}} at {{Finite Temperature}} by {{Interpolative Decomposition}}},
  author = {Takahashi, Hideaki and Borrelli, Raffaele},
  year = {2025},
  month = mar,
  journal = {Journal of Chemical Theory and Computation},
  volume = {21},
  number = {5},
  pages = {2206--2218},
  publisher = {American Chemical Society},
  issn = {1549-9618},
  doi = {10.1021/acs.jctc.4c01728}
}

```


## Authors

Hideaki Takahashi (hideaki.takahashi@unito.it)


## License

This project is distributed under the [MIT License](./LICENSE).
