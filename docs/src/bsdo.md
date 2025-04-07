# Bath Spectral Density Orthogonal 

## Theory

BSDO method uses Gauss quadrature to numerically integrate Eq. \eqref{eq:fdt1}.  The unique feature of the BSDO method is that it uses the SD as the weight of a quadrature, enabling to construct an efficient set of polynomial interpolants.  In this method, we consider a more general frequency interval $[\Omega_\mathrm{min},\Omega_\mathrm{max}]$.  Note that the choice of $\Omega_\mathrm{min}$ and $\Omega_\mathrm{max}$ affects the accuracy and convergence of the approximation, and it is desirable that the interval covers the entire frequency range of the QNSD i.e. $S_\beta(\omega)$ should be small enough for $\omega<\Omega_\mathrm{min}$ and $\omega>\Omega_\mathrm{max}$.

We first introduce the so-called hybridization function
```math
    \Lambda(z)=\int_{\Omega_\mathrm{min}}^{\Omega_\mathrm{max}} \mathrm{d} \omega \;\frac{S_\beta(\omega)}{z-\omega}, \quad z \in \mathbb{C}.
```
By the Sokhotski-Plemelj theorem, Eq. \eqref{eq:lambda} implies
```math
    S_\beta(\omega)=-\frac{1}{\pi} \lim_{\eta\rightarrow 0^+}\operatorname{Im} \Lambda(\omega+i \eta).
```
If $\Lambda(z)$ is approximated with a sum of rational functions as
```math
    \Lambda(z)\approx\sum_{k=1}^M\frac{g_k^2(\beta)}{z-\omega_k},
```
the BCF-QNSD relation can be expressed 
```math
    \begin{aligned}
    C(t) 
    &= -\frac{1}{\pi} \operatorname{Im}\int_{-\infty}^\infty \mathrm{d}\omega\; \lim_{\eta\rightarrow 0^+} \Lambda(\omega+i\eta) \mathrm{e}^{-i\omega t}\\
    &=\sum_{k=1}^M g_k^2(\beta)\mathrm{e}^{-i\omega_k t}.
    \end{aligned}
```
Thus we consider expressing $\Lambda(z)$ in the form of Eq. \eqref{eq:lambda_disc} using a Gauss quadrature.

## `bsdo_discr`

