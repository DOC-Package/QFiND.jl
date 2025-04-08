# BSDO Discretization

## Theory

BSDO method uses Gauss quadrature to numerically integrate the BCF-QNSD relation.  The unique feature of the BSDO method is that it uses the SD as the weight of a quadrature, enabling to construct an efficient set of polynomial interpolants.  In this method, we consider a frequency interval $[\Omega_\mathrm{min},\Omega_\mathrm{max}]$ where the interval covers the entire frequency range of the QNSD i.e. $S_\beta(\omega)$ should be small enough for $\omega<\Omega_\mathrm{min}$ and $\omega>\Omega_\mathrm{max}$.  Note that the choice of $\Omega_\mathrm{min}$ and $\Omega_\mathrm{max}$ affects the accuracy and convergence of the approximation.

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
    \Lambda(z)\approx\sum_{k=1}^M\frac{\lambda}{z-\omega_k},
```
the BCF-QNSD relation can be expressed 
```math
    \begin{aligned}
    C(t) 
    &= -\frac{1}{\pi^2} \operatorname{Im}\int_{-\infty}^\infty \mathrm{d}\omega\; \lim_{\eta\rightarrow 0^+} \Lambda(\omega+i\eta) \mathrm{e}^{-i\omega t}\\
    &=\sum_{k=1}^M \frac{1}{2} g_k^2(\beta)\mathrm{e}^{-i\omega_k t}, \quad \text{with  } g_k^2(\beta) = \frac{2\lambda}{\pi}.
    \end{aligned}
```

To achieve this, we express the integral in terms of the product of a weight function $w(\omega)\geq 0$ and a function $h(\omega,z)$,

```math
\Lambda(z)=\int_{\Omega_\mathrm{mix}}^{\Omega_\mathrm{max}} \mathrm{d} \omega \;\frac{S_\beta(\omega)}{z-\omega}
=\int_{\Omega_\mathrm{mix}}^{\Omega_\mathrm{max}} \mathrm{d} \omega\; w(\omega) h(\omega, z).
```

Now we consider a polynomial interpolant $h_M(\omega,z)$ of $h(\omega,z)$ with degree $M-1$

```math
h(\omega,z) = h_M(\omega,z)+r_M(\omega, z)
```

```math
h_M(\omega,z) =\sum_{k=1}^M h(\omega_k,z) l_k(\omega), \quad l_k(\omega_l)=\delta_{kl}
```

where $l_k(\omega)$ can be defined as the $(M-1)$-th order polynomial

```math
l_k(\omega)=\frac{\prod_{l\neq k}(\omega-\omega_l)}{\prod_{l\neq k}(\omega_k-\omega_l)}
```

and $r_M(\omega,z)$ is a remainder.  
If we choose the $M$ node points $\omega_k$, then

```math
\Lambda(z) =\int_{\Omega_\mathrm{mix}}^{\Omega_\mathrm{max}} \mathrm{d} \omega\; w(\omega) h(\omega, z)
\approx\sum_{k=1}^M W_k h(\omega_k, z)
```

with

```math
W_k =\int_{\Omega_\mathrm{mix}}^{\Omega_\mathrm{max}} \mathrm{d} \omega \; w(\omega) l_k(\omega),
```

which are referred to as Christoffel weights.

In the BSDO approach, the QNSD (originally the SD) is chosen as the weight of the polynomials, i.e.,

```math
w(\omega)=S_\beta(\omega)
```

and

```math
h(\omega,z)=\frac{1}{z-\omega}.
```

The nodes can be determined as the roots of the orthogonal polynomial $p_M(\omega)$ of degree $M$, obeying 

```math
\int_{\Omega_\mathrm{mix}}^{\Omega_\mathrm{max}} \mathrm{d} \omega \;S_\beta(\omega) p_M(\omega) p_l(\omega)=\delta_{kl}.
```

Such polynomials can be generated using the recurrence relation

```math
p_{k+1}(\omega) = (\omega-\alpha_k) p_k(\omega)-\zeta_k p_{k-1}(\omega),\quad k=0, \dots, M-1,
```

where

```math
p_0(\omega) = 1,\quad p_{-1}(\omega)=0.
```

The recurrence coefficients $\alpha_k,\zeta_k$ can be obtained by applying the Lanczos algorithm, [GraggHarrod1984NM,Gautschi1994ATMS,Gautschi2005JoCaAM] as the elements of a tridiagonal matrix

```math
\mathbf{A}_M = \begin{pmatrix}
\alpha_0 & \sqrt{\zeta_1} & 0 & \ldots \\
\sqrt{\zeta_1} & \alpha_1 & \sqrt{\zeta_2} & \ddots \\
0 & \sqrt{\zeta_2} & \alpha_2 & \ddots \\
\vdots & \ddots & \ddots & \ddots
\end{pmatrix} \in \mathbb{R}^{M \times M}.
```

By performing an eigenvalue decomposition,

```math
\mathbf{A}_M = \mathbf{U}\,\mathrm{diag}(\omega_1,\dots,\omega_M)\,\mathbf{U}^{\mathrm{T}},
```

we can obtain the nodes of $p_M$, leading to sample frequencies $\omega_k$, as the eigenvalues. The Christoffel weights are given by the square of the first element of each eigenvector

```math
g_k^2(\beta)= \frac{2\Delta_w}{\pi} \,[\mathbf{U}(1,k)]^2,
```

where

```math
\Delta_w=\int_{\Omega_\mathrm{mix}}^{\Omega_\mathrm{max}} \mathrm{d} \omega\; S_\beta(\omega)
```

is the normalization factor.

