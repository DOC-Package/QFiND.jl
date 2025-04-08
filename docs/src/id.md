# ID Discretization

## Theory

Our goal is finding an approximation $\tilde C(t)$ of the BCF of the form of 
```math
\tilde C(t) = \sum_k \frac{1}{2} g_k^2(\beta) e^{-i\omega_k t},\quad \omega_k,\,g_k(\beta)\in\mathbb{R},
```
within a specified time interval $[0,T]$, having a small number of summands, $M$. That is, the problem we want to solve can be cast in the following form

```math
\mathrm{Given: } 
C(t)=\frac{1}{\pi} \int_{-\infty}^{\infty} \mathrm{d} \omega \;S_\beta(\omega) \mathrm{e}^{-i \omega t},\quad t \in[0, T],
```

```math
\text{Find: }\min M,\ \{\omega_k, g_k\}_{k=1}^M \subset \mathbb{R} \times \mathbb{R},\quad \text{such that: }
```

```math
\left\|C(t)-\sum_{k=1}^M \frac{1}{2}g_k^2 \mathrm{e}^{-i \omega_k t}\right\|_2 <\delta,
```

where $\delta$ is a threshold that depends on the required accuracy level.

Rather than computing $C(t)$ and directly solving the above minimization problem, we first determine the frequencies $\omega_k$ by exploiting the low-rank structure using the ID theory. Subsequently, we obtain the coefficients $g_k$ by the so-called nonnegative least squares (NNLS) technique.

In the following, we present the fundamental theoretical framework and the detailed procedure of the methodology.

### Calculating the frequencies $\omega_k$

#### Initial discretization of $f(t,\omega)$

We begin by representing the BCF-QNSD relation in a form that is suitable for numerical analysis. First, we set a cutoff frequency $\Omega$ such that $S_\beta(\omega) = 0$ outside the interval $[-\Omega,\Omega]$ and a cutoff time $T$. Then, we rewrite the BCF-QNSD relation as

```math
C(t)=\int_{-\Omega}^{\Omega} \mathrm{d} \omega \;f(t,\omega),\quad t\in [0,T],
```

where
```math
f(t,\omega)\equiv \frac{1}{\pi}S_\beta(\omega) \mathrm{e}^{-i \omega t}.
```

We now define the sets  
$$\mathbf{T}=\{t_i\}_{i=1}^{m}$$  
and  
$$\boldsymbol{\Omega}=\{\omega_j\}_{j=1}^{n},$$  
where $t_i \in [0,T]$ and $\omega_j \in [-\Omega,\Omega]$, and create a dense rectangular grid $\mathbf{T}\times\boldsymbol{\Omega}$ on which $f(t,\omega)$ is discretized. In this study, we employ an equispaced fine grid.
Then the integral can be approximated for each $t_i$ using a quadrature

```math
C(t_i) = \sum_{j=1}^n f(t_i,\omega_j)w_j + e_i,
```

where $w_j$ are the weights of a quadrature and $e_i$ are approximation errors.

Using the information contained in the matrix $f(t_i,\omega_j)$, we obtain the subset of frequencies $\omega_k$ defining our discretized model, and use the actual values $C(t_i)$ to determine the system-bath coupling strength. The procedure for obtaining the frequency subset is based on the ID (Interpolative Decomposition) theory.

#### Interpolative Decomposition

In order to apply ID to our problem, we construct a real matrix $\mathbf{f}_{2m\times n}$ with entries $\operatorname{Re} f(t_i,\omega_j)$ and $\operatorname{Im} f(t_i,\omega_j)$:

```math
\mathbf{f}_{2m\times n} = \begin{pmatrix}
(\operatorname{Re}f(t_i,\omega_j))_{i,j=1}^{m,n} \\[5pt]
(\operatorname{Im}f(t_i,\omega_j))_{i,j=1}^{m,n}
\end{pmatrix} \in\mathbb{R}^{2m\times n}.
```

Then we can construct a rank $r$ ID of $\mathbf{f}_{2m\times n}$ as

```math
\mathbf{f}_{2m \times n} = \mathbf{B}_{2m\times r}\mathbf{P}_{r\times n} + \mathbf{E}_{2m\times n},
```

where

```math
\mathbf{B}_{2m\times r} = \begin{pmatrix}
(\operatorname{Re}f(t_i,\omega_k))_{i,k=1}^{m,r} \\[5pt]
(\operatorname{Im}f(t_i,\omega_k))_{i,k=1}^{m,r}
\end{pmatrix} \in\mathbb{R}^{2m\times r},
```

with $\omega_k \in \tilde{\boldsymbol{\Omega}} \subset \boldsymbol{\Omega}$, $\mathbf{P}_{r\times n}\in\mathbb{R}^{r\times n}$, and $\mathbf{E}_{2m \times n}\in\mathbb{R}^{2m \times n}$ is an error matrix satisfying $\|\mathbf{E}_{2m \times n}\|_2\leq\epsilon$. 

The subset $\tilde{\boldsymbol{\Omega}}$ corresponds to the selected columns in the ID. Entrywise, $f(t_i,\omega_j)$ can be written as

```math
f(t_i,\omega_j)=\sum_{k=1}^r f(t_i,\omega_k)P_{kj}+E_{ij}',
```

where $P_{kj}=(\mathbf{P}_{r\times n})_{kj}$ and $E_{ij}'$ represents the corresponding error.

Substituting this expression into the quadrature formula for $C(t_i)$, we obtain

```math
C(t_i)=\sum_{j=1}^n\sum_{k=1}^r f(t_i,\omega_k) P_{kj} w_j + \sum_{j=1}^n E'_{ij}w_j  + e_i.
```

This can be rearranged as

```math
C(t_i)= \sum_{k=1}^r f(t_i,\omega_k) z_k  + e'_i,
```

where 

```math
z_k = \sum_{j=1}^n P_{kj}w_j,\quad e'_i = e_i+\sum_{j=1}^n E'_{ij}w_j.
```

At this stage, a specific subset of $r$ frequencies $\omega_k$ has been selected from the initial set $\boldsymbol{\Omega}$ to represent the discretized bath. The final step involves determining the system-bath coupling parameters $g_k(\beta)$ by evaluating the weights $z_k$.

### Calculating the coupling strengths $g_k(\beta)$

Let us now define the vector $\mathbf{c}$ that stores the values of the BCF for each sample time

```math
\mathbf{c}=\begin{pmatrix}
(\operatorname{Re}C(t_{i}))_{i=1}^{m} \\[5pt]
(\operatorname{Im}C(t_{i}))_{i=1}^{m}
\end{pmatrix}  \in\mathbb{R}^{2m}.
```

Then, from the previous expression for $C(t_i)$, the coefficients $\mathbf{z}=(z_k)_{k=1}^{r}\in \mathbb{R}_{+}^{r}$ can be obtained by solving the overdetermined linear system

```math
\min_\mathbf{z}\left\|\mathbf{c}-\mathbf{B}\mathbf{z}\right\|_2^2 \quad\text{s.t.} \quad z_k \geq 0.
```

This problem can be solved using the NNLS technique\cite{Lawson1974}. Finally, we obtain the desired approximation

```math
\tilde{C}(t)=\sum_{k=1}^{M} f(t,\omega_k)z_k,\quad z_k > 0,
```

or equivalently,

```math
\tilde{C}(t)=\sum_{k=1}^{M} z_k S_\beta(\omega_k) \mathrm{e}^{-i \omega_k t},\quad z_k > 0.
```

For the Hamiltonian operator to be Hermitian, the weights $z_k$ must be positive. From the approximated BCF above, we can construct a temperature-dependent effective Hamiltonian (see Eq. \eqref{eq:tfdgenericham}) by setting

```math
g_k(\beta)=\sqrt{\frac{2}{\pi}z_k S_\beta(\omega_k)}.
```

Since the NNLS procedure may result in some $z_k$ being zero, in the general case $M\leq r$.

The overall accuracy of the method is determined by the accuracy of the ID and the NNLS. Note that ID induces a pivoting of the columns of $\mathbf{f}$; therefore, the $M$ frequencies $\{\omega_k\}$ do not correspond to the first $M$ frequencies of the initial discretization grid, but rather form a special subset of it.

