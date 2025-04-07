# Basic Theory

# Linear Dissipation in Thermo-field dynamics

## Linear dissipative model

Let us consider a system linearly coupled to a bosonic heat bath with inverse temperature $\beta$, which consists of a canonical distribution of harmonic oscillators. The Hamiltonian operator can be written as

$$
H = H_\mathrm{S} + \sum_k \omega_{k} a_k^\dagger a_k + V_\mathrm{SB}\sum_{k} \frac{g_k}{\sqrt{2}} \bigl(a_k^\dagger + a_k\bigr).
$$

Here, $H_\mathrm{S}$ is the system Hamiltonian and $V_{\mathrm{SB}}$ is the system part of the system–bath coupling. The symbols $a_{k}$, $a_{k}^\dagger$, $\omega_{k}$, and $g_{k}$ denote the annihilation operator, creation operator, frequency, and system–bath coupling constant for the $k$th mode of the bath, respectively.

From a dynamical point of view, the effect of the bosonic heat bath at an inverse temperature $\beta$ on the system is exclusively determined by the bath correlation function (BCF) [Weiss-book]. One may write the BCF as

$$
C(t)=\frac{1}{\pi} \int_{-\infty}^{\infty} \mathrm{d}\omega \; S_\beta(\omega) e^{-i \omega t},
$$

where the quantum noise spectral density (QNSD) is defined by

$$
S_\beta(\omega)\equiv \frac{1}{2}J(\omega)\Bigl[\coth\Bigl(\frac{\beta \omega}{2}\Bigr)+1\Bigr],
$$

with

$$
J(\omega)=\frac{\pi}{2}\sum_k g_k^2 \Bigl(\delta(\omega-\omega_k)-\delta(\omega+\omega_k)\Bigr).
$$

In the following, we will refer to the above relation as the **BCF-QNSD relation**.

Explicitly writing

$$
\begin{aligned}
C(t) &= \sum_k \frac{1}{2}\Biggl\{ \frac{g_k^2}{2}\Bigl[\coth\Bigl(\frac{\beta \omega_k}{2}\Bigr)+1\Bigr] e^{-i\omega_k t} \\
&\quad\quad - \frac{g_k^2}{2}\Bigl[-\coth\Bigl(\frac{\beta \omega_k}{2}\Bigr)+1\Bigr] e^{i\omega_k t} \Biggr\},
\end{aligned}
$$

we can state that when the BCF can be written in the above form, the corresponding dynamical system is described by the Hamiltonian given earlier and the bath is at the inverse temperature $\beta$.

## Thermo field dynamics

The latter statement can be transformed into a powerful theoretical and computational tool by exploiting the Thermo Field Dynamics (TFD) framework. TFD is a methodology that facilitates the treatment of quantum systems at nonzero temperatures using the wavefunction formalism [UmezawaEtAl1982, Suzuki1991IJMPB, TakahashiUmezawa1996IJMPB]. Here we briefly review the TFD approach and demonstrate its relation to the above BCF-QNSD relation.

We first consider the Hamiltonian of bosonic free particles

$$
H_\mathrm{B} = \sum_{k}\omega_k a_k^\dagger a_k,
$$

and label the eigenstates of the $k$th boson as $\ket{n_k}$. Then, we introduce the so-called *tilde space*, denoted by $\ket{\tilde{n}_k}$, which is the Hilbert space of a fictitious dynamical system identical to the original physical system. Additionally, the tensor product of the physical and tilde spaces is referred to as the *twin space* [Suzuki1991IJMPB, Schmutz1978ZPB]. A ket vector in the twin space is given by

$$
\ket{m_k,\tilde{n}_k} \equiv \ket{m_k}\otimes\ket{\tilde{n}_k},
$$

from which the identity vector is defined as

$$
\ket{I} \equiv \bigotimes_k \sum_{n_k} \ket{n_k,\tilde{n}_k}.
$$

Using the canonical distribution operator

$$
\rho_\mathrm{eq}=\frac{e^{-\beta H_\mathrm{B}}}{\mathrm{Tr}\bigl(e^{-\beta H_\mathrm{B}}\bigr)},
$$

and the definition of $\ket{I}$, the so-called **thermal vacuum state** is derived as

$$
\begin{aligned}
\ket{0(\beta)} &\equiv \rho_\mathrm{eq}^{1/2}\ket{I} \\
&=\prod_k \Bigl[1-e^{-\beta\omega_k}\Bigr]^{1/2} \exp\Bigl(e^{-\beta\omega_k/2}\, a_k^\dagger\tilde{a}_k^\dagger\Bigr)\ket{\mathbf{0}} \\
&= e^{-i G_\theta}\ket{\mathbf{0}},
\end{aligned}
$$

where $\ket{\mathbf{0}}=\bigotimes_{k}\ket{0_k,\tilde{0}_k}$ represents the vacuum state of the ensemble of physical and tilde bosons, and

$$
G_\theta=-i \sum_k \theta_k(\beta)\Bigl(a_k \tilde{a}_k-a_k^{\dagger} \tilde{a}_k^{\dagger}\Bigr),\quad \theta_k(\beta)=\operatorname{arctanh}\Bigl(e^{-\beta \omega_k /2}\Bigr).
$$

The transformation $e^{-i G_\theta}$ is known as the **thermal Bogoliubov transformation**. It is worth noting that the transformation from the exponential form in the second line to the compact form in the third line can be performed using the Baker-Campbell-Hausdorff formulas for the $\mathrm{su}(1,1)$ Lie algebra [Ban1993JOSABJ].

Let us now consider the case in which the initial density matrix of our model can be represented as the direct product of a pure state of the system, $|\psi_\mathrm{e}\rangle\langle\psi_\mathrm{e}|$, and the bath density matrix, i.e.,

$$
\rho(0)=\frac{1}{Z} |\psi_\mathrm{e}\rangle\langle\psi_\mathrm{e}| \, \rho_{\mathrm{B}},
$$

where $Z$ is a properly defined partition function. Under this assumption, it is possible to demonstrate [BorrelliGelin2016JCP, BorrelliGelin2017SR, BorrelliGelin2021WCMS, GelinBorrelli2017AdP, deVegaBanuls2015PRA] that the expectation value of an arbitrary operator $A$ acting on the Hilbert space of the physical system can be written as

$$
\langle A(t) \rangle = \langle A_\theta \rangle_{\psi_\theta(t)},
$$

where the wavefunction $\ket{\psi_\theta(t)}$ satisfies the Schrödinger equation

$$
\begin{aligned}
i \frac{\mathrm{d}}{\mathrm{d}t}\ket{\psi_\theta(t)} &= H_\theta\ket{\psi_\theta(t)},\\[1mm]
\ket{\psi_\theta(0)} &= \ket{\psi_\mathrm{e}}\otimes\ket{\mathbf{0}},
\end{aligned}
$$

with the thermal operators defined as

$$
H_\theta=\; e^{i G_\theta} \hat{H}\, e^{-i G_\theta},\qquad A_\theta=\; e^{i G_\theta} A\, e^{-i G_\theta}.
$$

The modified Hamiltonian operator $\hat{H}$ is defined by

$$
\hat{H}=H-\tilde{H}_\mathrm{B},
$$

where $\tilde{H}_{\mathrm{B}} = \sum_k \omega_k \tilde{a}_k^{\dagger} \tilde{a}_k$ is the free-boson Hamiltonian operator of the bosonic tilde space. The thermal Hamiltonian $H_\theta$, which governs the finite-temperature dynamics, is obtained by applying the Bogoliubov thermal transformation to $\hat{H}$:

$$
\begin{aligned}
H_\theta =\; e^{iG_\theta} \hat{H}\, e^{-iG_\theta} 
=&\; H_\mathrm{S}+\sum_k \omega_k\Bigl(a_k^{\dagger} a_k-\tilde{a}_k^{\dagger} \tilde{a}_k\Bigr) \\
&+\; V_\mathrm{SB}\sum_{k} \frac{g_{k}}{\sqrt{2}} \Bigl\{\Bigl(a_{k}+a_{k}^{\dagger}\Bigr) \cosh \bigl(\theta_k\bigr) \\
&\quad\quad + \Bigl(\tilde{a}_{k}+\tilde{a}_{k}^{\dagger}\Bigr) \sinh\bigl(\theta_k\bigr)\Bigr\}.
\end{aligned}
$$

Representing the operators and frequencies of the physical and tilde spaces with common symbols as

- $\{a_k,\tilde{a}_k\} \rightarrow \{a_k\}$,
- $\{\omega_k,-\omega_k\} \rightarrow \{\omega_k\}$,
- $\{g_k\cosh{\theta_k},\,g_k\sinh{\theta_k}\} \rightarrow \{g_k(\beta)\}$,

the TFD Hamiltonian can be expressed in the simpler form

$$
H_\theta =  H_\mathrm{S}+\sum_k \omega_k a_k^{\dagger} a_k +V_\mathrm{SB}\sum_{k} \frac{g_k(\beta)}{\sqrt{2}}\bigl(a_{k}+a_{k}^{\dagger}\bigr).
$$

## Connection between TFD and the BCF-QNSD relation

The parameters $g_{k} \cosh (\theta_k)$ and $g_{k} \sinh (\theta_k)$ entering the thermal Hamiltonian $H_\theta$ govern the coupling of the subsystem with the physical and tilde bosonic degrees of freedom. The spectral density (SD) for the thermal system–bath Hamiltonian can be written as

$$
J_\theta(\omega) = J_p(\omega) + J_t(\omega),
$$

where

$$
J_p(\omega) \equiv \frac{\pi}{2}\sum_k \Bigl(g_{k} \cosh \bigl(\theta_k\bigr)\Bigr)^2 \delta\bigl(\omega-\omega_k\bigr)\, \Theta(\omega),
$$

$$
J_t(\omega) \equiv \frac{\pi}{2}\sum_k \Bigl(g_{k} \sinh \bigl(\theta_k\bigr)\Bigr)^2 \delta\bigl(\omega+\omega_k\bigr)\, \Theta(-\omega),
$$

and $\Theta(\omega)$ is the Heaviside step function. These functions describe the system–bath couplings in the physical (subscript $p$) and tilde (subscript $t$) subspaces, respectively. As the temperature approaches zero, \$J_p(\omega) \rightarrow J(\omega)\$ and $J_t(\omega) \rightarrow 0$.

By using the relations

$$
\cosh^2(\theta_k) = \frac{1}{2}\Bigl[\coth\Bigl(\frac{\beta\omega_k}{2}\Bigr)+1\Bigr],
$$

$$
\sinh^2(\theta_k) = -\frac{1}{2}\Bigl[-\coth\Bigl(\frac{\beta\omega_k}{2}\Bigr)+1\Bigr],
$$

one immediately sees that the definitions of $S_\beta(\omega)$ and $J_\theta(\omega)$ are equivalent, namely

$$
J_\theta(\omega) = S_\beta(\omega).
$$

Therefore, we have demonstrated that the system–bath Hamiltonian at temperature $\beta$ and the extended zero-temperature Hamiltonian have the same correlation function and hence the same dynamical behavior.

In conclusion, whenever the BCF can be written in the form

$$
C(t) = \sum_k \frac{1}{2} g_k^2(\beta) e^{-i\omega_k t},\quad \omega_k,\,g_k(\beta)\in\mathbb{R},
$$

the system dynamics at temperature $\beta$ is described by the TFD Hamiltonian given above. Note that in the above equation the frequencies $\omega_k$ can be negative. Unlike the standard TFD formalism [deVegaBanuls2015PRA, BorrelliGelin2016JCP], we do not perform thermo-field doubling but instead directly determine a set of temperature-dependent parameters $(\omega_k, g_k(\beta))$.
