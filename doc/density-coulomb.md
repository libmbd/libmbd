We consider a system of $N$ distinguishable particles with masses $m_A$, each individually confined in a harmonic potential with frequency $\omega_A$ and center $\mathbf R_A$. The distinguishability is mandated by the individual (rather than global) harmonic potentials. Each particle has a charge $q_A$, which is compensated by the charge of an opposite sign located at $\mathbf R_A$. The compensating charges can be thought of as nuclei of infinite mass, and the harmonic potentials as representing a harmonic force between the particles and the nuclei. Physical motivation fur such a system is that each particle simulates all electrons of a given atom; the harmonic potential is a Taylor expansion of the true pseudopotential for the valence electrons.

The Hamiltonian of such a system is
$$
\begin{multline}
H=\sum_A\left(-\frac1{2m_A}\boldsymbol\nabla_A^2+\frac12m_A\omega_A^2\lvert\mathbf r_A-\mathbf R_A\rvert^2\right)\\
+\frac12\sum_{AB}q_Aq_B\left(\frac1{\lvert\mathbf r_A-\mathbf r_B\rvert}-\frac1{\lvert\mathbf R_A-\mathbf r_B\rvert}-\frac1{\lvert\mathbf r_A-\mathbf R_B\rvert}+\frac1{\lvert\mathbf R_A-\mathbf R_B\rvert}\right)
\label{eq:Hcoulomb}
\end{multline}
$$
where $\mathbf r_A$ is a position of the $A$-th particle. Assuming small displacements, $\mathbf r_A-\mathbf R_A$, of the particles with respect to the distances between them, $\mathbf R_A-\mathbf R_B$, we can expand the electrostatic terms around the equilibrium positions. At first order, this leads to the dipole potential:
$$
\begin{equation}
\begin{aligned}
\mathbf T_{AB}&=\boldsymbol\nabla_A\otimes\boldsymbol\nabla_B\frac1{|\mathbf r_A-\mathbf r_B|}\Bigg|_{\substack{\mathbf r_A=\mathbf R_A\\\mathbf r_B=\mathbf R_B}} \\
&=\frac{\lvert\mathbf R_A-\mathbf R_B\rvert^2-3(\mathbf R_A-\mathbf R_B)\otimes(\mathbf R_A-\mathbf R_B)}{\lvert\mathbf R_A-\mathbf R_B\rvert^5}
\end{aligned}
\label{eq:dip}
\end{equation}
$$
Inserting this expression into $\eqref{eq:Hcoulomb}$, we get
$$
\begin{multline}
H=\sum_A\left(-\frac1{2m_A}\boldsymbol\nabla_A^2+\frac12m_A\omega_A^2\lvert\mathbf r_A-\mathbf R_A\rvert^2\right)\\
+\frac12\sum_{AB}q_Aq_B(\mathbf r_A-\mathbf R_A)\cdot\mathbf T_{AB}\cdot(\mathbf r_B-\mathbf R_B)
\label{eq:MBDhamil}
\end{multline}
$$
Next, we define $3N$ generalized coordinates $\xi_{3(A-1)+\alpha}=\sqrt{m_A}(r_A^\alpha-R_A^\alpha)$, where $\alpha$ labels the Cartesian dimensions. This also defines a general correspondence between the particle index $A$ and coordinate index $i$, so, for instance, when we write $m_i$, we mean the mass of $A$-th particle that corresponds to the $i$-th generalized coordinate.

Using the relationship between the static polarizability $\alpha_0$, mass, frequency, and charge of a charged particle in a harmonic potential, $\alpha_0=q^2/m\omega^2$, and the generalized coordinates, the Hamiltonian can be rewritten as
$$
\begin{equation}
\begin{aligned}
H&=\sum_{i=1}^{3N}\left(-\frac12\frac{\partial^2}{\partial\xi_i^2}+\frac12\omega_i^2\xi_i^2\right)+\frac12\sum_{ij}\omega_i\omega_j\sqrt{\alpha_{0,i}\alpha_{0,j}}\xi_iT_{ij}\xi_j \\
&=\left(\sum_i-\frac12\frac{\partial^2}{\partial\xi_i^2}\right)+\frac12\boldsymbol\xi^\mathrm T\mathbf D\boldsymbol\xi, \qquad D_{ij}=\omega_i^2\delta_{ij}+\omega_i\omega_j\sqrt{\alpha_{0,i}\alpha_{0,j}}T_{ij}
\end{aligned}
\end{equation}
$$
Now, we search for a unitary transformation of $\boldsymbol\xi$ such that the dipole interaction term disappears in the new transformed coordinates. Since the kinetic term ($3N$-dimensional Laplacian) is invariant with respect to unitary transformations, we obtain the new coordinates by diagonalizing $\mathbf D$:
$$
H=\sum_i\left(-\frac12\frac{\partial^2}{\partial\tilde\xi_i^2}+\frac12\tilde\omega_i^2\tilde\xi_i^2,\right)
$$
Here, $\tilde\omega_i^2$ are the eigenvalues of $\mathbf D$ and $\tilde\xi_i$ are the transformed coordinates, $\tilde\xi_i=\sum_kC_{ki}\xi_k$, $\mathbf C$ being the column-wise matrix of the eigenvectors of $\mathbf D$.

But this is a Hamiltonian describing $3N$ uncoupled harmonic oscillators with frequencies $\tilde\omega_i$. The ground-state energy is
$$
E=\sum_i\frac{\tilde\omega_i}2
$$
The wavefunction is
$$
\Psi(\tilde\xi_i)=\prod_i\left(\frac{\tilde\omega_i}{\pi}\right)^\frac14\exp\left(-\frac12\tilde\omega_i\tilde\xi_i^2\right)
$$

### Particle density

Next, we want to calculate the particle density
$$
n(\mathbf r)=\langle\Psi|\sum_Aq_A\delta(\mathbf r-\mathbf r_A)|\Psi\rangle=\idotsint\mathrm d\mathbf r_1\cdots\mathrm d\mathbf r_N\sum_Aq_A\delta(\mathbf r-\mathbf r_A)\Psi(\mathbf r_B)^2
\label{eq:density}
$$
First, we transform the wavefunction back to $\boldsymbol\xi$ and gather the product of the exponentials,
$$
\Psi(\xi_i)=\left[\prod_i\left(\frac{\tilde\omega_i}{\pi}\right)^\frac14\right]\exp\Bigg(-\frac12\sum_{jk}\underbrace{\sum_iC_{ji}\tilde\omega_iC_{ki}}_{\Omega_{jk}}\xi_j\xi_k\Bigg)
\label{eq:wavefnc}
$$
In the following, we will use $\sum_{i\notin A}$ for a sum that skips the $A$-th particle and $\sum_{i\in A}$ for a sum over the three Cartesian coordinates of the $A$-th particle.

For a given $A$, we divide the sum over $jk$ according to the order of $\xi_{i\in A}$:
$$
\begin{equation}
\begin{aligned}
\sum_{jk}\Omega_{jk}\xi_j\xi_k&=\sum_{\substack{j\notin A\\k\notin A}}\Omega_{jk}\xi_j\xi_k+2\sum_{\substack{p\in A\\k\notin A}}\Omega_{pk}\xi_p\xi_k+2\sum_{\substack{p\in A\\q\in A}}\Omega_{pq}\xi_p\xi_q \\
&\equiv\boldsymbol\xi'^\mathrm T_A\boldsymbol\Omega''_A\boldsymbol\xi'_A+2\boldsymbol\xi_A'^\mathrm T\boldsymbol\Omega'_A\boldsymbol\xi_A+\boldsymbol\xi_A^\mathrm T\boldsymbol\Omega_A\boldsymbol\xi_A
\end{aligned}
\end{equation}
$$
[Completing the square](https://en.wikipedia.org/wiki/Completing_the_square#Formula) with respect to $\boldsymbol\xi'_A$, we get
$$
\sum_{jk}\Omega_{jk}\xi_j\xi_k=(\boldsymbol\xi'^\mathrm T_A-\mathbf h_A^\mathrm T)\boldsymbol\Omega''_A(\boldsymbol\xi'_A-\mathbf h_A)+\boldsymbol\xi_A^\mathrm T\boldsymbol\Omega_A\boldsymbol\xi_A-\boldsymbol\xi_A^\mathrm T\boldsymbol\Omega'^\mathrm T_A\boldsymbol\Omega_A''^{-1}\boldsymbol\Omega_A'\boldsymbol\xi_A
\label{eq:completesq}
$$
where $\mathbf h_A$ is some quantity that does not depend on $\boldsymbol\xi'_A$. We can now factor out the exponential and the $3N$-dimensional integral:
$$
n(\mathbf r)=\sum_Aq_A\left(\idotsint\mathrm d\mathbf r_1\cdots\mathrm d\mathbf r_{A-1}\mathrm d\mathbf r_{A+1}\mathrm d\cdots\mathbf r_N\right)\int\mathrm d\mathbf r_A\delta(\mathbf r-\mathbf r_A)\ldots
$$
First, we deal with the integrals in parentheses. Because $\mathbf h_A$ is constant there, and the integrals are over the whole space, $\mathbf h_A$ can be transformed away. Furthermore, we can rotate $\boldsymbol\Omega''_A$ into a new basis where it becomes diagonal, which factors the $3(N-1)$-dimensional integral into a product of $3(N-1)$ 1-dimensional integrals over gaussian functions of the form $\exp(-\bar\omega_{A,i}\bar\xi_{A,i}^2)$, where $\bar\omega_{A,i}$ are the eigenvalues of $\boldsymbol\Omega''_A$. (The factor of $\frac12$ disappears due to the square of the wavefunction.)

Second, the integral over $\mathbf r_A$ picks the value of the following function at point $\mathbf r$ via the $\delta$-function:
$$
\exp\big(-\boldsymbol\xi_A^\mathrm T(\underbrace{\boldsymbol\Omega_A-\boldsymbol\Omega'^\mathrm T_A\boldsymbol\Omega_A''^{-1}\boldsymbol\Omega_A'}_{\boldsymbol\Omega^{(A)}})\boldsymbol\xi_A\big)
$$
Combining $\eqref{eq:density}$, $\eqref{eq:wavefnc}$, $\eqref{eq:completesq}$, and the previous two paragraphs, and transforming from $\boldsymbol\xi_A$ back to $\mathbf r_A$, we get
$$
n(\mathbf r)=\sum_Aq_A\left(\frac{m_A}\pi\right)^\frac32\sqrt{\frac{\prod_{i=1}^{3N}\tilde\omega_i}{\prod_{i=1}^{3(N-1)}\bar\omega_{A,i}}}\exp\big(-m_A(\mathbf r-\mathbf R_A)^\mathrm T\boldsymbol\Omega^{(A)}(\mathbf r-\mathbf R_A)\big)
$$

### Coulomb interaction

We want to calculate a first-order perturbation correction to the dipole approximation, that is,
$$
E^{(1)}=\langle\Psi|V_\text{ee}-V_{\mathbf{pp}}|\Psi\rangle
$$
where $V_\text{ee}$ is the full Coulomb interaction and $V_\mathbf{pp}$ is the dipole interaction.

We start by calculating
$$
\langle\Psi|\frac12\sum_{AB}\frac{q_Aq_B}{\lvert\mathbf r_A-\mathbf r_B\rvert}|\Psi\rangle
$$
In analogy to the calculation of $n(\mathbf r)$, we rotate the $3(N-2)$ coordinates that do not participate in the Coulomb integral such that the integrals become integrals over gaussian functions, and then evaluate the remaining 6-dimensional integral over $\mathbf r_A$ and $\mathbf r_B$. To this end, we need to evaluate
$$
I=\iint\mathrm d\boldsymbol\xi_A\mathrm d\boldsymbol\xi_B \frac{\exp\big(-\boldsymbol\xi_{AB}^\mathrm T\boldsymbol\Omega^{(AB)}\boldsymbol\xi_{AB}\big)}{\lvert\mathbf r_A-\mathbf r_B\rvert}
\label{eq:intcoulomb}
$$
where $\boldsymbol\xi_{AB}$ is a 6-dimensional vector containing $\boldsymbol\xi_A$ and $\boldsymbol\xi_B$, and $\boldsymbol\Omega^{(AB)}$ is the equivalent of $\boldsymbol\Omega^{(A)}$ from the previous section.

We start by rewriting the Coulomb potential as
$$
\frac1{\lvert\mathbf r_A-\mathbf r_B\rvert}=\frac2{\sqrt\pi}\int_0^\infty\mathrm du\exp(-\lvert\mathbf r_A-\mathbf r_B\rvert^2u^2)
$$
Inserting into $\eqref{eq:intcoulomb}$, and transforming to $\mathbf r_A$, we obtain
$$
\begin{multline}
I=2\sqrt{\frac{m_Am_B}\pi}\iint\mathrm d\mathbf r_A\mathrm d\mathbf r_B\int_0^\infty\mathrm du \\
\times\exp\big[-(\mathbf r_{AB}-\mathbf R_{AB})^\mathrm T\boldsymbol\Omega_m'^{(AB)}(\mathbf r_{AB}-\mathbf R_{AB})-\mathbf r_{AB}^\mathrm T\mathbf U_2\mathbf r_{AB}\big]
\end{multline}
$$
where $\boldsymbol\Omega'^{(AB)}$ absorbed the masses and $\mathbf U_2$ is defined as
$$
\mathbf U_2=u^2\begin{pmatrix}
1&0&0&-1&0&0\\
0&1&0&0&-1&0\\
0&0&1&0&0&-1\\
-1&0&0&1&0&0\\
0&-1&0&0&1&0\\
0&0&-1&0&0&1
\end{pmatrix}
$$
Following only with the integrand, we rearrange terms, and complete the square with respect to $\mathbf r_{AB}$:
$$
\begin{multline}
\exp\big[-\mathbf r_{AB}^\mathrm T(\boldsymbol\Omega_m'^{(AB)}+\mathbf U_2)\mathbf r_{AB}+2\mathbf R_{AB}^\mathrm T\boldsymbol\Omega_m'^{(AB)}\mathbf r_{AB}-\mathbf R_{AB}^\mathrm T\boldsymbol\Omega_m'^{(AB)}\mathbf R_{AB}\big] \\
=\exp\big[-(\mathbf r_{AB}-\mathbf h_{AB})^\mathrm T(\boldsymbol\Omega_m'^{(AB)}+\mathbf U_2)(\mathbf r_{AB}-\mathbf h_{AB})\big] \\
\times\exp\big[-\mathbf R_{AB}^\mathrm T\big(\boldsymbol\Omega_m'^{(AB)}-\boldsymbol\Omega_m'^{(AB)}(\boldsymbol\Omega_m'^{(AB)}+\mathbf U_2)^{-1}\boldsymbol\Omega_m'^{(AB)}\big)\mathbf R_{AB}\big]
\end{multline}
$$
As in the density calculation, the first exponential can be shifted and rotated into a diagonal form, upon which the spatial integrals can be easily evaluated:
$$
\iint\mathrm d\mathbf r_A\mathrm d\mathbf r_B\exp\big[-\mathbf r_{AB}^\mathrm T(\boldsymbol\Omega_m'^{(AB)}+\mathbf U_2)\mathbf r_{AB}\big]=\frac{\pi^3}{\sqrt{\prod_{i=1}^6\lambda_{AB,i}(u)}}
$$
where $\lambda_{AB,i}(u)$ are the eigenvalues of $(\boldsymbol\Omega_m'^{(AB)}+\mathbf U_2)$. The remaining 1-dimensional integral over $u$ from 0 to $\infty$ has a finite integrand that decays exponentially to zero, and so can be readily evaluated numerically.

Putting everything together, we get
$$
\begin{multline}
\langle\Psi|\frac12\sum_{AB}\frac{q_Aq_B}{\lvert\mathbf r_A-\mathbf r_B\rvert}|\Psi\rangle
=\frac12\sum_{AB}q_Aq_B\sqrt{\frac{\prod_{i=1}^{3N}\tilde\omega_i}{\prod_{i=1}^{3(N-2)}\bar\omega_{AB,i}}} \\
\times\int_0^\infty\mathrm du\frac{\exp\big[-\mathbf R_{AB}^\mathrm T\big(\boldsymbol\Omega_m'^{(AB)}-\boldsymbol\Omega_m'^{(AB)}(\boldsymbol\Omega_m'^{(AB)}+\mathbf U_2)^{-1}\boldsymbol\Omega_m'^{(AB)}\big)\mathbf R_{AB}\big]}{\sqrt{\prod_{i=1}^6\lambda_{AB,i}(u)}}
\label{eq:coulombrr}
\end{multline}
$$

Swapping $\mathbf r_B$ for $\mathbf R_B$, we follow a similar path, this time with a 3-dimensional integral instead of the 6-dimensional integral. The result is:
$$
\langle\Psi|\sum_{AB}-\frac{q_Aq_B}{\lvert\mathbf r_A-\mathbf R_B\rvert}|\Psi\rangle
=\sum_{AB}-q_Aq_B\sqrt{\frac{\prod_{i=1}^{3N}\tilde\omega_i}{\prod_{i=1}^{3(N-1)}\bar\omega_{A,i}}}\int_0^\infty\mathrm du\frac{\exp\big[-\mathbf R_{AB}^\mathrm T\bar{\boldsymbol\Omega}_A\mathbf R_{AB}\big]}{\sqrt{\prod_{i=1}^3\lambda_{A,i}(u)}}
\label{eq:coulombrR}
$$
where
$$
\bar{\boldsymbol\Omega}_A=\begin{pmatrix}
\boldsymbol\Omega'^{(A)} & \mathbf0 \\
\mathbf 0 & u^2\mathbf I
\end{pmatrix}+\begin{pmatrix}
\boldsymbol\Omega'^{(A)} \\
u^2\mathbf I
\end{pmatrix}\big(\boldsymbol\Omega'^{(A)}+u^2\mathbf I)^{-1}\begin{pmatrix}
\boldsymbol\Omega'^{(A)} & u^2\mathbf I
\end{pmatrix}
$$
The nucleus–nucleus term reduces trivially:
$$
\langle\Psi|\frac12\sum_{AB}\frac{q_Aq_B}{\lvert\mathbf R_A-\mathbf R_B\rvert}|\Psi\rangle=\frac12\sum_{AB}\frac{q_Aq_B}{\lvert\mathbf R_A-\mathbf R_B\rvert}
\label{eq:coulombRR}
$$
For calculation of $E^{(1)}$, we are missing the last piece: $\langle\Psi|V_\mathbf{pp}|\Psi\rangle$. First, we transform the dipole potential to the coupled basis and gather the prefactors:


$$
\tilde T_{ij}=\sum_{kl}C_{ki}C_{lj}\omega_k\omega_l\sqrt{\alpha_{0,k}\alpha_{0,l}}T_{kl}
$$
Then,
$$
\begin{multline}
\langle\Psi|V_\mathbf{pp}|\Psi\rangle=\langle\Psi|\frac12\sum_{ij}\tilde\xi_i\tilde\xi_j\tilde T_{ij}|\Psi\rangle \\
=\frac12\sum_{i\neq j}\tilde T_{ij}\left(\frac{\tilde\omega_i\tilde\omega_j}{\pi^2}\right)^\frac14\int\mathrm d\tilde\xi_i\tilde\xi_i\exp\left(-\frac12\tilde\omega_i\tilde\xi_i^2\right)\int\mathrm d\tilde\xi_j\tilde\xi_j\exp\left(-\frac12\tilde\omega_j\tilde\xi_j^2\right) \\
+\frac12\sum_{i}\tilde T_{ii}\sqrt{\frac{\tilde\omega_i}{\pi}}\int\mathrm d\tilde\xi_i\tilde\xi_i^2\exp\left(-\tilde\omega_i\tilde\xi_i^2\right)=\sum_i\frac{\tilde T_{ii}}{4\tilde\omega_i}
\label{eq:dipoleterm}
\end{multline}
$$
where the $i\neq j$ terms disappear because the integrands are odd functions.

Putting $\eqref{eq:coulombrr}$, $\eqref{eq:coulombrR}$, $\eqref{eq:coulombRR}$, and $\eqref{eq:dipoleterm}$ together, we have now all terms necessary to calculate the first-order correction to the dipole approximation in $\eqref{eq:MBDhamil}$.