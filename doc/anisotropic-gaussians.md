### Coulomb interaction of two anisotropic Gaussian charge densities

We want to calculate
$$
I_1(\mathbf K_1,\mathbf K_2)=\frac{\sqrt{\det\mathbf K_1\mathbf K_2}}{\pi^3}\iint\mathrm d\mathbf r_1\mathrm d\mathbf r_2\frac{\mathrm e^{-(\mathbf r_1-\mathbf R_1)^\mathrm T\mathbf K_1(\mathbf r_1-\mathbf R_1)}\mathrm e^{-(\mathbf r_2-\mathbf R_2)^\mathrm T\mathbf K_2(\mathbf r_2-\mathbf R_2)}}{\lvert\mathbf r_1-\mathbf r_2\rvert}
$$
which is normalized such that
$$
\lim_{k\rightarrow\infty} I_1(k\mathbf I,k\mathbf I)=\frac1{\lvert\mathbf R_1-\mathbf R_2\rvert}
$$
We generalize to a more general problem by introducing
$$
\mathbf r=\begin{bmatrix}\mathbf r_1\\\mathbf r_2\end{bmatrix},\qquad
\mathbf R=\begin{bmatrix}\mathbf R_1\\\mathbf R_2\end{bmatrix},\qquad
\mathbf K=\begin{bmatrix}\mathbf K_1&\mathbf 0\\\mathbf 0&\mathbf K_2\end{bmatrix}
$$
and
$$
I_2(\mathbf K)=\frac{\sqrt{\det\mathbf K}}\pi\iint\mathrm d\mathbf r_1\mathrm d\mathbf r_2\frac{\mathrm e^{-(\mathbf r-\mathbf R)^\mathrm T\mathbf K(\mathbf r-\mathbf R)}}{\lvert\mathbf r_1-\mathbf r_2\rvert}
$$
Using
$$
\frac1{\lvert\mathbf r_1-\mathbf r_2\rvert}=\frac2{\sqrt\pi}\int_0^\infty\mathrm du\exp(-\lvert\mathbf r_1-\mathbf r_2\rvert^2u^2)
$$
we transform to
$$
I_2(\mathbf K)=\frac{2\sqrt{\det\mathbf K}}{\pi^\frac72}\iint\mathrm d\mathbf r_1\mathrm d\mathbf r_2\int_0^\infty\mathrm du\exp\left[-(\mathbf r-\mathbf R)^\mathrm T\mathbf K(\mathbf r-\mathbf R)-\lvert\mathbf r_1-\mathbf r_2\rvert^2u^2\right]
\label{eq:I2}
$$
Next, we introduce
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
so that
$$
\lvert\mathbf r_1-\mathbf r_2\rvert^2u^2=\mathbf r^\mathrm T\mathbf U_2\mathbf r
$$
so that we can rewrite $\eqref{eq:I2}$ as
$$
I_2(\mathbf K)=\frac{2\sqrt{\det\mathbf K}}{\pi^\frac72}\int\mathrm d\mathbf r\int_0^\infty\mathrm du\exp\left[-(\mathbf r-\mathbf R)^\mathrm T\mathbf K(\mathbf r-\mathbf R)-\mathbf r^\mathrm T\mathbf U_2\mathbf r\right]
$$
Next, we rearrange the terms and [complete the square](https://en.wikipedia.org/wiki/Completing_the_square#Formula):
$$
\begin{multline}
I_2(\mathbf K)=\frac{2\sqrt{\det\mathbf K}}{\pi^\frac72}\int\mathrm d\mathbf r\int_0^\infty\mathrm du\exp\left[-\mathbf r^\mathrm T(\mathbf K+\mathbf U_2)\mathbf r+2\mathbf R^\mathrm T\mathbf K\mathbf r-\mathbf R^\mathrm T\mathbf K\mathbf R\right] \\
=\frac{2\sqrt{\det\mathbf K}}{\pi^\frac72}\int\mathrm d\mathbf r\int_0^\infty\mathrm du\\
\times\exp\left[-(\mathbf r-\mathbf h)^\mathrm T(\mathbf K+\mathbf U_2)(\mathbf r-\mathbf h)-\mathbf R^\mathrm T\mathbf K\mathbf R+\mathbf R^\mathrm T\mathbf K(\mathbf K+\mathbf U_2)^{-1}\mathbf K\mathbf R\right]
\end{multline}
$$
The first term in the exponential is a 6-dimensional Gaussian, the integral of which is $\sqrt{\pi^3/\det(\mathbf K+\mathbf U_2)}$. This leaves us with a 1-dimensional integral over $u$:
$$
I_2(\mathbf K)=\frac{2}{\sqrt\pi}\int_0^\infty\mathrm du\sqrt{\frac{\det\mathbf K}{\det(\mathbf K+\mathbf U_2)}}\exp\left[-\mathbf R^\mathrm T\left(\mathbf K-\mathbf K(\mathbf K+\mathbf U_2)^{-1}\mathbf K\right)\mathbf R\right]
$$
For verification, we consider the special case of $\mathbf K_i=\mathbf I/2\sigma_i^2$. Then,
$$
I_1(\sigma_1,\sigma_2)=\frac2{\sqrt\pi}\int_0^\infty\mathrm du\left(1+2u^2(\sigma_1^2+\sigma_2^2)\right)^{-\frac32}\exp\left[-\frac{u^2\lvert\mathbf R_1-\mathbf R_2\rvert^2}{1+2u^2(\sigma_1^2+\sigma_2^2)}\right]
$$
Substituting $v^2=u^2/(1+2u^2(\sigma_1^2+\sigma_2^2))$, we obtain
$$
\begin{equation}
\begin{aligned}
I_1(\sigma_1,\sigma_2)&=\frac2{\sqrt\pi}\int_0^{1/\sqrt{2(\sigma_1^2+\sigma_2^2)}}\mathrm dv\exp\left(-v^2\lvert\mathbf R_1-\mathbf R_2\rvert^2\right) \\
&=\operatorname{erf}\Bigg[\frac{\lvert\mathbf R_1-\mathbf R_2\rvert}{\sqrt{2(\sigma_1^2+\sigma_2^2)}}\Bigg]\frac1{\lvert\mathbf R_1-\mathbf R_2\rvert}
\end{aligned}
\end{equation}
$$
which is the expected result.

### Dipole potential

Next, we want to take the dipole derivative, $\boldsymbol\nabla_{\mathbf R_1}\otimes\boldsymbol\nabla_{\mathbf R_2}$, to obtain the dipole potential. We will use
$$
I_2=\int_0^\infty\mathrm du\,i_2(u)
$$

$$
\boldsymbol\nabla_{\mathbf R_1}\otimes\boldsymbol\nabla_{\mathbf R_2}i_2(\mathbf K)=\left[-2\mathbf K_{12}+4(\mathbf K_{11}\mathbf R_1+\mathbf K_{12}\mathbf R_2)\otimes(\mathbf K_{12}\mathbf R_1+\mathbf K_{22}\mathbf R_2)\right]i_2(\mathbf K)
$$
where
$$
\begin{bmatrix}\mathbf K_{11}&\mathbf K_{12}\\\mathbf K_{12}&\mathbf K_{22}\end{bmatrix}=\mathbf K-\mathbf K(\mathbf K+\mathbf U_2)^{-1}\mathbf K
$$
