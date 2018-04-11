# MBD forces

$$
E_\text{MBD}=\frac12\operatorname{Tr}\big(\sqrt{\mathbf Q})-3\sum_i\frac{\omega_i}2,\qquad\mathbf Q_{ij}=\omega_i^2\delta_{ij}\mathbf I+\omega_i\omega_j\sqrt{\alpha_{0,i}\alpha_{0,j}}\mathbf T_{ij}
$$

$$
\mathbf Q\equiv\mathbf C\boldsymbol\Lambda\mathbf C^\text T,\qquad\boldsymbol\Lambda\equiv\operatorname{diag}(\{\tilde\omega_i^2\}),\qquad\operatorname{Tr}\big(\sqrt{\mathbf Q}\big)=\sum_i\tilde\omega_i
$$

$$
\partial E_\text{MBD}=\frac14\operatorname{Tr}\big(\mathbf C\boldsymbol\Lambda^{-\frac12}\mathbf C^\text T\partial\mathbf Q)-3\sum_i\frac{\partial\omega_i}2
$$

$$
\begin{equation}
\begin{aligned}
\partial\mathbf Q_{ij}=&2\delta_{ij}\omega_i\partial\omega_i\mathbf I
+\omega_i\omega_j\sqrt{\alpha_{0,i}\alpha_{0,j}}\mathbf T_{ij}\left(\frac{\partial\omega_i}{\omega_i}+\frac{\partial\omega_j}{\omega_j}+\frac12\frac{\partial\alpha_{0,i}}{\alpha_{0,i}}+\frac12\frac{\partial\alpha_{0,j}}{\alpha_{0,j}}\right)\\
&+\omega_i\omega_j\sqrt{\alpha_{0,i}\alpha_{0,j}}\partial\mathbf T_{ij}
\end{aligned}
\end{equation}
$$

$$
\begin{gather}
T_{(ij)ab}\equiv\frac{\partial^2}{\partial R_a\partial R_b}\frac1R=\frac{-3R_aR_b+R^2\delta_{ab}}{R^5},\qquad \mathbf R\equiv \mathbf R_{ij}\equiv\mathbf R_j-\mathbf R_i \\
\frac{\partial\mathbf T_{ij}}{\partial\mathbf R_k}=\frac{\partial\mathbf T}{\partial\mathbf R}(\delta_{jk}-\delta_{ik}),\qquad\frac{\partial T_{ab}}{\partial R_c}=-3\left(\frac{R_a\delta_{bc}+R_b\delta_{ca}+R_c\delta_{ab}}{R^5}-\frac{5R_aR_bR_c}{R^7}\right)
\end{gather}
$$

$$
\begin{equation}
\begin{aligned}
f_{(ij)}&=\frac1{1+\exp\big({-}a(\eta-1)\big)},\qquad\eta=\frac{R_{(ij)}}{S_{\text{vdW}(ij)}}\equiv\frac{R_{(ij)}}{\beta R_{\text{vdW}(ij)}} \\
\frac{\mathrm df}{\mathrm dR_c}&=\frac a{2+2\cosh\big(a(\eta-1)\big)}\frac{\mathrm d\eta}{\mathrm dR_c},\qquad\frac{\mathrm d\eta}{\mathrm dR_c}=\frac{R_c}{RS_\text{vdW}}-\frac{R}{S_\text{vdW}^2}\frac{\mathrm dS_\text{vdW}}{\mathrm dR_c}
\end{aligned}
\end{equation}
$$

$$
\begin{equation}
\begin{gathered}
T^\text{GG}_{(ij)ab}\equiv\frac{\partial^2}{\partial R_a\partial R_b}\frac{\operatorname{erf}(\zeta)}R=\big(\operatorname{erf}(\zeta)-\Theta(\zeta)\big)T_{ab}+2\zeta^2\Theta(\zeta)\frac{R_aR_b}{R^5},
\\\Theta(\zeta)=\frac{2\zeta}{\sqrt\pi}\exp(-\zeta^2),\qquad \zeta=\frac{R_{(ij)}}{\sigma_{(ij)}}
\end{gathered}
\end{equation}
$$

$$
\begin{multline}
\frac{\mathrm d T_{ab}^\text{GG}}{\mathrm dR_c}=2\zeta\Theta(\zeta)\left(T_{ab}+(3-2\zeta^2)\frac{R_aR_b}{R^5}\right)\frac{\mathrm d\zeta}{\mathrm dR_c}\\
+\big(\operatorname{erf}(\zeta)-\Theta(\zeta)\big)\frac{\partial T_{ab}}{\partial R_c}-2\zeta^2\Theta(\zeta)\left(\frac13\frac{\partial T_{ab}}{\partial R_c}+\frac{R_c\delta_{ab}}{R^5}\right),\qquad\frac{\mathrm d\zeta}{\mathrm dR_c}=\frac{R_c}{R\sigma}-\frac R{\sigma^2}\frac{\mathrm d\sigma}{\mathrm dR_c}
\end{multline}
$$

$$
\begin{equation}
\begin{gathered}
\sigma_i(u)=\left(\frac13\sqrt{\frac2\pi}\alpha_i(u)\right)^{\frac13},\qquad\partial\sigma_i=\sigma_i\frac{\partial\alpha_i}{3\alpha_i}\\
\sigma_{ij}(u)=\sqrt{\sigma_i(u)^2+\sigma_j(u)^2},\qquad\partial\sigma_{ij}=\frac{\sigma_i\partial\sigma_i+\sigma_j\partial\sigma_j}{\sigma_{ij}}
\end{gathered}\tag{sigma}
\end{equation}
$$

$$
\begin{equation}
\begin{gathered}
\bar\alpha_i=\tfrac13\operatorname{Tr}\big(\textstyle\sum_j\boldsymbol{\bar\alpha}_{ij}\big),\qquad \boldsymbol{\bar\alpha}=(\boldsymbol\alpha^{-1}+\mathbf T_\text{GG})^{-1} \\
\partial\boldsymbol{\bar\alpha}=-\boldsymbol{\bar\alpha}(\partial\boldsymbol\alpha^{-1}+\partial\mathbf T_\text{GG})\boldsymbol{\bar\alpha}
\end{gathered}
\end{equation}
$$

$$
\omega=\frac{4C_6}{3\alpha_{0}^2},\qquad\partial\omega=\omega\left(\frac{\partial C_6}{C_6}-\frac{2\partial\alpha_0}{\alpha_0}\right)
$$

$$
\bar X=X\left(\frac{\bar\alpha_0}{\alpha_0}\right)^q,\qquad\partial\bar X=\bar X\left(\frac{\partial X}{X}+q\frac{\partial\bar\alpha_0}{\bar\alpha_0}-q\frac{\partial\alpha_0}{\alpha_0}\right)
$$

$$
\alpha(\mathrm iu)=\frac{\alpha_0}{1+u^2/\omega^2},\qquad
\partial\alpha(\mathrm iu)=\alpha(\mathrm iu)\left(\frac{\partial\alpha_0}{\alpha_0}+\frac2\omega\frac{\partial\omega}{1+\omega^2/u^2}\right)
$$

$$
\bar C_6=\frac3\pi\int_0^\infty\mathrm du\,\bar\alpha(u)^2,\qquad\partial\bar C_6=\frac6\pi\int_0^\infty\mathrm du\,\bar\alpha(u)\partial\bar\alpha(u)
$$

