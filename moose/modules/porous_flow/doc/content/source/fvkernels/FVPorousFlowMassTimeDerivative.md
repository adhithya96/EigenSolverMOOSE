# FVPorousFlowMassTimeDerivative

!syntax description /FVKernels/FVPorousFlowMassTimeDerivative

This `FVKernel` implements the strong form of
\begin{equation*}
  \frac{\partial}{\partial t}\left(\phi\sum_{\beta}S_{\beta}\rho_{\beta}\chi_{\beta}^{\kappa}\right)
\end{equation*}
where all parameters are defined in the [nomenclature](/nomenclature.md).

!alert note
Presently, a first-order accurate implicit Euler time derivative is hard-coded.

!syntax parameters /FVKernels/FVPorousFlowMassTimeDerivative

!syntax inputs /FVKernels/FVPorousFlowMassTimeDerivative

!syntax children /FVKernels/FVPorousFlowMassTimeDerivative
