# PINSFVMomentumPressureFlux


This object adds the $\nabla (\epsilon p)$ pressure gradient term of the
incompressible porous media Navier Stokes momentum equation as a surface term
using the divergence theorem. Note that the porosity gradient that originates
from the integration by part is added by
`PINSFVMomentumPressurePorosityGradient`.

The pressure is a Lagrange Multiplier that ensures the incompressibility constraint.

!syntax parameters /FVKernels/PINSFVMomentumPressureFlux

!syntax inputs /FVKernels/PINSFVMomentumPressureFlux

!syntax children /FVKernels/PINSFVMomentumPressureFlux
