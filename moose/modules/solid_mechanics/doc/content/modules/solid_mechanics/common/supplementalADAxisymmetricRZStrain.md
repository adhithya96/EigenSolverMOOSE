Axisymmetric (cylindrical) materials are included in Solid Mechanics for
revolved geometries and assume symmetrical loading. These 'strain calculator'
materials compute the strain within the appropriate coordinate system and rely
on specialized AxisymmetricRZ kernels to handle the stress divergence. This
material supplies material properties with all derivatives required to form an
exact Jacobian.

!alert warning title=Symmetry Assumed About the $z$-axis
The axis of symmetry must lie along the $z$-axis in a $\left(r, z, \theta
\right)$ or cylindrical coordinate system. This symmetry orientation is required
for the calculation of the residual [ADStressDivergenceRZTensors](/ADStressDivergenceRZTensors.md)
for the residual equation and the germane discussion.

The `AxisymmetricRZ` material is appropriate for a 2D simulation and assumes
symmetry revolved about the z-axis. A 2D formulation of an appropriate
simulation problem can reduce the simulation run time while preserving key
physics. Axisymmetric simulations are appropriate to problems in which a solid
is generated by [revolving a planar area about an axis](https://en.wikipedia.org/wiki/Axial_symmetry)
in the same plane.

!alert note title=Use `RZ` Coordinate Type
The coordinate type in the `[Problem]` block of the input file must be set to
`coord_type = RZ`.

## Axisymmetric Strain Formulation

The axisymmetric model employs the cylindrical coordinates, $r$, $z$, and
$\theta$, where the planar cross section formed by the $r$ and $z$ axes is
rotated about the axial $z$ axis, along the length of the cylinder, in the
$\theta$ direction. The cylindrical coordinate system strain tensor for
axisymmetric problems has the form

\begin{equation}
\begin{bmatrix}
\epsilon_{rr} & \epsilon_{rz} & 0 \\
\epsilon_{zr} & \epsilon_{zz} & 0 \\
0 & 0 & \epsilon_{\theta \theta}
\end{bmatrix}
\end{equation}

where the value of the strain $\epsilon_{\theta \theta}$ depends on the
displacement and position in the radial direction

\begin{equation}
\epsilon_{\theta \theta} = \frac{u_r}{X_r}.
\end{equation}

Although axisymmetric problems solve for 3D stress and strain fields, the
problem is mathematically 2D. Using an appropriate set of geometry and boundary
conditions, these types of problems have strain and stress fields which are not
functions of the out of plane coordinate variable.  In the cylindrical
coordinate axisymmetric system, the values of stress and strain in the $\theta$
direction do not depend on the $\theta$ coordinate.

!alert note title=Notation Order Change
The axisymmetric system changes the order of the displacement vector from $(u_r,
u_{\theta}, u_z)$, usually seen in textbooks, to $(u_r, u_z, u_{\theta})$. Take
care to follow this convention in your input files and when adding eigenstrains
or extra stresses.