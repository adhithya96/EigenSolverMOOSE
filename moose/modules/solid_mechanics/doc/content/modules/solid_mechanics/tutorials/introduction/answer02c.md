> If you created a large strain version of the input, try and convert it to use
> MOOSE's automatic differentiation system. A few places to look at:
>
> - [!param](/Physics/SolidMechanics/QuasiStatic/QuasiStaticSolidMechanicsPhysics/use_automatic_differentiation) in the solid mechanics quasi-static physics
> - [!param](/BCs/Pressure/PressureAction/use_automatic_differentiation) in the Pressure BC action
> - [ADDirichletBC](ADDirichletBC.md)
> - [ADComputeIsotropicElasticityTensor](ComputeIsotropicElasticityTensor.md)
> - [ADComputeFiniteStrainElasticStress](ADComputeFiniteStrainElasticStress.md)

Here is the converted input:

!listing modules/solid_mechanics/tutorials/introduction/mech_step02a.i

## Input file

### SolidMechanics QuasiStatic Physics

Adding `use_automatic_differentiation = true` here causes the action to build
the automatic differentiation (AD) enabled versions of the materials, kernels,
and output objects it sets up.

### `BCs`

We replaced `DirichletBC` with the AD-enabled `ADDirichletBC` object. Note that
in general it is fine to mix AD and non-AD objects, although you might then end
up with a less than perfect Jacobian. Also keep in mind that you cannot use AD
and non-AD versions of *material properties* (such as the stiffness tensor or the
strain) interchangeably. If an object requests an AD property you need to use the
AD-enabled version of the material to provide it.

In the `Pressure` action we also supplied `use_automatic_differentiation = true`
to have the action build the AD-enabled versions of the individual boundary
condition objects that act on the displacement variables in the problem.

### `Materials`

As mentioned above, when using AD-enabled kernels (added through the quasi static
physics syntax), we must supply them with AD-enabled material properties. That's why we
select the AD-enabled objects to compute stiffness tensor and stress. (The
AD-enabled strain calculator is automatically added by the quasi-static physics.)

Note that in [the first exercise](solid_mechanics/tutorials/introduction/answer02a.md)
we switched to a large strain formulation. In case you didn't make that change
use `ADComputeLinearElasticStress` here and make sure you do not have `strain =
FINITE` in the quasi-static physics.

### `Executioner`

We select `Newton` as the solve type. MOOSE actually sets up the appropriate
preconditioning block for us automatically in this case, which is why we
removed it from the input file here.

> What did you observe when you ran the converted example?

You should see a substantial reduction in linear iterations. This is an
indication that the new Jacobian matrix generated by automatic differentiation
is more accurate than the hand coded Jacobians in the non-AD version.

> Rerun the problem again with a Young's modulus of `1e8`.

Again each non-linear iteration converges with just two linear iterations, but
the problem is still exhibiting a lot of non-linear steps and even some cut time
steps due to the nonlinearity of the large deformation formulation.