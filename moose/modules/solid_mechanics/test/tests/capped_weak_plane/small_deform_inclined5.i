# Plastic deformation, shear failure, with inclined normal direction = (1, 0, 0)
# With Young = 10, poisson=0.25 (Lame lambda=4, mu=4)
# applying the following
# deformation to the xmax surface of a unit cube:
# disp_x = 5*t/6
# disp_y = 6*t
# disp_z = 8*t
# should yield trial stress:
# stress_xx = 10*t
# stress_xz = 32*t
# stress_xy = 24*t (so q_trial = 40*t)
# Use tan(friction_angle) = 0.5 and tan(dilation_angle) = 1/6, and cohesion=20,
# the system should return to p=0, q=20, ie stress_xx=0, stress_zx=16,
# stress_yx=12 on the first time step (t=1)
[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -0.5
  zmax = 0.5
[]


[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics/SolidMechanics/QuasiStatic]
  [./all]
    add_variables = true
    incremental = true
    generate_output = 'stress_xx stress_xy stress_xz stress_yy stress_yz stress_zz plastic_strain_xx plastic_strain_xy plastic_strain_xz plastic_strain_yy plastic_strain_yz plastic_strain_zz strain_xx strain_xy strain_xz strain_yy strain_yz strain_zz'
  [../]
[]

[BCs]
  [./bottomx]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./bottomy]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0.0
  [../]
  [./bottomz]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0.0
  [../]

  [./topx]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = right
    function = 5*t/6
  [../]
  [./topy]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = right
    function = 6*t
  [../]
  [./topz]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = right
    function = 8*t
  [../]
[]

[AuxVariables]
  [./f_shear]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_tensile]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./f_compressive]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./intnl_shear]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./intnl_tensile]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./iter]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./ls]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./f_shear]
    type = MaterialStdVectorAux
    property = plastic_yield_function
    index = 0
    variable = f_shear
  [../]
  [./f_tensile]
    type = MaterialStdVectorAux
    property = plastic_yield_function
    index = 1
    variable = f_tensile
  [../]
  [./f_compressive]
    type = MaterialStdVectorAux
    property = plastic_yield_function
    index = 2
    variable = f_compressive
  [../]
  [./intnl_shear]
    type = MaterialStdVectorAux
    property = plastic_internal_parameter
    index = 0
    variable = intnl_shear
  [../]
  [./intnl_tensile]
    type = MaterialStdVectorAux
    property = plastic_internal_parameter
    index = 1
    variable = intnl_tensile
  [../]
  [./iter]
    type = MaterialRealAux
    property = plastic_NR_iterations
    variable = iter
  [../]
  [./ls]
    type = MaterialRealAux
    property = plastic_linesearch_needed
    variable = ls
  [../]
[]

[Postprocessors]
  [./stress_xx]
    type = PointValue
    point = '0 0 0'
    variable = stress_xx
  [../]
  [./stress_xy]
    type = PointValue
    point = '0 0 0'
    variable = stress_xy
  [../]
  [./stress_xz]
    type = PointValue
    point = '0 0 0'
    variable = stress_xz
  [../]
  [./stress_yy]
    type = PointValue
    point = '0 0 0'
    variable = stress_yy
  [../]
  [./stress_yz]
    type = PointValue
    point = '0 0 0'
    variable = stress_yz
  [../]
  [./stress_zz]
    type = PointValue
    point = '0 0 0'
    variable = stress_zz
  [../]
  [./strainp_xx]
    type = PointValue
    point = '0 0 0'
    variable = plastic_strain_xx
  [../]
  [./strainp_xy]
    type = PointValue
    point = '0 0 0'
    variable = plastic_strain_xy
  [../]
  [./strainp_xz]
    type = PointValue
    point = '0 0 0'
    variable = plastic_strain_xz
  [../]
  [./strainp_yy]
    type = PointValue
    point = '0 0 0'
    variable = plastic_strain_yy
  [../]
  [./strainp_yz]
    type = PointValue
    point = '0 0 0'
    variable = plastic_strain_yz
  [../]
  [./strainp_zz]
    type = PointValue
    point = '0 0 0'
    variable = plastic_strain_zz
  [../]
  [./straint_xx]
    type = PointValue
    point = '0 0 0'
    variable = strain_xx
  [../]
  [./straint_xy]
    type = PointValue
    point = '0 0 0'
    variable = strain_xy
  [../]
  [./straint_xz]
    type = PointValue
    point = '0 0 0'
    variable = strain_xz
  [../]
  [./straint_yy]
    type = PointValue
    point = '0 0 0'
    variable = strain_yy
  [../]
  [./straint_yz]
    type = PointValue
    point = '0 0 0'
    variable = strain_yz
  [../]
  [./straint_zz]
    type = PointValue
    point = '0 0 0'
    variable = strain_zz
  [../]
  [./f_shear]
    type = PointValue
    point = '0 0 0'
    variable = f_shear
  [../]
  [./f_tensile]
    type = PointValue
    point = '0 0 0'
    variable = f_tensile
  [../]
  [./f_compressive]
    type = PointValue
    point = '0 0 0'
    variable = f_compressive
  [../]
  [./intnl_shear]
    type = PointValue
    point = '0 0 0'
    variable = intnl_shear
  [../]
  [./intnl_tensile]
    type = PointValue
    point = '0 0 0'
    variable = intnl_tensile
  [../]
  [./iter]
    type = PointValue
    point = '0 0 0'
    variable = iter
  [../]
  [./ls]
    type = PointValue
    point = '0 0 0'
    variable = ls
  [../]
[]


[UserObjects]
  [./coh]
    type = SolidMechanicsHardeningConstant
    value = 20
  [../]
  [./tanphi]
    type = SolidMechanicsHardeningConstant
    value = 0.5
  [../]
  [./tanpsi]
    type = SolidMechanicsHardeningConstant
    value = 0.166666666667
  [../]
  [./t_strength]
    type = SolidMechanicsHardeningConstant
    value = 100
  [../]
  [./c_strength]
    type = SolidMechanicsHardeningConstant
    value = 100
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    fill_method = symmetric_isotropic
    C_ijkl = '4 4'
  [../]
  [./admissible]
    type = ComputeMultipleInelasticStress
    inelastic_models = stress
    perform_finite_strain_rotations = false
  [../]
  [./stress]
    type = CappedWeakInclinedPlaneStressUpdate
    normal_vector = '1 0 0'
    cohesion = coh
    tan_friction_angle = tanphi
    tan_dilation_angle = tanpsi
    tensile_strength = t_strength
    compressive_strength = c_strength
    tip_smoother = 0
    smoothing_tol = 0
    yield_function_tol = 1E-5
  [../]
[]

[Executioner]
  end_time = 1
  dt = 1
  type = Transient
[]

[Outputs]
  file_base = small_deform_inclined5
  csv = true
[]

