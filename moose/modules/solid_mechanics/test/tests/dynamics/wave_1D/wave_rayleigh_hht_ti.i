# Wave propogation in 1D using HHT time integration in the presence of Rayleigh damping
#
# The test is for an 1D bar element of length 4m  fixed on one end
# with a sinusoidal pulse dirichlet boundary condition applied to the other end.
# alpha, beta and gamma are HHT  time integration parameters
# eta and zeta are mass dependent and stiffness dependent Rayleigh damping
# coefficients, respectively.
# The equation of motion in terms of matrices is:
#
# M*accel + (eta*M+zeta*K)*((1+alpha)*vel-alpha*vel_old)
# +(1+alpha)*K*disp-alpha*K*disp_old = 0
#
# Here M is the mass matrix, K is the stiffness matrix
#
# The displacement at the first, second, third and fourth node at t = 0.1 are
# -7.787499960311491942e-02, 1.955566679096475483e-02 and -4.634888180231294501e-03, respectively.

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 4
  nz = 1
  xmin = 0.0
  xmax = 0.1
  ymin = 0.0
  ymax = 4.0
  zmin = 0.0
  zmax = 0.1
[]


[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[AuxVariables]
  [./vel_x]
  [../]
  [./accel_x]
  [../]
  [./vel_y]
  [../]
  [./accel_y]
  [../]
  [./vel_z]
  [../]
  [./accel_z]
  [../]
[]

[Kernels]
  [./DynamicSolidMechanics]
    displacements = 'disp_x disp_y disp_z'
    hht_alpha = -0.3
    stiffness_damping_coefficient = 0.1
  [../]
  [./inertia_x]
    type = InertialForce
    variable = disp_x
    eta=0.1
    alpha = -0.3
  [../]
  [./inertia_y]
    type = InertialForce
    variable = disp_y
    eta=0.1
    alpha = -0.3
  [../]
  [./inertia_z]
    type = InertialForce
    variable = disp_z
    eta = 0.1
    alpha = -0.3
  [../]

[]

[AuxKernels]
  [./accel_x] # These auxkernels are only to check output
    type = TestNewmarkTI
    displacement = disp_x
    variable = accel_x
    first = false
  [../]
  [./accel_y]
    type = TestNewmarkTI
    displacement = disp_y
    variable = accel_y
    first = false
  [../]
  [./accel_z]
    type = TestNewmarkTI
    displacement = disp_z
    variable = accel_z
    first = false
  [../]
  [./vel_x]
    type = TestNewmarkTI
    displacement = disp_x
    variable = vel_x
  [../]
  [./vel_y]
    type = TestNewmarkTI
    displacement = disp_y
    variable = vel_y
  [../]
  [./vel_z]
    type = TestNewmarkTI
    displacement = disp_z
    variable = vel_z
  [../]
[]


[BCs]
  [./top_y]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value=0.0
  [../]
  [./top_x]
   type = DirichletBC
    variable = disp_x
    boundary = top
    value=0.0
  [../]
  [./top_z]
    type = DirichletBC
    variable = disp_z
    boundary = top
    value=0.0
  [../]
  [./right_x]
    type = DirichletBC
    variable = disp_x
    boundary = right
    value=0.0
  [../]
  [./right_z]
    type = DirichletBC
    variable = disp_z
    boundary = right
    value=0.0
  [../]
  [./left_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value=0.0
  [../]
  [./left_z]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value=0.0
  [../]
  [./front_x]
    type = DirichletBC
    variable = disp_x
    boundary = front
    value=0.0
  [../]
  [./front_z]
    type = DirichletBC
    variable = disp_z
    boundary = front
    value=0.0
  [../]
  [./back_x]
    type = DirichletBC
    variable = disp_x
    boundary = back
    value=0.0
  [../]
  [./back_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value=0.0
  [../]
  [./bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value=0.0
  [../]
  [./bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = bottom
    value=0.0
  [../]
  [./bottom_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = bottom
    function = displacement_bc
  [../]
[]

[Materials]
  [./Elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    fill_method = symmetric_isotropic
    C_ijkl = '1 0'
  [../]

  [./strain]
    type = ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
  [../]

  [./stress]
    type = ComputeLinearElasticStress
    block = 0
  [../]

  [./density]
    type = GenericConstantMaterial
    block = 0
    prop_names = 'density'
    prop_values = '1'
  [../]

[]

[Executioner]
  type = Transient
  start_time = 0
  end_time = 6.0
  l_tol = 1e-12
  nl_rel_tol = 1e-12
  dt = 0.1

  [./TimeIntegrator]
    type = NewmarkBeta
    beta = 0.422
    gamma = 0.8
  [../]
[]


[Functions]
  [./displacement_bc]
    type = PiecewiseLinear
    data_file = 'sine_wave.csv'
    format = columns
  [../]

[]

[Postprocessors]
  [./_dt]
    type = TimestepSize
  [../]
  [./disp_1]
    type = NodalVariableValue
    nodeid = 1
    variable = disp_y
  [../]
  [./disp_2]
    type = NodalVariableValue
    nodeid = 3
    variable = disp_y
  [../]
  [./disp_3]
    type = NodalVariableValue
    nodeid = 10
    variable = disp_y
  [../]
  [./disp_4]
    type = NodalVariableValue
    nodeid = 14
    variable = disp_y
  [../]
[]

[Outputs]
  file_base = 'wave_rayleigh_hht_out'
  exodus = true
  perf_graph = true
[]
