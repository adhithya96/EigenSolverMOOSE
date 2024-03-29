[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
  uniform_refine = 3
[]

[Variables]
  [./u]
  [../]
[]

[AuxVariables]
  [./ssm]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./conv]
    type = Convection
    variable = u
    velocity = '1 0 0'
  [../]
[]

[AuxKernels]
  [./ssm]
    type = MaterialRealAux
    variable = ssm
    property = diffusivity
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  [../]
  [./right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  [../]
[]

[Materials]
  [./ssm]
    type = SpatialStatefulMaterial
    block = 0
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  num_steps = 4
  dt = 1
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_abs_tol = 1e-14
[]

[Adaptivity]
  marker = box
  [./Markers]
    [./box]
      type = BoxMarker
      bottom_left = '0.2 0.2 0'
      top_right = '0.4 0.4 0'
      inside = refine
      outside = coarsen
    [../]
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
