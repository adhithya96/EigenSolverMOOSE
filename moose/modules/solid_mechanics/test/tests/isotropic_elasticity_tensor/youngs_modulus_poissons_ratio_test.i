[Mesh]
  type = GeneratedMesh
  dim = 3
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
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
  [./stress_11]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Physics/SolidMechanics/QuasiStatic]
  [./all]
    strain = SMALL
    add_variables = true
  [../]
[]

[AuxKernels]
  [./stress_11]
    type = RankTwoAux
    variable = stress_11
    rank_two_tensor = stress
    index_j = 1
    index_i = 1
  [../]
[]

[BCs]
  [./bottom]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  [../]
  [./left]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  [../]
  [./back]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  [../]
  [./top]
    type = DirichletBC
    variable = disp_y
    boundary = top
    value = 0.001
  [../]
[]

[Materials]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    poissons_ratio = 0.1
    youngs_modulus = 1e6
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  l_max_its = 20
  nl_max_its = 10
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
[]
