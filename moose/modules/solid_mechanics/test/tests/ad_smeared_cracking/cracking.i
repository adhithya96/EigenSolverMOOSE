#
# Simple pull test for cracking.
# The stress increases for two steps and then drops to zero.

[Mesh]
  file = cracking_test.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Functions]
  [./displ]
    type = PiecewiseLinear
    x = '0 1 2 3  4'
    y = '0 1 0 -1 0'
  [../]
[]

[Physics/SolidMechanics/QuasiStatic]
  [./all]
    strain = FINITE
    add_variables = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_zx'
    use_automatic_differentiation = true
  [../]
[]

[BCs]
  [./pull]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = 4
    function = displ
  [../]
  [./left]
    type = ADDirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]
  [./bottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = 2
    value = 0.0
  [../]
  [./back]
    type = ADDirichletBC
    variable = disp_z
    boundary = 3
    value = 0.0
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 2.8e7
    poissons_ratio = 0
  [../]
  [./elastic_stress]
    type = ADComputeSmearedCrackingStress
    cracking_stress = 1.68e6
    softening_models = abrupt_softening
  [../]
  [./abrupt_softening]
    type = ADAbruptSoftening
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  petsc_options_iname = '-ksp_gmres_restart -pc_type -sub_pc_type'
  petsc_options_value = '101                asm      lu'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 100
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  l_tol = 1e-5
  start_time = 0.0
  end_time = 0.1
  dt = 0.025
[]

[Outputs]
  exodus = true
[]
