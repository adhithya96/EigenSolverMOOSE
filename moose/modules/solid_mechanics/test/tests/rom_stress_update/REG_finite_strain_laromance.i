[GlobalParams]
  displacements = 'disp_r disp_z'
[]

[Problem]
  coord_type = RZ
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 1
  xmax = 2
  nx = 50
  ny = 50
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = FINITE
    incremental = true
    add_variables = true
    eigenstrain_names = 'thermal'
    use_automatic_differentiation = false
  []
[]

[AuxVariables]
  [temp]
    initial_condition = 1000.0
  []
[]

[AuxKernels]
  [cooling]
    type = FunctionAux
    variable = temp
    function = '1000-10*t*x'
  []
[]

[BCs]
  [bottom_fix]
    type = DirichletBC
    variable = disp_z
    boundary = bottom
    value = 0.0
  []
  [left_fix]
    type = DirichletBC
    variable = disp_r
    boundary = left
    value = 0.0
  []
[]

[Materials]
  [eigenstrain]
    type = ComputeThermalExpansionEigenstrain
    eigenstrain_name = 'thermal'
    stress_free_temperature = 1000
    thermal_expansion_coeff = 1e-4 #1e-4
    temperature = temp
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 3.30e11
    poissons_ratio = 0.3
  []
  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = rom_stress_prediction
  []
  [rom_stress_prediction]
    type = SS316HLAROMANCEStressUpdateTest
    temperature = temp
    initial_cell_dislocation_density = 6.0e12
    initial_wall_dislocation_density = 4.4e11
    outputs = all
  []
[]

[Postprocessors]
  [lin_its]
    type = NumLinearIterations
  []
  [total_lin_its]
    type = CumulativeValuePostprocessor
    postprocessor = lin_its
  []
  [nl_its]
    type = NumNonlinearIterations
  []
  [total_nl_its]
    type = CumulativeValuePostprocessor
    postprocessor = nl_its
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-8
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  line_search = 'none'
  end_time = 10
  dt = 1

  automatic_scaling = true
[]

[Outputs]
  # print_linear_converged_reason = false
  # print_nonlinear_converged_reason = false
  # print_linear_residuals = false
  perf_graph = true
[]
