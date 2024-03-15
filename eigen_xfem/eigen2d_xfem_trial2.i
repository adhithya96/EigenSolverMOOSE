[GlobalParams]
    displacements = 'ux uy'
    order = FIRST
    family = LAGRANGE
[]  

[XFEM]
    geometric_cut_userobjects = 'cut_mesh'
    qrule = volfrac
    output_cut_plane = true
[]

[UserObjects]
  [./level_set_cut_uo]
    type = LevelSetCutUserObject
    level_set_var = ls
    heal_always = true
  [../]
[]

  
[Mesh]
    [First]
        type = FileMeshGenerator
        file = square.msh
    []
[]

[AuxVariables]
  [./ls]
    order = FIRST
    family = LAGRANGE
  [../]
[]  

[AuxKernels]
  [./ls_function]
    type = FunctionAux
    variable = ls
    function = ls_func
  [../]
  [cut_id]
    type = CutSubdomainIDAux
    variables = 'ux uy'
    cut = 'level_set_cut_uo'
    cut_subdomains = '1 2'
  []
[]

[Variables]
  [./ux]
  [../]
  [./uy]
  [../]
[]

[Functions]
  [./ls_func]
    type = ParsedFunction
    expression = 'x*x + y*y - 25'
  [../]
[]

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./a_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./a_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./a_strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./b_strain_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./b_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./b_strain_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'ux uy'
  [../]
  [eigen1]
      type = CoefReaction
      variable = ux
      coefficient = -1.0
      extra_vector_tags = 'eigen'
  []
  [eigen2]
      type = CoefReaction
      variable = uy
      coefficient = -1.0
      extra_vector_tags = 'eigen'
  []
[]

[AuxKernels]
  [./stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 0
    variable = stress_xx
  [../]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy
  [../]
  [./stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 0
    index_j = 1
    variable = stress_xy
  [../]
  [./a_strain_xx]
    type = RankTwoAux
    rank_two_tensor = A_total_strain
    index_i = 0
    index_j = 0
    variable = a_strain_xx
  [../]
  [./a_strain_yy]
    type = RankTwoAux
    rank_two_tensor = A_total_strain
    index_i = 1
    index_j = 1
    variable = a_strain_yy
  [../]
  [./a_strain_xy]
    type = RankTwoAux
    rank_two_tensor = A_total_strain
    index_i = 0
    index_j = 1
    variable = a_strain_xy
  [../]
  [./b_strain_xx]
    type = RankTwoAux
    rank_two_tensor = B_total_strain
    index_i = 0
    index_j = 0
    variable = b_strain_xx
  [../]
  [./b_strain_yy]
    type = RankTwoAux
    rank_two_tensor = B_total_strain
    index_i = 1
    index_j = 1
    variable = b_strain_yy
  [../]
  [./b_strain_xy]
    type = RankTwoAux
    rank_two_tensor = B_total_strain
    index_i = 0
    index_j = 1
    variable = b_strain_xy
  [../]
[]

[Constraints]
  [./dispx_constraint]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = ux
    alpha = 1e8
    geometric_cut_userobject = 'level_set_cut_uo'
  [../]
  [./dispy_constraint]
    type = XFEMSingleVariableConstraint
    use_displaced_mesh = false
    variable = uy
    alpha = 1e8
    geometric_cut_userobject = 'level_set_cut_uo'
  [../]
[]
  
[BCs]
    [fixed1]
        type = DirichletBC
        variable = ux
        boundary = 'boundary'
        value = 0
    []
    [fixed2]
        type = DirichletBC
        variable = uy
        boundary = 'boundary'
        value = 0
    []
[]

[Materials]
    [./elasticity_tensor_A]
      type = ComputeIsotropicElasticityTensor
      base_name = A
      youngs_modulus = 1e9
      poissons_ratio = 0.3
    [../]
    [./strain_A]
      type = ComputeSmallStrain
      base_name = A
    [../]
    [./stress_A]
      type = ComputeLinearElasticStress
      base_name = A
    [../]
    [./elasticity_tensor_B]
      type = ComputeIsotropicElasticityTensor
      base_name = B
      youngs_modulus = 1e-6
      poissons_ratio = 1e-6
    [../]
    [./strain_B]
      type = ComputeSmallStrain
      base_name = B
    [../]
    [./stress_B]
      type = ComputeLinearElasticStress
      base_name = B
    [../]
    [./combined_stress]
      type = XFEMCutSwitchingMaterialRankTwoTensor
      base_names = 'A B'
      level_set_var = ls
      prop_name = stress
      geometric_cut_userobject = 'level_set_cut_uo'
      cut_subdomain_ids =  '1 2'
    [../]
    [./combined_elasticity_tensor]
      type = XFEMCutSwitchingMaterialRankFourTensor
      base_names = 'A B'
      level_set_var = ls
      prop_name = elasticity_tensor
      geometric_cut_userobject = 'level_set_cut_uo'
      cut_subdomain_ids =  '1 2'
    [../]
[]
  

[VectorPostprocessors]
    [eigen]
      type = Eigenvalues
      inverse_eigenvalue = true
    []
[]
  
#[Postprocessors]
#    [fluxintegral]
#      type = ElementIntegralVariablePostprocessor
#      variable = u
#      execute_on = linear
#    []
#[]

[Problem]
    type = EigenProblem
[]

#Linear  eigenvalue solver
[Executioner]
    type = Eigenvalue
    solve_type = arnoldi
    n_eigen_pairs = 4
[]

[Outputs]
    file_base =  eigen_circularhole_out
    exodus =  true
    csv = true
[]