[Mesh]
    [1D]
        type = GeneratedMeshGenerator
        elem_type = EDGE
        dim =  1
        xmin = 0
        xmax = 10
        nx = 10
    []
[]
  
[Variables]
    [u]
        order = FIRST
        family = LAGRANGE
    []
[]

[Kernels]
    [diff]
        type = Diffusion
        variable = u
    []
[]

[BCs]

[]


[Constraints]
    [bloch]
      type = BlochDirichlet
      variable = u
      a = 10
      k = 3.14
      primary = '0'
      secondary_node_ids = '10'
      penalty = 1e+06
    []
  []
  
[Executioner]
    type = Steady
    solve_type = NEWTON
    [./Quadrature]
        type = GAUSS_LOBATTO
        order = SECOND
    [../]

    petsc_options_iname =  '-pc_type  -pc_hypre_type'
    petsc_options_value = 'hypre  boomeramg'
[]
  
[Outputs]
    exodus = true
[]