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
    [fixed]
        type = DirichletBC
        variable = u
        boundary = left
        value = 0
    []
    [load]
        type = NeumannBC
        variable = u
        boundary = right
        value = 200
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
