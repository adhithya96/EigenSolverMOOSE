[Mesh]
    [First]
        type = FileMeshGenerator
        elem_type = TwoLines.msh
    []
    [left]
        type = SubdomainBoundingBoxGenerator
        input = file_mesh
        block_id = 0
        bottom_left = '0 0 0'
        top_right = '5 0 0'
    []
    [right]
        type = SubdomainBoundingBoxGenerator
        input = file_mesh
        block_id = 1
        bottom_left = '5 0 0'
        top_right = '10 0 0'
    []
[]
  
[Variables]
    [u]
        order = FIRST
        family = LAGRANGE
    []
[]

[Kernels]
    [diff1]
        type = MatDiffusion
        variable = u
        diffusivity = diffusivity
    []
[]

[BCs]
    [fixed]
        type =  DirichletBC
        variable = u
        boundary = 'left'
        value = 0
    []
    [load]
        type = NeumannBC
        variable = u
        boundary = 'right'
        value = 0
    []
[]

[Materials]
    [nm]
      type = GenericConstantMaterial
      block = 0
      prop_names = 'diffusivity'
      prop_values = 100
    []
[]
  

[Constraints]

[]
  
[Executioner]
    type = Eigenvalue
    eigen_problem_type = gen_non_hermitian
    which_eigen_pairs = SMALLEST_MAGNITUDE
    n_eigen_pairs = 5
    n_basis_vectors = 18
    solve_type = jacobi_davidson
    petsc_options = '-eps_view'
[]
  
[Outputs]
    exodus = true
[]