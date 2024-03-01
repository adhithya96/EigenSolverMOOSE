[Mesh]
    [First]
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
    [diff1]
        type = Diffusion
        variable = u
    []
    [eigen1]
        type = MassEigenKernel
        variable = u
    []
[]

[BCs]

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