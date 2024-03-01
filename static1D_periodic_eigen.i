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
    [ur]
        order = FIRST
        family = LAGRANGE
    []
    [uim]
        order = FIRST
        family = LAGRANGE
    []
[]

[Kernels]
    [diff1]
        type = Diffusion
        variable = ur
    []
    [diff2]
        type = Diffusion
        variable = uim
    []
    [eigen1]
        type = MassEigenKernel
        variable = ur
    []
    [eigen2]
        type = MassEigenKernel
        variable = uim
    []
[]

[BCs]
    [BlochU]
        type = BlochDirichletBC
        variable = ur
        lattice length = 1
        wave number = 3.14 / 2
        translation vector = 1 0 0
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