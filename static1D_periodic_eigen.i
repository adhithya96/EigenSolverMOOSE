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
        type = MatDiffusion
        variable = ur
        diffusivity = diffusivity 
    []
    [diff2]
        type = MatDiffusion
        variable = uim
        diffusivity = diffusivity 
    []
    [eigen1]
        type = CoefReaction
        variable = ur
        coefficient = -1.0
        extra_vector_tags = 'eigen'
    []
    [eigen2]
        type = CoefReaction
        variable = uim
        coefficient = -1.0
        extra_vector_tags = 'eigen'
    []
[]

[BCs]
    [Blochur]
        type = BlochDirichletBCReal
        variable = ur
        ur = 'ur'
        uim = 'uim'
        lattice_length = 10.0
        wave_number = 1.732
        boundary = 'left'
    []
    [Blochuim]
        type = BlochDirichletBCImag
        variable = uim
        ur = 'ur'
        uim = 'uim'
        lattice_length = 10
        wave_number = 1.732
        boundary ='left'
    []
[]

[Materials]
    [nm]
      type = GenericConstantMaterial
      prop_names = 'diffusivity'
      prop_values = 0.333333333333333333
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
