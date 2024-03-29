[Mesh]
    [First]
        type = GeneratedMeshGenerator
        elem_type = EDGE3
        dim =  1
        xmin = 0
        xmax = 10
        nx = 10
    []
[]
  
[Variables]
    [ur]
        order = SECOND
        family = LAGRANGE
    []
    [uim]
        order = SECOND
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
        coefficient = 1.0
        wave_number = 1.0
        boundary = 'right'
    []
    [Blochuim]
        type = BlochDirichletBCImag
        variable = uim
        ur = 'ur'
        uim = 'uim'
        lattice_length = 10.0
        coefficient  =  1.0
        wave_number = 1.0
        boundary ='right'
    []
[]

[Materials]
    [nm]
      type = GenericConstantMaterial
      prop_names = 'diffusivity'
      prop_values = 0
    []
[]

[Problem]
    type = EigenProblem
    error_on_jacobian_nonzero_reallocation = false
[]

[Executioner]
    type = Eigenvalue
    solve_type = krylovschur
    n_eigen_pairs = 4
    [./Quadrature]
        type = GAUSS_LOBATTO
        order = SECOND
    [../]
[]
  
[Outputs]
    exodus = true
[]
