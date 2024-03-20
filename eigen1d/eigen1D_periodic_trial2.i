[GlobalParams]
    displacements = 'ur uim'
[]  

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

[Modules/TensorMechanics/Master]
    [./all]
      strain = SMALL
      planar_formulation = PLANE_STRAIN
      add_variables = true
      generate_output = 'stress_xx stress_xy stress_yy stress_zz strain_xx strain_xy strain_yy strain_zz'
    [../]
[]

[Kernels]
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
        wave_number = 1.0
        boundary = 'right'
    []
    [Blochuim]
        type = BlochDirichletBCImag
        variable = uim
        ur = 'ur'
        uim = 'uim'
        lattice_length = 10.0
        wave_number = 1.0
        boundary ='right'
    []
[]

[Materials]
    [./linear_stress]
        type = ComputeLinearElasticStress
    [../]
    [./elasticity_tensor]
        type = ComputeIsotropicElasticityTensor
        poissons_ratio = 0.3
        youngs_modulus = 1e6
    [../]
[]

[Problem]
    type = EigenProblem
    error_on_jacobian_nonzero_reallocation = false
[]

[Executioner]
    type = Eigenvalue
    solve_type = lanczos
    n_eigen_pairs = 4
    [./Quadrature]
        type = GAUSS_LOBATTO
        order = SECOND
    [../]
[]
  
[Outputs]
    exodus = true
[]
