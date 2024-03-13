[GlobalParams]
    displacements = 'ux uy'
[]  

[XFEM]
    geometric_cut_userobjects = 'cut_mesh'
    qrule = volfrac
    output_cut_plane = true
[]
  

[Mesh]
    [First]
        type = FileMeshGenerator
        file = square.msh
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
  

[UserObjects]
    [cut_mesh]
      type = InterfaceMeshCut2DUserObject
      mesh_file = circle_surface.e
      interface_velocity_function = vel_func
      #heal_always = true
    []
[]

[Functions]
    [vel_func]
      type = ConstantFunction
      value = 0
    []
[]

[AuxVariables]
    [ls]
      order = FIRST
      family = LAGRANGE
    []
[]  

[Kernels]
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
    [ls]
      type = MeshCutLevelSetAux
      mesh_cut_user_object = cut_mesh
      variable = ls
      execute_on = 'TIMESTEP_BEGIN'
    []
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
    [./linear_stress]
        type = ComputeLinearElasticStress
    [../]
    [./elasticity_tensor]
        type = ComputeIsotropicElasticityTensor
        poissons_ratio = 0.3
        youngs_modulus = 1e10
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