[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  # This test currently diffs when run in parallel with DistributedMesh enabled,
  # most likely due to the fact that it uses some geometric search stuff.
  # For more information, see #2121.
  parallel_type = replicated
[]

[Variables]
  [u]
  []
[]

[AuxVariables]
  [t]
  []
  [u_from_sub]
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = u
  []
  [td]
    type = TimeDerivative
    variable = u
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = u
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = u
    boundary = right
    value = 1
  []
[]

[Executioner]
  type = Transient
  num_steps = 10
  dt = 0.01

  solve_type = 'PJFNK'

  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]

[MultiApps]
  [sub]
    type = TransientMultiApp
    app_type = MooseTestApp
    execute_on = timestep_end
    positions = '0 0 0'
    input_files = sub.i
    reset_apps = 0
    reset_time = 0.05
  []
[]

[Transfers]
  [t_from_sub]
    type = MultiAppNearestNodeTransfer
    from_multi_app = sub
    source_variable = t
    variable = t
  []
  [u_from_sub]
    type = MultiAppNearestNodeTransfer
    from_multi_app = sub
    source_variable = u
    variable = u_from_sub
  []
  [u_to_sub]
    type = MultiAppNearestNodeTransfer
    to_multi_app = sub
    source_variable = u
    variable = u_from_master
  []
[]
