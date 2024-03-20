[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
[]

[AuxVariables]
  [./a]
    family = SCALAR
    order = FIFTH
  [../]
[]

[Variables]
  [./dummy]
  [../]
[]

[Kernels]
  [./dummy]
    type = Diffusion
    variable = dummy
  [../]
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[MultiApps]
  [./sub]
    type = TransientMultiApp
    positions = '0 0 0'
    input_files = 'sub_wrong_order.i'
  [../]
[]

[Transfers]
  [./to_sub]
    type = MultiAppScalarToAuxScalarTransfer
    to_multi_app = sub
    source_variable = 'a'
    to_aux_scalar = 'b'
  [../]
[]

[Outputs]
    exodus = true
[]
