[Mesh]
  [accg]
    type = AdvancedConcentricCircleGenerator
    num_sectors = 9
    ring_radii = '1 2'
    ring_intervals = '2 2'
    ring_block_ids = '10 15 20'
    ring_block_names = 'inner_tri inner outer'
    external_boundary_id = 100
    external_boundary_name = 'ext'
    create_outward_interface_boundaries = false
  []
  [fpg]
    type = FlexiblePatternGenerator
    inputs = 'accg'
    boundary_type = HEXAGON
    boundary_size = ${fparse 16.0*sqrt(3.0)}
    boundary_sectors = 10
    extra_positions = '0.0 6.0 0.0
                       -3.0 0.0 0.0
                       3.0 0.0 0.0
                       -6.0 -6.0 0.0
                       0.0 -6.0 0.0
                       6.0 -6.0 0.0'
    extra_positions_mg_indices = '0 0 0 0 0 0'
    desired_area = 1.0
  []
[]

[Problem]
  solve = false
[]

[Postprocessors]
  [background]
    type = VolumePostprocessor
    block = 0
  []
  [circle1]
    type = VolumePostprocessor
    block = '10 15'
  []
  [circle2]
    type = VolumePostprocessor
    block = '20'
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]

[Outputs]
  [csv]
    type = CSV
    execute_on = FINAL
    file_base = 'custom_pattern'
  []
[]
