# Base input for testing between-multiapp transfers. It has the following complexities:
# - multiapps may not be run with the same number of ranks
# - both nodal and elemental variables
# - transfers between mixes of nodal and elemental variables
# Tests derived from this input may add or remove complexities through command line arguments

[Problem]
  solve = false
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
[]

# This application use at most 3 processes
[MultiApps/ma1]
  type = TransientMultiApp
  input_files = sub_between_diffusion1.i
  max_procs_per_app = 3
[]

# This application will use as many processes as the main app
[MultiApps/ma2]
  type = TransientMultiApp
  input_files = sub_between_diffusion2.i
[]

[Transfers]
  # Nodal to nodal variables
  [app1_to_2_nodal_nodal]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = ma1
    to_multi_app = ma2
    source_variable = sent_nodal
    variable = received_nodal
  []
  [app2_to_1_nodal_nodal]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = ma2
    to_multi_app = ma1
    source_variable = sent_nodal
    variable = received_nodal
  []

  # Elemental to elemental variables
  [app1_to_2_elem_elem]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = ma1
    to_multi_app = ma2
    source_variable = sent_elem
    variable = received_elem
  []
  [app2_to_1_elem_elem]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = ma2
    to_multi_app = ma1
    source_variable = sent_elem
    variable = received_elem
  []

  # Elemental to nodal variables
  [app1_to_2_elem_nodal]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = ma1
    to_multi_app = ma2
    source_variable = sent_elem
    variable = received_nodal
  []
  [app2_to_1_elem_nodal]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = ma2
    to_multi_app = ma1
    source_variable = sent_elem
    variable = received_nodal
  []

  # Nodal to elemental variables
  [app1_to_2_nodal_elem]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = ma1
    to_multi_app = ma2
    source_variable = sent_nodal
    variable = received_elem
  []
  [app2_to_1_nodal_elem]
    type = MultiAppGeneralFieldNearestLocationTransfer
    from_multi_app = ma2
    to_multi_app = ma1
    source_variable = sent_nodal
    variable = received_elem
  []
[]

[Executioner]
  type = Transient
  num_steps = 1
[]
