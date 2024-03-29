[Mesh]
  [outer_bdy]
    type = PolyLineMeshGenerator
    points = '-1.0 0.0 0.0
              0.0 -1.0 0.0
              1.0 0.0 0.0
              0.0 2.0 0.0'
    loop = true
  []
  [hole_1]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 3
    ny = 3
    xmin = -0.5
    xmax = -0.3
    ymin = -0.1
    ymax = 0.1
  []
  [hole_1_name]
    type = RenameBlockGenerator
    input = hole_1
    old_block = 0
    new_block = hole
  []
  [hole_2]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 3
    ny = 3
    xmin = 0.3
    xmax = 0.5
    ymin = -0.1
    ymax = 0.1
  []
  [hole_2_name_1]
    type = RenameBlockGenerator
    input = hole_2
    old_block = 0
    new_block = 1
  []
  [hole_2_name_2]
    type = RenameBlockGenerator
    input = hole_2_name_1
    old_block = 1
    new_block = hole
  []
  [triang]
    type = XYDelaunayGenerator
    boundary = 'outer_bdy'
    holes = 'hole_1_name
             hole_2'
    stitch_holes = 'true
                    false'
    refine_holes = 'false
                    false'
    add_nodes_per_boundary_segment = 3
    refine_boundary = false
    desired_area = 0.05
    output_subdomain_name = "triangles"
  []
[]
