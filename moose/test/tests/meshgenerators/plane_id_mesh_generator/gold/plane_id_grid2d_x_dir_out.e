CDF      
      
len_string     !   len_line   Q   four      	time_step          len_name   !   num_dim       	num_nodes      $   num_elem      
num_el_blk        num_node_sets         num_side_sets         num_el_in_blk1        num_nod_per_el1       num_side_ss1      num_side_ss2      num_side_ss3      num_side_ss4      num_nod_ns1       num_nod_ns2       num_nod_ns3       num_nod_ns4       num_elem_var      num_info  1         api_version       @�
=   version       @�
=   floating_point_word_size            	file_size               int64_status             title         &plane_reporting_id_grid2d_x_dir_out.e      maximum_name_length                     
time_whole                            q4   	eb_status                             	   eb_prop1               name      ID              	   	ns_status         	                    	   ns_prop1      	         name      ID              	$   	ss_status         
                    	4   ss_prop1      
         name      ID              	D   coordx                            	T   coordy                            
t   eb_names                       $      �   ns_names      	                 �      �   ss_names      
                 �      <   
coor_names                         D      �   node_num_map                    �         connect1                  	elem_type         QUAD4        �      �   elem_num_map                    d      $   elem_ss1                          �   side_ss1                          �   elem_ss2                          �   side_ss2                          �   elem_ss3                          �   side_ss3                          �   elem_ss4                              side_ss4                             node_ns1                          (   node_ns2                          @   node_ns3                          X   node_ns4                          p   name_elem_var                          $      �   vals_elem_var1eb1                          �      q<   elem_var_tab                             �   info_records                      `�      �                                                                 @      @              @      @      @$      @$      @)      @)      @.      @.      @              @      @$      @)      @.      @              @      @$      @)      @.      @              @      @$      @)      @.      @              @      @$      @)      @.                      @      @              @              @              @              @      @      @      @      @      @      @      @$      @$      @$      @$      @$      @$      @)      @)      @)      @)      @)      @)      @.      @.      @.      @.      @.      @.                                          bottom                           right                            top                              left                             bottom                           right                            left                             top                                                                                                                             	   
                                                                      !   "   #   $                                          	   
      	         
                                          
         
                                                                                                                                                       !            "   !         #   "         $   #                           	   
                                                                                 
                                                                                                   	                     $          !   "   #   $                   plane_id                               ####################        ?�      ?�      @       @               ?�      ?�   # Created by MOOSE #       ?�      ?�      @       @                             ####################                                                             ### Command Line Arguments ###                                                    ../../../../../reactor-opt -i plane_reporting_id_grid2d.i### Version Info ###   Framework Information:                                                           MOOSE Version:           git commit 47c71914fe on 2021-11-06                     LibMesh Version:                                                                 PETSc Version:           3.15.1                                                  SLEPc Version:           3.15.1                                                  Current Time:            Mon Nov  8 14:15:19 2021                                Executable Timestamp:    Mon Nov  8 11:28:27 2021                                                                                                                                                                                                  ### Input File ###                                                                                                                                                []                                                                                 inactive                       = (no_default)                                    initial_from_file_timestep     = LATEST                                          initial_from_file_var          = INVALID                                         allow_negative_qweights        = 1                                               custom_blocks                  = (no_default)                                    custom_orders                  = (no_default)                                    element_order                  = AUTO                                            order                          = AUTO                                            side_order                     = AUTO                                            type                           = GAUSS                                         []                                                                                                                                                                [AuxKernels]                                                                                                                                                        [./set_plane_id]                                                                   inactive                     = (no_default)                                      isObjectAction               = 1                                                 type                         = ElemExtraIDAux                                    block                        = INVALID                                           boundary                     = INVALID                                           check_boundary_restricted    = 1                                                 control_tags                 = AuxKernels                                        enable                       = 1                                                 execute_on                   = 'LINEAR TIMESTEP_END'                             extra_id_name                = plane_id                                          prop_getter_suffix           = (no_default)                                      seed                         = 0                                                 use_displaced_mesh           = 0                                                 variable                     = plane_id                                        [../]                                                                          []                                                                                                                                                                [AuxVariables]                                                                                                                                                      [./plane_id]                                                                       family                       = MONOMIAL                                          inactive                     = (no_default)                                      initial_condition            = INVALID                                           isObjectAction               = 1                                                 order                        = CONSTANT                                          scaling                      = INVALID                                           type                         = MooseVariableBase                                 initial_from_file_timestep   = LATEST                                            initial_from_file_var        = INVALID                                           block                        = INVALID                                           components                   = 1                                                 control_tags                 = AuxVariables                                      eigen                        = 0                                                 enable                       = 1                                                 fv                           = 0                                                 outputs                      = INVALID                                           use_dual                     = 0                                               [../]                                                                          []                                                                                                                                                                [Executioner]                                                                      auto_preconditioning           = 1                                               inactive                       = (no_default)                                    isObjectAction                 = 1                                               type                           = Steady                                          accept_on_max_fixed_point_iteration = 0                                          accept_on_max_picard_iteration = 0                                               auto_advance                   = INVALID                                         automatic_scaling              = INVALID                                         compute_initial_residual_before_preset_bcs = 0                                   compute_scaling_once           = 1                                               contact_line_search_allowed_lambda_cuts = 2                                      contact_line_search_ltol       = INVALID                                         control_tags                   = (no_default)                                    custom_abs_tol                 = 1e-50                                           custom_pp                      = INVALID                                         custom_rel_tol                 = 1e-08                                           direct_pp_value                = 0                                               disable_fixed_point_residual_norm_check = 0                                      disable_picard_residual_norm_check = 0                                           enable                         = 1                                               fixed_point_abs_tol            = 1e-50                                           fixed_point_algorithm          = picard                                          fixed_point_force_norms        = 0                                               fixed_point_max_its            = 1                                               fixed_point_min_its            = 1                                               fixed_point_rel_tol            = 1e-08                                           l_abs_tol                      = 1e-50                                           l_max_its                      = 10000                                           l_tol                          = 1e-05                                           line_search                    = default                                         line_search_package            = petsc                                           max_xfem_update                = 4294967295                                      mffd_type                      = wp                                              n_max_nonlinear_pingpong       = 100                                             nl_abs_div_tol                 = 1e+50                                           nl_abs_step_tol                = 0                                               nl_abs_tol                     = 1e-50                                           nl_div_tol                     = 1e+10                                           nl_forced_its                  = 0                                               nl_max_funcs                   = 10000                                           nl_max_its                     = 50                                              nl_rel_step_tol                = 0                                               nl_rel_tol                     = 1e-08                                           num_grids                      = 1                                               off_diagonals_in_auto_scaling  = 0                                               outputs                        = INVALID                                         petsc_options                  = INVALID                                         petsc_options_iname            = INVALID                                         petsc_options_value            = INVALID                                         picard_abs_tol                 = 1e-50                                           picard_custom_pp               = INVALID                                         picard_force_norms             = 0                                               picard_max_its                 = 1                                               picard_rel_tol                 = 1e-08                                           relaxation_factor              = 1                                               relaxed_variables              = (no_default)                                    resid_vs_jac_scaling_param     = 0                                               restart_file_base              = (no_default)                                    scaling_group_variables        = INVALID                                         skip_exception_check           = 0                                               snesmf_reuse_base              = 1                                               solve_type                     = INVALID                                         splitting                      = INVALID                                         time                           = 0                                               transformed_postprocessors     = (no_default)                                    transformed_variables          = (no_default)                                    update_xfem_at_timestep_begin  = 0                                               verbose                        = 0                                             []                                                                                                                                                                [Mesh]                                                                             displacements                  = INVALID                                         inactive                       = (no_default)                                    use_displaced_mesh             = 1                                               include_local_in_ghosting      = 0                                               output_ghosting                = 0                                               block_id                       = INVALID                                         block_name                     = INVALID                                         boundary_id                    = INVALID                                         boundary_name                  = INVALID                                         construct_side_list_from_node_list = 0                                           ghosted_boundaries             = INVALID                                         ghosted_boundaries_inflation   = INVALID                                         isObjectAction                 = 1                                               second_order                   = 0                                               skip_deletion_repartition_after_refine = 0                                       skip_partitioning              = 0                                               type                           = FileMesh                                        uniform_refine                 = 0                                               allow_renumbering              = 1                                               build_all_side_lowerd_mesh     = 0                                               centroid_partitioner_direction = INVALID                                         construct_node_list_from_side_list = 1                                           control_tags                   = INVALID                                         dim                            = 1                                               enable                         = 1                                               final_generator                = INVALID                                         ghosting_patch_size            = INVALID                                         max_leaf_size                  = 10                                              nemesis                        = 0                                               parallel_type                  = DEFAULT                                         partitioner                    = default                                         patch_size                     = 40                                              patch_update_strategy          = never                                           skip_refine_when_use_split     = 1                                                                                                                                [./grid]                                                                           inactive                     = (no_default)                                      isObjectAction               = 1                                                 type                         = CartesianMeshGenerator                            control_tags                 = Mesh                                              dim                          = 2                                                 dx                           = '5 10'                                            dy                           = '5 10'                                            dz                           = INVALID                                           enable                       = 1                                                 ix                           = '1 4'                                             iy                           = '1 4'                                             iz                           = INVALID                                           show_info                    = 0                                                 subdomain_id                 = INVALID                                         [../]                                                                                                                                                             [./plane_id_gen]                                                                   inactive                     = (no_default)                                      isObjectAction               = 1                                                 type                         = PlaneIDGenerator                                  control_tags                 = Mesh                                              enable                       = 1                                                 id_name                      = plane_id                                          input                        = grid                                              layer_axis                   = x                                                 layer_coordiantes            = '0 5 15'                                          num_sublayers                = '1 2'                                             show_info                    = 0                                                 tolerance                    = 0.0001                                          [../]                                                                          []                                                                                                                                                                [Mesh]                                                                                                                                                              [./grid]                                                                         [../]                                                                                                                                                             [./plane_id_gen]                                                                 [../]                                                                          []                                                                                                                                                                [Mesh]                                                                                                                                                              [./grid]                                                                         [../]                                                                                                                                                             [./plane_id_gen]                                                                 [../]                                                                          []                                                                                                                                                                [Outputs]                                                                          append_date                    = 0                                               append_date_format             = INVALID                                         checkpoint                     = 0                                               color                          = 1                                               console                        = 1                                               controls                       = 0                                               csv                            = 0                                               dofmap                         = 0                                               execute_on                     = TIMESTEP_END                                    exodus                         = 1                                               file_base                      = plane_reporting_id_grid2d_x_dir_out             gmv                            = 0                                               gnuplot                        = 0                                               hide                           = INVALID                                         inactive                       = (no_default)                                    interval                       = 1                                               json                           = 0                                               nemesis                        = 0                                               output_if_base_contains        = INVALID                                         perf_graph                     = 0                                               perf_graph_live                = 1                                               perf_graph_live_mem_limit      = 100                                             perf_graph_live_time_limit     = 5                                               print_linear_converged_reason  = 1                                               print_linear_residuals         = 1                                               print_mesh_changed_info        = 0                                               print_nonlinear_converged_reason = 1                                             print_perf_log                 = 0                                               show                           = INVALID                                         solution_history               = 0                                               sync_times                     = (no_default)                                    tecplot                        = 0                                               vtk                            = 0                                               xda                            = 0                                               xdr                            = 0                                               xml                            = 0                                             []                                                                                                                                                                [Problem]                                                                          inactive                       = (no_default)                                    isObjectAction                 = 1                                               name                           = 'MOOSE Problem'                                 type                           = FEProblem                                       library_name                   = (no_default)                                    library_path                   = (no_default)                                    object_names                   = INVALID                                         register_objects_from          = INVALID                                         block                          = INVALID                                         control_tags                   = (no_default)                                    coord_type                     = XYZ                                             default_ghosting               = 0                                               enable                         = 1                                               error_on_jacobian_nonzero_reallocation = INVALID                                 extra_tag_matrices             = INVALID                                         extra_tag_solutions            = INVALID                                         extra_tag_vectors              = INVALID                                         force_restart                  = 0                                               fv_bcs_integrity_check         = 1                                               ignore_zeros_in_jacobian       = 0                                               kernel_coverage_check          = 1                                               material_coverage_check        = 1                                               material_dependency_check      = 1                                               near_null_space_dimension      = 0                                               null_space_dimension           = 0                                               parallel_barrier_messaging     = 0                                               previous_nl_solution_required  = 0                                               restart_file_base              = INVALID                                         rz_coord_axis                  = Y                                               skip_additional_restart_data   = 0                                               skip_nl_system_check           = 0                                               solve                          = 0                                               transpose_null_space_dimension = 0                                               use_nonlinear                  = 1                                             []                                                                                  ?�              ?�      ?�      @       @               ?�      ?�      @       @               ?�      ?�      @       @               ?�      ?�      @       @               ?�      ?�      @       @       