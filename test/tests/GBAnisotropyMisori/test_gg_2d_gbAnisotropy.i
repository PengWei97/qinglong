# TODO - materials中的模型使用考虑copperThinFilm2023中的初始模型

my_filename = 'cs_gbAnisotropy'
my_filename2 = 'cs_gbAnisotropy'

my_interval = 5
my_num_adaptivity = 3

my_length_scale = 1e-8
my_time_scale = 1.0
my_GBmob0 = 5.0e-13 # m^4/(J·s)
my_wGB = 5 # nm

my_grain_num = 20
my_rand_seed = 40
my_nx = 20
my_max = 200

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = ${my_nx} # Number of elements in the x-direction
  ny = ${my_nx} # Number of elements in the y-direction

  xmin = 0    # minimum x-coordinate of the mesh
  xmax = ${my_max} # 1000 maximum x-coordinate of the mesh 2000-400 400 1600
  ymin = 0    # minimum y-coordinate of the mesh
  ymax = ${my_max} # 1000 maximum y-coordinate of the mesh

  elem_type = QUAD4  # Type of elements used in the mesh
  uniform_refine = 0 # Initial uniform refinement of the mesh

  parallel_type = distributed # Periodic BCs distributed replicated
[]


[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 12 # Number of order parameters used
  var_name_base = gr # Base name of grains
  grain_num = ${my_grain_num} #Number of grains

  length_scale = ${my_length_scale}
  time_scale = ${my_time_scale}
[]

[Variables]
  # Variable block, where all variables in the simulation are declared
  [./PolycrystalVariables]
  [../]
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = grn_20_gbAnisotropy.tex
  [../]
  [./voronoi]
    type = PolycrystalVoronoi
    rand_seed = ${my_rand_seed}
    coloring_algorithm = jp
    int_width = ${my_wGB}
  [../]
  [grain_tracker]
    type = GrainTracker
    connecting_threshold = 0.05
    compute_var_to_feature_map = true
    flood_entity_type = ELEMENTAL
    execute_on = 'initial timestep_begin'
  []
  [./term]
    type = Terminator
    expression = 'grain_tracker < 5'
  [../]
[]

[ICs]
  [./PolycrystalICs]
    [./PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle1]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle2]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle3]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  # Kernel block, where the kernels defining the residual equations are set up.
  [./PolycrystalKernel]
    # Custom action creating all necessary kernels for grain growth.  All input parameters are up in GlobalParams
  [../]
[]

[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
  [../]
  [./var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    execute_on = timestep_end
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
  [../]
  [./euler_angle1]
    type = OutputEulerAngles
    variable = euler_angle1
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    execute_on = 'initial timestep_end'
  [../]
  [./euler_angle2]
    type = OutputEulerAngles
    variable = euler_angle2
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'Phi'
    execute_on = 'initial timestep_end'
  [../]
  [./euler_angle3]
    type = OutputEulerAngles
    variable = euler_angle3
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi2'
    execute_on = 'initial timestep_end'
  [../]
[]

[BCs]
  # Boundary Condition block
  [Periodic]
    [All]
      auto_direction = 'x y' # Makes problem periodic in the x and y directions
    []
  []
[]

[Materials]
  [./CuGrGranisotropic]
    type = GBAnisotropyMisori

    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    gb_energy_anisotropy = true
    gb_mobility_anisotropy = true

    T = 673.15 # K
    wGB = ${my_wGB}
  
    GBsigma_HAGB = 0.708
    GBmob_HAGB = ${my_GBmob0}
    Q_HAGB = 0.23

    # rate1_HABvsLAB_mob = 0.6 # rate_HABvsLAB + 1
    # rate2_HABvsLAB_mob = 0.4
    # rate1_HABvsLAB_sigma = 0.8 # rate_HABvsLAB + 1
    # rate2_HABvsLAB_sigma = 0.2

    output_properties = 'kappa_op L mu misori_angle twinning_type'
    outputs = my_exodus
  [../]
[]

[VectorPostprocessors]
  [./grain_volumes] 
    type = FeatureDataVectorPostprocessor 
    flood_counter = grain_tracker # The FeatureFloodCount UserObject to get values from.
    execute_on = 'initial timestep_end'
    output_centroids = true
  [../]
[]

[Postprocessors]
  # Scalar postprocessors
  [./dt]
    # Outputs the current time step
    type = TimestepSize
  [../]
  [./dofs]
    type = NumDOFs
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
  [./avg_grain_volumes]
    type = AverageGrainVolume
    feature_counter = grain_tracker
    execute_on = 'initial timestep_end'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'

  l_max_its = 20 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 12 # Max number of nonlinear iterations
  nl_abs_tol = 1e-11 # Relative tolerance for nonlinear solves
  nl_rel_tol = 1e-10 # Absolute tolerance for nonlinear solves
  dtmin = 1e-5
  start_time = 0.0
  num_steps = 10
  # end_time = 3500

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.05
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
    initial_adaptivity = ${my_num_adaptivity} # Number of times mesh is adapted to initial condition
    refine_fraction = 0.7 # Fraction of high error that will be refined
    coarsen_fraction = 0.1 # Fraction of low error that will coarsened
    max_h_level = ${my_num_adaptivity} # Max number of refinements used, starting from initial mesh (before uniform refinement)
  [../]
[]

[Outputs]
  print_linear_residuals = false
  [./my_checkpoint]
    type = Checkpoint
    file_base = ./${my_filename}/out_${my_filename}
    interval = ${my_interval}
    additional_execute_on = 'FINAL'
  [../]  
  [my_exodus]
    type = Nemesis
    file_base = ./ex_${my_filename2}/out_${my_filename} 
    interval = ${my_interval}
    additional_execute_on = 'FINAL'
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    type = CSV
  [../]
[]