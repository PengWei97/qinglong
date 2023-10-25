# This simulation predicts GB migration of a 2D copper polycrystal with 100 grains represented with 8 order parameters
# Mesh adaptivity and time step adaptivity are used
# An AuxVariable is used to calculate the grain boundary locations
# Postprocessors are used to record time step and the number of grains

[Mesh]
  # Mesh block.  Meshes can be read in or automatically generated
  type = GeneratedMesh
  dim = 2 # Problem dimension
  nx = 10 # Number of elements in the x-direction
  ny = 10 # Number of elements in the y-direction
  xmax = 250 # maximum x-coordinate of the mesh
  ymax = 250 # maximum y-coordinate of the mesh
  elem_type = QUAD4 # Type of elements used in the mesh
  # uniform_refine = 2 # Initial uniform refinement of the mesh
  parallel_type = distributed
[]

[GlobalParams]
  # Parameters used by several kernels that are defined globally to simplify input file
  op_num = 6 # Number of order parameters used
  var_name_base = gr # Base name of grains
[]

[Modules]
  [PhaseField]
    [GrainGrowth]
    []
  []
[]

[UserObjects]
  [voronoi]
    type = PolycrystalVoronoi
    grain_num = 6 # Number of grains
    rand_seed = 5
    int_width = 7
  []
  [grain_tracker]
    type = GrainTracker

    flood_entity_type = ELEMENTAL
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_end'
  []
[]

[ICs]
  [PolycrystalICs]
    [PolycrystalColoringIC]
      polycrystal_ic_uo = voronoi
    []
  []
[]

[AuxVariables]
  # Dependent variables
  [./sigma]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [unique_grains]
    order = CONSTANT
    family = MONOMIAL
  []
  [var_indices]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  # AuxKernel block, defining the equations used to calculate the auxvars
  [bnds_aux]
    # AuxKernel that calculates the GB term
    type = BndsCalcAux
    variable = bnds
    execute_on = 'initial timestep_end'
  []
  [unique_grains]
    type = FeatureFloodCountAux
    variable = unique_grains
    flood_counter = grain_tracker
    field_display = UNIQUE_REGION
    execute_on = 'initial timestep_end'
  []
  [var_indices]
    type = FeatureFloodCountAux
    variable = var_indices
    flood_counter = grain_tracker
    field_display = VARIABLE_COLORING
    execute_on = 'initial timestep_end'
  []
  [./sigma]
    type = MaterialRealAux
    variable = sigma
    property = sigma
    execute_on = timestep_end
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
  [CuGrGr]
    # Material properties
    type = GBEvolution
    T = 450 # Constant temperature of the simulation (for mobility calculation)
    wGB = 14 # Width of the diffuse GB
    GBmob0 = 2.5e-6 #m^4(Js) for copper from Schoenfelder1997
    Q = 0.23 #eV for copper from Schoenfelder1997
    GBenergy = 0.708 #J/m^2 from Schoenfelder1997
  []
[]

[Postprocessors]
  # Scalar postprocessors
  [dt]
    # Outputs the current time step
    type = TimestepSize
  []
[]

[VectorPostprocessors]
  [./grain_volumes]
    type = FeatureMatePropVectorPostprocessor
    flood_counter = grain_tracker
    # mat_prop = sigma
    execute_on = 'initial timestep_begin'
  [../]
[]

[Executioner]
  type = Transient # Type of executioner, here it is transient with an adaptive time step
  scheme = bdf2 # Type of time integration (2nd order backward euler), defaults to 1st order backward euler

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  # Uses newton iteration to solve the problem.
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'

  l_max_its = 50 # Max number of linear iterations
  l_tol = 1e-4 # Relative tolerance for linear solves
  nl_max_its = 10 # Max number of nonlinear iterations

  num_steps = 3
  # end_time = 4000

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 20 # Initial time step.  In this simulation it changes.
    optimal_iterations = 6 # Time step will adapt to maintain this number of nonlinear iterations
  []

  # [Adaptivity]
  #   # Block that turns on mesh adaptivity. Note that mesh will never coarsen beyond initial mesh (before uniform refinement)
  #   initial_adaptivity = 2 # Number of times mesh is adapted to initial condition
  #   refine_fraction = 0.8 # Fraction of high error that will be refined
  #   coarsen_fraction = 0.05 # Fraction of low error that will coarsened
  #   max_h_level = 2 # Max number of refinements used, starting from initial mesh (before uniform refinement)
  # []
[]

[Outputs]
  exodus = true # Exodus file will be outputted
  csv = true

  print_linear_residuals = false
[]
