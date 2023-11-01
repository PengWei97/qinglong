my_filename = "t0_gg_elastic"
my_filename2 = "t0_gg_elastic"

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 20
  ny = 20
  xmax = 180
  ymax = 180
  elem_type = QUAD4
[]

[GlobalParams]
  op_num = 5
  var_name_base = gr
  grain_num = 5

  displacements = 'disp_x disp_y'

  length_scale = 1.0e-9
  pressure_scale = 1.0e6
[]

[Variables]
  [./PolycrystalVariables] # action
  [../]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
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
  [./elastic_strain22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./unique_grains]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./var_indices]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./C1111]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./euler_angle]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
    # ACGrGrPoly ACInterface TimeDerivative
  [../]
  [./PolycrystalElasticDrivingForce]
      # ACGrGrElasticDrivingForce
  [../]
  [./TensorMechanics]
    use_displaced_mesh = true
    displacements = 'disp_x disp_y'
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
    execute_on = timestep_end
  [../]
  [./elastic_strain22]
    type = RankTwoAux
    variable = elastic_strain22
    rank_two_tensor = elastic_strain
    index_i = 1
    index_j = 1
    execute_on = timestep_end
  [../]
  [./stress22]
    type = RankTwoAux
    variable = stress22
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
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
  [./C1111]
    type = RankFourAux
    variable = C1111
    rank_four_tensor = elasticity_tensor
    index_l = 0
    index_j = 0
    index_k = 0
    index_i = 0
    execute_on = timestep_end
  [../]
  [./vonmises_stress]
    type = RankTwoScalarAux
    variable = vonmises_stress
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    execute_on = timestep_end
  [../]
  [./euler_angle]
    type = OutputEulerAngles
    variable = euler_angle
    euler_angle_provider = euler_angle_file
    grain_tracker = grain_tracker
    output_euler_angle = 'phi1'
    execute_on = 'initial timestep_end'
  [../]
[]

[Functions]
  [./tdisp]
    type = ParsedFunction
    expression = 'if (t<90,0.1*t,9.0)'
  [../]
[]

[BCs]
  # [./Periodic]
  #   [./All]
  #     auto_direction = 'x'
  #     # variable = 'gr0 gr1 gr2 gr3 gr4 gr5 gr6 gr7 gr8'
  #     variable = 'gr0 gr1 gr2 gr3 gr4'
  #   [../]
  # [../]
  [./top_displacement]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = top
    function = tdisp
  [../]
  [./x_anchor]
    type = DirichletBC
    variable = disp_x
    boundary = left # 'left right'
    value = 0.0
  [../]
  [./y_anchor]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
[]

[UserObjects]
  [./euler_angle_file]
    type = EulerAngleFileReader
    file_name = grn_36_rand_2D.tex
  [../]
  [./voronoi]
    type = PolycrystalVoronoi
    coloring_algorithm = jp

    int_width = 5
    rand_seed = 10
  [../]
  [./grain_tracker]
    type = GrainTrackerElasticity # GrainTrackerElasticity
    threshold = 0.2
    compute_var_to_feature_map = true
    execute_on = 'initial timestep_begin'

    C_ijkl = '1.27e5 0.708e5 0.708e5 1.27e5 0.708e5 1.27e5 0.7355e5 0.7355e5 0.7355e5'
    fill_method = symmetric9

    flood_entity_type = ELEMENTAL
    euler_angle_provider = euler_angle_file
  [../]
[]

[VectorPostprocessors]
  [./grain_volumes]
    type = FeatureVolumeVectorPostprocessor # FeatureMatePropVectorPostprocessor
    flood_counter = grain_tracker
    # mat_prop = elastic_energy
    execute_on = 'INITIAL TIMESTEP_BEGIN' # TIMESTEP_BEGIN ~ Must exist, otherwise an error will be reported
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    block = 0
    T = 500 # K
    wGB = 5 # nm
    GBmob0 = 2.5e-6 # m^4/(Js) from Schoenfelder 1997
    Q = 0.23 # Migration energy in eV
    # GBMobility = 0.0
    GBenergy = 0.708 # GB energy in J/m^2

    # length_scale = 1.0e-9
    # time_scale = 1.0
  [../]
  [./ElasticityTensor]
    type = ComputePolycrystalElasticityTensor # ComputePolycrystalElasticityTensor
    grain_tracker = grain_tracker
  [../]
  [./strain]
    type = ComputeFiniteStrain # ComputeSmallStrain
    block = 0
    displacements = 'disp_x disp_y'
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress # ComputeLinearElasticStress
    block = 0
  [../]
[]

[Postprocessors]
  [./elastic_strain22]
    type = ElementAverageValue
    variable = elastic_strain22
  [../]
  [./stress22]
    type = ElementAverageValue
    variable = stress22
  [../]
  [./ngrains]
    type = FeatureFloodCount
    variable = bnds
    threshold = 0.7
  [../]
  [./dofs]
    type = NumDOFs
  [../]
  [./dt]
    type = TimestepSize
  [../]
  [./run_time]
    type = PerfGraphData
    section_name = "Root"
    data_type = total
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    coupled_groups = 'disp_x,disp_y'
  [../]
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart -pc_hypre_boomeramg_strong_threshold'
  petsc_options_value = 'hypre boomeramg 31 0.7'
  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 25
  nl_rel_tol = 1.0e-7
  start_time = 0.0
  # dt = 0.001
  # dtmax = 0.1
  # end_time = 100
  num_steps = 3
  
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    growth_factor = 1.2
    cutback_factor = 0.8
    optimal_iterations = 8
  [../]
  [./Adaptivity]
    initial_adaptivity = 2
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 3
  [../]
[]

[Outputs]
  [my_exodus]
    file_base = ./ex_${my_filename2}/out_${my_filename} 
    interval = 2
    type = Nemesis
    additional_execute_on = 'FINAL'
  [../]
  [./csv]
    file_base = ./csv_${my_filename}/out_${my_filename}
    # interval = 5
    type = CSV
  [../]
  [pgraph]
    type = PerfGraphOutput
    interval = 5
    level = 2                     # Default is 1
    heaviest_branch = false        # Default is false
    heaviest_sections = 7         # Default is 0
    execute_on = 'TIMESTEP_END'
    additional_execute_on = 'FINAL'  # Default is "final"
  []
  print_linear_residuals = false
[]

# mpiexec -np 35 ~/projects/qinglong/qinglong-opt -i coupled_elastic_energy_0.i > 01.log
# gdb --args ~/projects/qinglong/qinglong-dbg -i coupled_elastic_energy_cp.i
