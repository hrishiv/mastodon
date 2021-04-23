[Mesh]
  type = FileMesh
  #file = singleonlysgnew.e
   file = simplesupported.e
   allow_renumbering = false

[]



[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./rot_x]
  [../]
  [./rot_y]
  [../]
  [./rot_z]
  [../]
[]

[AuxVariables]
  [./vel_x]
  [../]
  [./vel_y]
  [../]
  [./vel_z]
  [../]
  [./accel_x]
  [../]
  [./accel_y]
  [../]
  [./accel_z]
  [../]
  [./rot_vel_x]
  [../]
  [./rot_vel_y]
  [../]
  [./rot_vel_z]
  [../]
  [./rot_accel_x]
  [../]
  [./rot_accel_y]
  [../]
  [./rot_accel_z]
  [../]
  []

  [AuxKernels]
    [./accel_x]
      type = NewmarkAccelAux
      variable = accel_x
      displacement = disp_x
      velocity = vel_x
      beta = 0.25
      execute_on = 'timestep_end'
    [../]
    [./vel_x]
      type = NewmarkVelAux
      variable = vel_x
      acceleration = accel_x
      gamma = 0.5
      execute_on = 'timestep_end'
    [../]
    [./accel_y]
      type = NewmarkAccelAux
      variable = accel_y
      displacement = disp_y
      velocity = vel_y
      beta = 0.25
      execute_on = 'timestep_end'
    [../]
    [./vel_y]
      type = NewmarkVelAux
      variable = vel_y
      acceleration = accel_y
      gamma = 0.5
      execute_on = 'timestep_end'
    [../]
    [./accel_z]
      type = NewmarkAccelAux
      variable = accel_z
      displacement = disp_z
      velocity = vel_z
      beta = 0.25
      execute_on = 'timestep_end'
    [../]
    [./vel_z]
      type = NewmarkVelAux
      variable = vel_z
      acceleration = accel_z
      gamma = 0.5
      execute_on = 'timestep_end'
    [../]
    [./rot_vel_x]
      type = NewmarkVelAux
      variable = rot_vel_x
      acceleration = rot_accel_x
      gamma = 0.5
      execute_on = timestep_end
    [../]
    [./rot_accel_y]
      type = NewmarkAccelAux
      variable = rot_accel_y
      displacement = rot_y
      velocity = rot_vel_y
      beta = 0.25
      execute_on = timestep_end
    [../]
    [./rot_vel_y]
      type = NewmarkVelAux
      variable = rot_vel_y
      acceleration = rot_accel_y
      gamma = 0.5
      execute_on = timestep_end
    [../]
    [./rot_accel_z]
      type = NewmarkAccelAux
      variable = rot_accel_z
      displacement = rot_z
      velocity = rot_vel_z
      beta = 0.25
      execute_on = timestep_end
    [../]
    [./rot_vel_z]
      type = NewmarkVelAux
      variable = rot_vel_z
      acceleration = rot_accel_z
      gamma = 0.5
      execute_on = timestep_end
    [../]
  []



[Kernels]
  [./solid_disp_x]
    type = StressDivergenceBeam
    block = '1'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 0
    variable = disp_x
  [../]
  [./solid_disp_y]
    type = StressDivergenceBeam
    block = '1'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 1
    variable = disp_y
  [../]
  [./solid_disp_z]
    type = StressDivergenceBeam
    block = '1'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 2
    variable = disp_z
  [../]
  [./solid_rot_x]
    type = StressDivergenceBeam
    block = '1'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 3
    variable = rot_x
  [../]
  [./solid_rot_y]
    type = StressDivergenceBeam
    block = '1'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 4
    variable = rot_y
  [../]
  [./solid_rot_z]
    type = StressDivergenceBeam
    block = '1'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 5
    variable = rot_z
  [../]
  [./inertial_force_x]
    type = InertialForceBeam
    block = 1
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    rotational_velocities = 'rot_vel_x rot_vel_y rot_vel_z'
    rotational_accelerations = 'rot_accel_x rot_accel_y rot_accel_z'
    beta = 0.25
    gamma = 0.5
    area = 65.67
    Iy = 343.174
    Iz = 343.174
    component = 0
    variable = disp_x
  [../]
  [./inertial_force_y]
    type = InertialForceBeam
    block = 1
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    rotational_velocities = 'rot_vel_x rot_vel_y rot_vel_z'
    rotational_accelerations = 'rot_accel_x rot_accel_y rot_accel_z'
    beta = 0.25
    gamma = 0.5
    area = 65.67
    Iy = 343.174
    Iz = 343.174
    component = 1
    variable = disp_y
  [../]
  [./inertial_force_z]
    type = InertialForceBeam
    block = 1
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    rotational_velocities = 'rot_vel_x rot_vel_y rot_vel_z'
    rotational_accelerations = 'rot_accel_x rot_accel_y rot_accel_z'
    beta = 0.25
    gamma = 0.5
    area = 65.67
    Iy = 343.174
    Iz = 343.174
    component = 2
    variable = disp_z
  [../]
  [./inertial_force_rot_x]
    type = InertialForceBeam
    block = 1
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    rotational_velocities = 'rot_vel_x rot_vel_y rot_vel_z'
    rotational_accelerations = 'rot_accel_x rot_accel_y rot_accel_z'
    beta = 0.25
    gamma = 0.5
    area = 65.67
    Iy = 343.174
    Iz = 343.174
    component = 3
    variable = rot_x
  [../]
  [./inertial_force_rot_y]
    type = InertialForceBeam
    block = 1
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    rotational_velocities = 'rot_vel_x rot_vel_y rot_vel_z'
    rotational_accelerations = 'rot_accel_x rot_accel_y rot_accel_z'
    beta = 0.25
    gamma = 0.5
    area = 65.67
    Iy = 343.174
    Iz = 343.174
    component = 4
    variable = rot_y
  [../]
  [./inertial_force_rot_z]
    type = InertialForceBeam
    block = 1
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    velocities = 'vel_x vel_y vel_z'
    accelerations = 'accel_x accel_y accel_z'
    rotational_velocities = 'rot_vel_x rot_vel_y rot_vel_z'
    rotational_accelerations = 'rot_accel_x rot_accel_y rot_accel_z'
    beta = 0.25
    gamma = 0.5
    area = 65.67
    Iy = 343.174
    Iz = 343.174
    component = 5
    variable = rot_z
  [../]

[]



[Materials]
  [./elasticity_fuel]
    type = ComputeElasticityBeam
    youngs_modulus = 143367.38  #N/mm2 or MPa Zircaloy4
    poissons_ratio = 0.0
    block = '1'
  [../]
  [./stress_beams]
    type = ComputeBeamResultants
    block = '1'
  [../]
  [./strain]
  type = ComputeIncrementalBeamStrain
  block = '1'
  displacements = 'disp_x disp_y disp_z'
  rotations = 'rot_x rot_y rot_z'
  area = 65.67
  Iy = 343.174
  Iz = 343.174
  y_orientation = '0.0 1.0 0.0'
[../]
  [./density1]
  type = GenericConstantMaterial
  block = 1
  prop_names = 'density'
  prop_values = '8.27e-5'
[../]
[]

[BCs]

  [./dis_x]
    type = DirichletBC
    boundary = '10'
    variable = disp_x
    value = 0.0
  [../]
  [./dis_y]
    type = DirichletBC
    boundary = '10'
    variable = disp_y
    value = 0.0
  [../]
  [./dis_z]
    type = DirichletBC
    boundary = '10'
    variable = disp_z
    value = 0.0
  [../]
  # [./rot_x]
  #   type = DirichletBC
  #   boundary = '10'
  #   variable = rot_x
  #   value = 0.0
  # [../]
  # [./rot_y]
  #   type = DirichletBC
  #   boundary = '10'
  #   variable = rot_y
  #   value = 0.0
  # [../]
  # [./rot_z]
  #   type = DirichletBC
  #   boundary = '10'
  #   variable = rot_z
  #   value = 0.0
  # [../]
  # [./disp_y]
  #   type = PresetAcceleration
  #   variable = disp_y
  #   velocity = vel_y
  #   acceleration = accel_y
  #   beta = 0.25
  #   function = accel_y
  #   boundary = '5'
  # [../]
  # [./dispy]
  #   type = FunctionDirichletBC
  #   variable = disp_y
  #   function = disp
  #   boundary = 11
  # [../]
[]

[Functions]
  [./disp] #linear
    type = PiecewiseLinear
    x = '0 1 2 3 4 10'
    y = '0 20 0 0 0 0'
  [../]

[]


[NodalKernels]
[./force]
  type = UserForcingFunctionNodalKernel
  function = disp
  boundary = 11
  variable = disp_y
[../]
[]
# [ICs]
#   [./initial_accel]
#     type = FunctionIC
#     function = accel_y
#     variable = accel_y
#     boundary = 5
#   [../]
# []



[Preconditioning]
  [./smp]
    type = SMP
    full = true
    # petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type'
    # petsc_options_value = 'gmres asm lu NONZERO'

  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201                hypre    boomeramg      4'
   line_search = none
    # petsc_options = '-ksp_snes_ew'
    # petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
    # petsc_options_value = '201                hypre    boomeramg      4'
# scheme = explicit-euler
  end_time = 10
  dt = 0.005
  dtmin = 0.000125
  nl_abs_tol = 5e-5
    nl_rel_tol = 1e-8
  l_tol = 1e-8
  l_max_its = 100
  nl_max_its = 50
  timestep_tolerance = 1e-8
  [./TimeIntegrator]
    type = NewmarkBeta
    beta = 0.25
    gamma = 0.5
  [../]
    #  automatic_scaling = true
[]

[Postprocessors]
  [./disp_tip]
    type = PointValue
    variable = disp_y
    point = '1354.175 0 0'
  [../]

[]





[Outputs]
  csv = true
  exodus = true
  file_base = 20simplesupport
  print_linear_residuals = false
  perf_graph = true
#  file_base = singleroddynamics
[]

# [Debug]
#   show_var_residual_norms = true
#   []
