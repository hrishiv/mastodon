
[Mesh]
  type = FileMesh
  file = cant_contactscaled.e
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_x]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]
  [./rot_y]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]
  [./rot_z]
    order = FIRST
    family = LAGRANGE
    block = 1
  [../]
[]

[AuxVariables]
  [./spring_forcex]
    order = FIRST
    family = LAGRANGE
  [../]
  [./spring_forcey]
    order = FIRST
    family = LAGRANGE
  [../]
  [./spring_forcez]
    order = FIRST
    family = LAGRANGE
  [../]
  []

[Modules/TensorMechanics/LineElementMaster]

    add_variables = true
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    use_displaced_mesh = false
    block = 1
  #  dynamic_consistent_inertia = true

  [./block1]
    area = 1e4
    Iy = 8.33e6
    Iz = 8.33e6
    y_orientation = '0.0 1.0 0.0'
  [../]

[]

[Kernels]
  [./spring_disp_x]
    type = StressDivergenceGapContactSpring
    block = '2'
    displacements = 'disp_x disp_y disp_z'
    component = 0
    variable = disp_x
    save_in = spring_forcex
  [../]
  [./spring_disp_y]
    type = StressDivergenceGapContactSpring
    block = '2'
    displacements = 'disp_x disp_y disp_z'
    component = 1
    variable = disp_y
    save_in = spring_forcey
  [../]
  [./spring_disp_z]
    type = StressDivergenceGapContactSpring
    block = '2'
    displacements = 'disp_x disp_y disp_z'
    component = 2
    variable = disp_z
    save_in = spring_forcez
  [../]
[]

[BCs]
  [./fixx1]
    type = DirichletBC
    variable = disp_x
    boundary = '1 3'
    value = 0.0
  [../]
  [./fixy1]
    type = DirichletBC
    variable = disp_y
    boundary = '1 3'
    value = 0.0
  [../]
  [./fixz1]
    type = DirichletBC
    variable = disp_z
    boundary = '1 3'
    value = 0.0
  [../]
  [./fixr1]
    type = DirichletBC
    variable = rot_x
    boundary = '1'
    value = 0.0
  [../]
  [./fixr2]
    type = DirichletBC
    variable = rot_y
    boundary = '1'
    value = 0.0
  [../]
  [./fixr3]
    type = DirichletBC
    variable = rot_z
    boundary = '1'
    value = 0.0
  [../]
[]

[NodalKernels]
  [./force_y2]
    type = ConstantRate
    variable = disp_y
    boundary = '2'
    rate = -7.5e4
  [../]
  # [./force_x]
  #   type = ConstantRate
  #   variable = disp_x
  #   boundary = 2
  #   rate = 30
  #
  # [../]
[]

[Materials]
  [./elasticity_1]
    type = ComputeElasticityBeam
    youngs_modulus = 2.10e5 #n/mm2
    poissons_ratio = 0
    block = 1
  [../]
  [./stress]
    type = ComputeBeamResultants
    block = 1
    outputs = exodus
    output_properties = 'forces moments'
  [../]
  [./gap_spring_test]
    type = ComputeGapContactSpringElasticity
    block = 2
    y_orientation = '0.0 0.0 1.0'
    displacements = 'disp_x disp_y disp_z'
    Kn = 2.1e8 #4252.5
    Kt = 0
    mu = 0
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_atol -snes_rtol -snes_max_it -ksp_atol -ksp_rtol -sub_pc_factor_shift_type'
    petsc_options_value = 'gmres asm lu 1E-6 1E-6 100 1E-6 1E-6 NONZERO'
  [../]
[]

[Executioner]
  type = Transient
 solve_type = PJFNK
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  start_time = 0.0
  l_abs_tol = 1e-10
  end_time = 1
#2.05 #3.8
  dt = 1
  dtmin = 1
  timestep_tolerance = 1e-6
[]

[Postprocessors]
  [./disp_y]
    type = PointValue
    point = '1000 0.0 0.0'
    variable = disp_y
  [../]
  [./disp_x]
    type = PointValue
    point = '1000 0.0 0.0'
    variable = disp_x
  [../]
  [./disp_z]
    type = PointValue
    point = '1000 0.0 0.0'
    variable = disp_z
  [../]
  [./forces_x]
    type = PointValue
    point = '1000 0.0 0.0'
    variable = forces_x
  [../]
  [./forces_y]
    type = PointValue
    point = '1000 0.0 0.0'
    variable = forces_y
  [../]
  [./forces_z]
    type = PointValue
    point = '1000 0.0 0.0'
    variable = forces_z
  [../]
  [./springforce_x]
    type = PointValue
    point = '1000 -1.0 0.0'
    variable = spring_forcex
  [../]
  [./springforce_y]
    type = PointValue
    point = '1000 -1.0 0.0'
    variable = spring_forcey
  [../]
  [./springforce_z]
    type = PointValue
    point = '1000 -1.0 0.0'
    variable = spring_forcez
  [../]
[]

[Outputs]
  exodus = true
[]
