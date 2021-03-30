[Mesh]
  type = FileMesh
  file = contactatboth.e
[]

[Variables]
  [./disp_x]
    block = '1 2 3 4 5 6 7 8 9'
  [../]
  [./disp_y]
    block = '1 2 3 4 5 6 7 8 9'
  [../]
  [./disp_z]
    block = '1 2 3 4 5 6 7 8 9'
  [../]
  [./rot_x]
    block = '1 2 3 4 5 6 8'
  [../]
  [./rot_y]
      block = '1 2 3 4 5 6'
  [../]
  [./rot_z]
    block = '1 3 4 5 6 8'
  [../]
[]

[AuxVariables]
  [./spring_forcex]
    order = FIRST
    family = LAGRANGE
    block = '3 4 5 6'
  [../]
  [./spring_forcey]
    order = FIRST
    family = LAGRANGE
    block = '3 4 5 6'
  [../]
  [./spring_forcez]
    order = FIRST
    family = LAGRANGE
    block = '3 4 5 6'
  [../]
  [./contact_forcex]
    order = FIRST
    family = LAGRANGE
    block = '7 9'
  [../]
  [./contact_forcey]
    order = FIRST
    family = LAGRANGE
    block = '7 9'
  [../]
  [./contact_forcez]
    order = FIRST
    family = LAGRANGE
    block = '7 9'
  [../]
  []


[Kernels]
  [./spring_disp_x]
    type = StressDivergenceSpring
    block = '3 4 5 6'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 0
    variable = disp_x
    save_in = spring_forcex
  [../]
  [./spring_disp_y]
    type = StressDivergenceSpring
    block = '3 4 5 6'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 1
    variable = disp_y
    save_in = spring_forcey
  [../]
  [./spring_disp_z]
    type = StressDivergenceSpring
    block = '3 4 5 6'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 2
    variable = disp_z
    save_in = spring_forcez
  [../]
  [./spring_rot_x]
    type = StressDivergenceSpring
    block = '3 4 5 6'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 3
    variable = rot_x
  [../]
  [./spring_rot_y]
    type = StressDivergenceSpring
    block = '3 4 5 6'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 4
    variable = rot_y
  [../]
  [./spring_rot_z]
    type = StressDivergenceSpring
    block = '3 4 5 6'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 5
    variable = rot_z
  [../]
  [./shell_x]
    type = ADStressDivergenceShell
    block = '2 8'
    component = 0
    variable = disp_x
    through_thickness_order = SECOND
  [../]
  [./shell_y]
    type = ADStressDivergenceShell
    block = '2 8'
    component = 1
    variable = disp_y
    through_thickness_order = SECOND
  [../]
  [./shell_z]
    type = ADStressDivergenceShell
    block = '2 8'
    component = 2
    variable = disp_z
    through_thickness_order = SECOND
  [../]
  [./shell_rotx]
    type = ADStressDivergenceShell
    block = '2 8'
    component = 3
    variable = rot_x
    through_thickness_order = SECOND
  [../]
  [./shell_roty]
    type = ADStressDivergenceShell
    block = 2
    component = 4
    variable = rot_y
    through_thickness_order = SECOND
  [../]
  [./shell_rotz]
    type = ADStressDivergenceShell
    block = 8
    component = 4
    variable = rot_z
    through_thickness_order = SECOND
  [../]
  [./contact_disp_x]
    type = StressDivergenceGapContactSpring
    block = '7 9'
    displacements = 'disp_x disp_y disp_z'
    component = 0
    variable = disp_x
    save_in = contact_forcex
  [../]
  [./contact_disp_y]
    type = StressDivergenceGapContactSpring
    block = '7 9'
    displacements = 'disp_x disp_y disp_z'
    component = 1
    variable = disp_y
    save_in = contact_forcey
  [../]
  [./contact_disp_z]
    type = StressDivergenceGapContactSpring
    block = '7 9'
    displacements = 'disp_x disp_y disp_z'
    component = 2
    variable = disp_z
    save_in = contact_forcez
  [../]
[]


[Modules/TensorMechanics/LineElementMaster]
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'


  [./block_2] #fuel rod
    block = 1 #'1 4 6'
    area = 65.67
    Iy = 343.17
    Iz = 343.17
    y_orientation = '0.0 1.0 0.0'
  [../]
  strain_type = FINITE

[]

[Materials]
  [./elasticity_fuel]
    type = ComputeElasticityBeam
    youngs_modulus = 143367.38  #N/mm2 or MPa Zircaloy4
    poissons_ratio =   0.332# 0.332
    block = 1#'1 3 4 5 6'
  [../]
  [./stress_beams]
    type = ComputeBeamResultants
    block = 1#'1 3 4 5 6'
  [../]
  [./dimplev]
    type = LinearSpring
    block = 3
    y_orientation = '0.0 0.0 1.0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    kx = 413.0
    ky = 0
    kz = 0
    krx = 0
    kry = 0
    krz = 0
  [../]
  [./dimpleh]
    type = LinearSpring
    block = 4
    y_orientation = '1.0 0.0 0.0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    kx = 413.0
    ky = 0
    kz = 0
    krx = 0
    kry = 0
    krz = 0
  [../]
  [./springv]
    type = LinearSpring
    block = 5
    y_orientation = '0.0 0.0 1.0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    kx = 49
    ky = 0
    kz = 0
    krx = 0
    kry = 0
    krz = 0
  [../]
  [./springh]
    type = LinearSpring
    block = 6
    y_orientation = '1.0 0.0 0.0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    kx = 49
    ky = 0
    kz = 0
    krx = 0
    kry = 0
    krz = 0
  [../]
  [./elasticityshell]
    type = ADComputeIsotropicElasticityTensorShell
    youngs_modulus = 91650
    poissons_ratio = 0.292
    block = '2 8'
    through_thickness_order = SECOND
  [../]
  [./strainshell]
    # type = ADComputeIncrementalShellStrain
    type = ADComputeFiniteShellStrain
    block = '2'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y'
    thickness = 0.6
    through_thickness_order = SECOND
  [../]
  [./strainshell2]
    # type = ADComputeIncrementalShellStrain
    type = ADComputeFiniteShellStrain
    block = '8'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_z'
    thickness = 0.6
    through_thickness_order = SECOND
  [../]
  [./stressshell]
    type = ADComputeShellStress
    block = '2 8'
    through_thickness_order = SECOND
  [../]
  [./gap_spring_test]
    type = ComputeGapContactSpringElasticity
    block = '7 9'
    y_orientation = '0.0 0.0 1.0'
    displacements = 'disp_x disp_y disp_z'
    Kn = 1e5
    Kt = 0
    mu = 0
  [../]
[]

[BCs]
  [./dis_x]
    type = DirichletBC
    boundary = ' 1 12 13'
    variable = disp_x
    value = 0.0
  [../]
  [./dis_y]
    type = DirichletBC
    boundary = '1 12 13'
    variable = disp_y
    value = 0.0
  [../]
  [./dis_z]
    type = DirichletBC
    boundary = '1 12 13'
    variable = disp_z
    value = 0.0
  [../]
  [./rot_x]
    type = DirichletBC
    boundary = '1 101'
    variable = rot_x
    value = 0.0
  [../]
  [./rot_y]
    type = DirichletBC
    boundary = ' 1 101'
    variable = rot_y
    value = 0.0
  [../]
  [./rot_z]
    type = DirichletBC
    boundary = '1 101'
    variable = rot_z
    value = 0.0
  [../]
[]

[NodalKernels]
  [./force_y2]
    type = ConstantRate
    variable = disp_y
    boundary = '5'
    rate = -25
  [../]
  [./force_x]
    type = ConstantRate
    variable = disp_y
    boundary = '100'
    rate = -125
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
    # petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -snes_atol -snes_rtol -snes_max_it -ksp_atol -ksp_rtol -sub_pc_factor_shift_type'
    # petsc_options_value = 'gmres asm lu 1E-6 1E-6 10 1E-6 1E-6 NONZERO'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
   petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-ksp_gmres_restart -pc_type -pc_hypre_type -pc_hypre_boomeramg_max_iter'
  petsc_options_value = '201                hypre    boomeramg      4'
    line_search = none
  end_time = 1
  dt = 1
  dtmin = 1
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-6
  l_tol = 1e-6
  l_max_its = 50
  nl_max_its = 1000
  timestep_tolerance = 1e-8
[]

[Postprocessors]
  [./disp_tip]
    type = PointValue
    variable = disp_y
    point = '521.97 0 0'
  [../]
  [./disp_tipx]
    type = PointValue
    variable = disp_x
    point = '521.97 0 0'
  [../]
  [./disp_y]
    type = PointValue
    variable = disp_y
  #  point = '238.11 0 0'
    point = '260.985 -6.3 6.3'
  [../]
  [./disp_y2]
    type = PointValue
    variable = disp_y
  #  point = '260.985 0 0'
  point = '260.985 -6.3 0'
  [../]
  [./disp_y3]
    type = PointValue
    variable = disp_y
  #  point = '283.86 0 0'
  point = '260.985 -6.3 -6.3'
  [../]
  [./disp_y_1]
    type = PointValue
    variable = disp_y
    point = '238.11 6.3 0'
  [../]
  [./disp_y2_1]
    type = PointValue
    variable = disp_y
    point = '260.985 -6.3 0'
  [../]
  [./disp_y3_1]
    type = PointValue
    variable = disp_y
    point = '283.86 6.3 0'
  [../]
  [./springforce_y]
    type = PointValue
    point = '238.11 6.3 0'
    variable = spring_forcey
  [../]
  [./springforce_y2]
    type = PointValue
    point = '260.985 0 0'
    variable = spring_forcey
  [../]
  [./springforce_y3]
    type = PointValue
    point = '283.86 6.3 0'
    variable = spring_forcey
  [../]
  [./sg_2]
    type = PointValue
    variable = disp_y
    point = '289.56 -6.3 0'
  [../]
  [./sg_21]
    type = PointValue
    variable = disp_y
    point = '289.56 -6.3 6.3'
  [../]
  [./sg_22]
    type = PointValue
    variable = disp_y
    point = '289.56 -6.3 -6.3'
  [../]
  [./sg_1]
    type = PointValue
    variable = disp_y
    point = '232.41 -6.3 0'
  [../]
  [./sg_12]
    type = PointValue
    variable = disp_y
    point = '232.41 -6.3 6.3'
  [../]
  [./sg_13]
    type = PointValue
    variable = disp_y
    point = '232.41 -6.3 -6.3'
  [../]

  [./cdisp_y]
    type = PointValue
    variable = disp_y
  #  point = '238.11 0 0'
    point = '289.56 -6.3 6.3'
  [../]
  [./contactforcetip]
    type = PointValue
    point = '521.97 -3.2 0'
    variable = contact_forcey
  [../]
  [./contactforce0]
    type = PointValue
    point = '521.97 -3.2 0'
    variable = contact_forcey
  [../]
  [./contactforce1]
    type = PointValue
    point = '289.56 -6.3 6.3'
    variable = contact_forcey
  [../]
  [./contactforce2]
    type = PointValue
    point = '289.56 -6.3 0'
    variable = contact_forcey
  [../]
  [./contactforce3]
    type = PointValue
    point = '283.86 -6.3 6.3'
    variable = contact_forcey
  [../]
  [./contactforce4]
    type = PointValue
    point = '283.86 -6.3 0'
    variable = contact_forcey
  [../]
  [./contactforce5]
    type = PointValue
    point = '260.985 -6.3 6.3'
    variable = contact_forcey
  [../]
  [./contactforce6]
    type = PointValue
    point = '260.985 -6.3 0'
    variable = contact_forcey
  [../]

[]





[Outputs]
  csv = true
  exodus = true
  file_base = papercontact
[]

[Debug]
  show_var_residual_norms = true
  []
