# Test for linear spring material model
# A sawtooth force history is applied at one node of a spring element in each
# of the 6 degrees of freedom. The displacements and rotations are calculated
# at this node. The other node is fixed in all directions.
# This is a static analysis and therefore, the inertia kernels are ommitted.

[Mesh]
  type = GeneratedMesh
  xmin = 0
  xmax = 1
  nx = 1
  dim = 1
  displacements = 'disp_x disp_y disp_z'
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
  [../]
  [./rot_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_z]
    order = FIRST
    family = LAGRANGE
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
  [./spring_momentx]
    order = FIRST
    family = LAGRANGE
  [../]
  [./spring_momenty]
    order = FIRST
    family = LAGRANGE
  [../]
  [./spring_momentz]
    order = FIRST
    family = LAGRANGE
  [../]
  []

[Kernels]
  [./spring_disp_x]
    type = StressDivergenceSpring
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 0
    variable = disp_x
    save_in = spring_forcex
  [../]
  [./spring_disp_y]
    type = StressDivergenceSpring
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 1
    variable = disp_y
    save_in = spring_forcey
  [../]
  [./spring_disp_z]
    type = StressDivergenceSpring
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 2
    variable = disp_z
    save_in = spring_forcez
  [../]
  [./spring_rot_x]
    type = StressDivergenceSpring
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 3
    variable = rot_x
    save_in = spring_momentx
  [../]
  [./spring_rot_y]
    type = StressDivergenceSpring
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 4
    variable = rot_y
    save_in = spring_momenty
  [../]
  [./spring_rot_z]
    type = StressDivergenceSpring
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 5
    variable = rot_z
    save_in = spring_momentz
  [../]
[]

[BCs]
  [./fixx1]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./fixy1]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0.0
  [../]
  [./fixz1]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0.0
  [../]
  [./fixr1]
    type = DirichletBC
    variable = rot_x
    boundary = left
    value = 0.0
  [../]
  [./fixr2]
    type = DirichletBC
    variable = rot_y
    boundary = left
    value = 0.0
  [../]
  [./fixr3]
    type = DirichletBC
    variable = rot_z
    boundary = left
    value = 0.0
  [../]
[]

[NodalKernels]
  [./force_x]
    type = UserForcingFunctionNodalKernel
    variable = disp_x
    boundary = right
    function = force_x
  [../]
  [./force_y]
    type = UserForcingFunctionNodalKernel
    variable = disp_y
    boundary = right
    function = force_y
  [../]
  [./force_z]
    type = UserForcingFunctionNodalKernel
    variable = disp_z
    boundary = right
    function = force_z
  [../]
  [./moment_x]
    type = UserForcingFunctionNodalKernel
    variable = rot_x
    boundary = right
    function = moment_x
  [../]
  [./moment_y]
    type = UserForcingFunctionNodalKernel
    variable = rot_y
    boundary = right
    function = moment_y
  [../]
  [./moment_z]
    type = UserForcingFunctionNodalKernel
    variable = rot_z
    boundary = right
    function = moment_z
  [../]
[]

[Functions]
  [./force_x]
    type = PiecewiseLinear
    x = '0.0 1.0 1.5 2.0 3.0 4.0 5 6 7' # time
    y = '0.0 10.0 15.0 10.0 0.0 -10 -20 -10 0'  # force
  [../]
  [./force_y]
    type = PiecewiseLinear
    x = '0.0 1.0 1.5 2.0 3.0 4.0 5 6 7' # time
    y = '0.0 10.0 15.0 10.0 0.0 -10 -20 -10 0'  # force
  [../]
  [./force_z]
    type = PiecewiseLinear
    x = '0.0 1.0 1.5 2.0 3.0 4.0 5 6 7' # time
    y = '0.0 10.0 15.0 10.0 0.0 -10 -20 -10 0'  # force
  [../]
  [./moment_x]
    type = PiecewiseLinear
    x = '0.0 1.0 1.5 2.0 3.0 4.0 5 6 7' # time
    y = '0.0 10.0 15.0 10.0 0.0 -10 -20 -10 0'  # moment
  [../]
  [./moment_y]
    type = PiecewiseLinear
    x = '0.0 1.0 1.5 2.0 3.0 4.0 5 6 7' # time
    y = '0.0 10.0 15.0 10.0 0.0 -10 -20 -10 0'  # moment
  [../]
  [./moment_z]
    type = PiecewiseLinear
    x = '0.0 1.0 1.5 2.0 3.0 4.0 5 6 7' # time
    y = '0.0 10.0 15.0 10.0 0.0 -10 -20 -10 0'  # moment
  [../]
[]

[Materials]
  [./linear_spring_test]
    type = IsotropicSpring
    block = 0
    y_orientation = '0.0 1.0 0.0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    kx = 2.0
    ky = 2.5
    kz = 1.0
    krx = 2.0
    kry = 2.5
    krz = 1.0
    yield_force = '10 12.5 15'
    yield_moments = '10 12.5 15'
    hardening_constant_force = '1 1 1'
    hardening_constant_moment = '2 2 2'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
 solve_type = NEWTON
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  start_time = 0.0
  end_time = 7
  dt = 0.005
  dtmin = 0.005
  timestep_tolerance = 1e-6
[]

[Postprocessors]
  [./disp_x]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = disp_x
  [../]
  [./springforce_x]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = spring_forcex
  [../]
  [./disp_y]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = disp_y
  [../]
  [./springforce_y]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = spring_forcey
  [../]
  [./disp_z]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = disp_z
  [../]
  [./springforce_z]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = spring_forcez
  [../]
  [./rot_x]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = rot_x
  [../]
  [./springmoment_x]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = spring_momentx
  [../]
  [./rot_y]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = rot_y
  [../]
  [./springmoment_y]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = spring_momenty
  [../]
  [./rot_z]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = rot_z
  [../]
  [./springmoment_z]
    type = PointValue
    point = '1.0 0.0 0.0'
    variable = spring_momentz
  [../]
[]

[Outputs]
  file_base = isotropicspring
  exodus = true
  csv = true
[]
