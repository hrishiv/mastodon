# Test for fuelcompression spring material model
# A  displacement time history is applied at one node of a spring element
# The other node is fixed in all directions.
# The force deformation behavior of the spring with complete deteriotation is obtained for monotonic loading
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
  [./dispx]
    type = FunctionDirichletBC
    variable = disp_x
    function = disp
    boundary = right
  [../]
[]

[Functions]
[./disp] #linear
  type = PiecewiseLinear
  x = '0 3.5 7 7.5 8 11.5' # time
  y = '0 -3.5 0 0.5 0 -3.5'  # force
[../]
[]

[Materials]
  [./fuelspring_spring_test]
    type = FuelCompressionSpring
    block = 0
    y_orientation = '0.0 1.0 0.0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    k_trans =  2.0
    yield_force = 4
     th_disp = 3
     Mahin_mod = true
     complete_det = true
     k_sg = 10
     alpha_s = 0.25
     component = 0
     tension_factor = 0
     recovery_factor = 0
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
  end_time = 11.5
  dt = 0.005
  dtmin = 0.001
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

[]

[Outputs]
  exodus = true
  csv = true
  file_base = test6
[]
