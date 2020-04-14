// MASTODON includes
#include "IsotropicSpring.h"

// MOOSE includes
#include "MooseMesh.h"
#include "Assembly.h"
#include "NonlinearSystem.h"
#include "MooseVariable.h"
#include "RankTwoTensor.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("MastodonApp", IsotropicSpring);

template <>
InputParameters
validParams<IsotropicSpring>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription(
      "Compute the deformations, forces and stiffness matrix of a two-noded spring element.");
  params.addRequiredCoupledVar(
      "rotations",
      "The rotation variables appropriate for the simulation geometry and coordinate system.");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacement variables appropriate for the simulation geometry and coordinate system.");
  params.addRequiredParam<RealGradient>("y_orientation",
                                        "Orientation of the y direction along "
                                        "which, Ky is provided. This should be "
                                        "perpendicular to the axis of the spring.");
  params.addRequiredCoupledVar("kx", "Axial stiffness of the spring");
  params.addRequiredCoupledVar("ky", "Shear stiffness in the y direction of the spring.");
  params.addRequiredCoupledVar("kz", "Shear stiffness in the z direction of the spring.");
  params.addRequiredCoupledVar("krx", "Torsional stiffness of the spring.");
  params.addRequiredCoupledVar("kry", "Rotational stiffness in the y direction of the spring.");
  params.addRequiredCoupledVar("krz", "Rotational stiffness in the z direction of the spring.");
  params.set<MooseEnum>("constant_on") = "ELEMENT"; // set _qp to 0
  params.addRequiredParam<Real>("yield_force",
                               "Yield force after which plastic strain starts accumulating");
 params.addParam<Real>("hardening_constant", 0.0, "Hardening slope");
 params.addParam<FunctionName>("hardening_function",
                               "Engineering stress as a function of plastic strain");
 params.addParam<Real>(
     "absolute_tolerance", 1e-10, "Absolute convergence tolerance for Newton iteration");
 params.addParam<Real>(
     "relative_tolerance", 1e-8, "Relative convergence tolerance for Newton iteration");
  return params;
}

IsotropicSpring::IsotropicSpring(const InputParameters & parameters)
  : Material(parameters),
    _nrot(coupledComponents("rotations")),
    _ndisp(coupledComponents("displacements")),
    _rot_num(3),
    _disp_num(3),
    _deformations(declareProperty<RealVectorValue>("deformations")),
    _deformations_old(getMaterialPropertyOld<RealVectorValue>("deformations")),
    _rotations(declareProperty<RealVectorValue>("rotations")),
    _rotations_old(getMaterialPropertyOld<RealVectorValue>("deformations")),
    _kx(coupledValue("kx")),
    _ky(coupledValue("ky")),
    _kz(coupledValue("kz")),
    _krx(coupledValue("krx")),
    _kry(coupledValue("kry")),
    _krz(coupledValue("krz")),
    _spring_forces_global(declareProperty<RealVectorValue>("global_forces")),
      _spring_forces_global_old(getMaterialPropertyOld<RealVectorValue>("global_forces")),
    _spring_moments_global(declareProperty<RealVectorValue>("global_moments")),
      _spring_moments_global_old(getMaterialPropertyOld<RealVectorValue>("global_moments")),
    _kdd(declareProperty<RankTwoTensor>("displacement_stiffness_matrix")),
    _krr(declareProperty<RankTwoTensor>("rotation_stiffness_matrix")),
    _total_global_to_local_rotation(
        declareProperty<RankTwoTensor>("total_global_to_local_rotation")),
    _yield_force(getParam<Real>("yield_force")),
    _hardening_constant(getParam<Real>("hardening_constant")),
    _hardening_function(isParamValid("hardening_function") ? &getFunction("hardening_function")
                                                           : NULL),
    _absolute_tolerance(parameters.get<Real>("absolute_tolerance")),
    _relative_tolerance(parameters.get<Real>("relative_tolerance")),
    // _total_stretch_old(getMaterialPropertyOld<Real>(_base_name + "total_stretch")),
    _plastic_deformation(declareProperty<RealVectorValue>("_plastic_deformation")),
    _plastic_deformation_old(getMaterialPropertyOld<RealVectorValue>("_plastic_deformation")),
    _plastic_rotation(declareProperty<RealVectorValue>("_plastic_rotation")),
    _plastic_rotation_old(getMaterialPropertyOld<RealVectorValue>("_plastic_rotation")),
    // _stress_old(getMaterialPropertyOld<Real>(_base_name + "axial_stress")),
    _hardening_variable(declareProperty<RealVectorValue>("hardening_variable")),
    _hardening_variable_old(getMaterialPropertyOld<RealVectorValue>("hardening_variable")),
    _elastic_deformation(declareProperty<RealVectorValue>("elastic_deformation")),
      _elastic_rotation(declareProperty<RealVectorValue>("elastic_rotation")),
    _max_its(1000)
{
  // Checking for consistency between length of the provided displacements and rotations vector
  if (_ndisp != _nrot)
    mooseError("IsotropicSpring: The number of variables supplied in 'displacements' "
               "and 'rotations' input parameters must be equal.");

  // Fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    MooseVariable * disp_variable = getVar("displacements", i);
    _disp_num[i] = disp_variable->number(); // Displacement variable numbers in MOOSE

    MooseVariable * rot_variable = getVar("rotations", i);
    _rot_num[i] = rot_variable->number(); // Rotation variable numbers in MOOSE
  }

  if (!parameters.isParamSetByUser("hardening_constant") && !isParamValid("hardening_function"))
   mooseError("PlasticTruss: Either hardening_constant or hardening_function must be defined");

 if (parameters.isParamSetByUser("hardening_constant") && isParamValid("hardening_function"))
   mooseError("PlasticTruss: Only the hardening_constant or only the hardening_function can be "
              "defined but not both");
}

void
IsotropicSpring::initQpStatefulProperties()
{
  _deformations[_qp].zero();
_spring_forces_global[_qp].zero();
 _elastic_deformation[_qp] = 0.0;
}

void
IsotropicSpring::computeQpProperties()
{
  // Compute initial orientation and length of the spring in global coordinate system
  // Fetch the two nodes of the link element
  std::vector<const Node *> node;
  for (unsigned int i = 0; i < 2; ++i)
    node.push_back(_current_elem->node_ptr(i));
  RealGradient x_orientation;
  for (unsigned int i = 0; i < _ndisp; ++i)
    x_orientation(i) = (*node[1])(i) - (*node[0])(i);
  x_orientation /= x_orientation.norm(); // Normalizing with length to get orientation

  // Get y orientation of the spring in global coordinate system
  RealGradient y_orientation = getParam<RealGradient>("y_orientation");
  Real dot = x_orientation * y_orientation;

  // Check if x and y orientations are perpendicular
  if (abs(dot) > 1e-4)
    mooseError("Error in IsotropicSpring: y_orientation should be perpendicular to "
               "the axis of the beam.");

  // Calculate z orientation in the global coordinate system as a cross product of the x and y
  // orientations
  RealGradient z_orientation = x_orientation.cross(y_orientation);

  // Calculate the rotation matrix from global to spring local configuration at t = 0
  _original_global_to_local_rotation(0, 0) = x_orientation(0);
  _original_global_to_local_rotation(0, 1) = x_orientation(1);
  _original_global_to_local_rotation(1, 0) = y_orientation(0);
  _original_global_to_local_rotation(1, 1) = y_orientation(1);
  _original_global_to_local_rotation(0, 2) = x_orientation(2);
  _original_global_to_local_rotation(1, 2) = y_orientation(2);
  _original_global_to_local_rotation(2, 0) = z_orientation(0);
  _original_global_to_local_rotation(2, 1) = z_orientation(1);
  _original_global_to_local_rotation(2, 2) = z_orientation(2);

  // _qp = 0
  computeDeformations();

  // computing spring forces
  computeForces();

  // computing spring stiffness matrix
  computeStiffnessMatrix();
}

void
IsotropicSpring::computeTotalRotation()
{
  _qp = 0;

  // Currently this forumation is limited to small deformations in the spring,
  // namely, it is assumed that there are no rigid body rotations in the spring,
  // and that the total rotation matrix from global to local coordinates
  // (calculated below) remains the same as the one at t = 0 throughout the
  // duration of the simulation.
  _total_global_to_local_rotation[_qp] = _original_global_to_local_rotation;
}

void
IsotropicSpring::computeDeformations()
{
  // fetch the two end nodes for _current_elem
  std::vector<const Node *> node;
  for (unsigned int i = 0; i < 2; ++i)
    node.push_back(_current_elem->node_ptr(i));

  // Fetch the solution for the two end nodes at time t
  NonlinearSystemBase & nonlinear_sys = _fe_problem.getNonlinearSystemBase();
  const NumericVector<Number> & sol = *nonlinear_sys.currentSolution();

  // std::cout<<" \n solution " << sol;

  // Calculating global displacements and rotations at the end nodes
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _global_disp0(i) = sol(node[0]->dof_number(nonlinear_sys.number(), _disp_num[i], 0));
    _global_disp1(i) = sol(node[1]->dof_number(nonlinear_sys.number(), _disp_num[i], 0));
    _global_rot0(i) = sol(node[0]->dof_number(nonlinear_sys.number(), _rot_num[i], 0));
    _global_rot1(i) = sol(node[1]->dof_number(nonlinear_sys.number(), _rot_num[i], 0));
  }

  // Convert spring nodal displacements and rotations from global coordinate system to local
  // coordinate system First, compute total rotation
  computeTotalRotation();
  _local_disp0 = _total_global_to_local_rotation[_qp] * _global_disp0;
  _local_disp1 = _total_global_to_local_rotation[_qp] * _global_disp1;
  _local_rot0 = _total_global_to_local_rotation[_qp] * _global_rot0;
  _local_rot1 = _total_global_to_local_rotation[_qp] * _global_rot1;

  // Calculating spring deformations and rotations in the local
  // coordinate system. Deformations and rotations are assumed to be constant
  // through the length of the spring.
  _deformations[_qp] = _local_disp1 - _local_disp0;
    // std::cout<<" \n deformation " << _deformations[_qp](0);
  _rotations[_qp] = _local_rot1 - _local_rot0;
}

void
IsotropicSpring::computeForces()
// spring forces = Kdd * deformations
// spring moments = Krr * rotations
{
  for (unsigned int i = 0; i < 2; ++i)
{
Real deformation_increment = _deformations[_qp](i) - _deformations_old[_qp](i);
    // std::cout<<" \n deformation increment " << deformation_increment ;
Real trial_force = _spring_forces_global_old[_qp](i) + _kx[_qp] * deformation_increment;

_hardening_variable[_qp](i) = _hardening_variable_old[_qp](i);
// std::cout<<" \n hardening variable " << _hardening_variable[_qp] ;
_plastic_deformation[_qp](i) = _plastic_deformation_old[_qp](i);
// std::cout<<" \n plastic deformation " << _plastic_deformation[_qp] ;

Real yield_condition = std::abs(trial_force) - _hardening_variable[_qp](i) - _yield_force;
// std::cout<<" \n yield condition " << yield_condition;

Real iteration = 0;
Real plastic_deformation_increment = 0.0;
Real elastic_deformation_increment = deformation_increment;
// std::cout<<" \n elastic deformation increment " << elastic_deformation_increment;

if (yield_condition > 0.0)
{
  Real residual = std::abs(trial_force) - _hardening_variable[_qp](i) - _yield_force -
                  _kx[_qp] * plastic_deformation_increment;

  Real reference_residual =
      std::abs(trial_force) - _kx[_qp] * plastic_deformation_increment;

  while (std::abs(residual) > _absolute_tolerance ||
         std::abs(residual / reference_residual) > _relative_tolerance)
  {
    _hardening_variable[_qp](i) = computeHardeningValue(plastic_deformation_increment,i);
    Real hardening_slope = computeHardeningDerivative(plastic_deformation_increment,i);

    Real scalar = (std::abs(trial_force) - _hardening_variable[_qp](i) - _yield_force -
                   _kx[_qp] * plastic_deformation_increment) /
                  (_kx[_qp] + hardening_slope);
    plastic_deformation_increment += scalar;

    residual = std::abs(trial_force) - _hardening_variable[_qp](i) - _yield_force -
               _kx[_qp] * plastic_deformation_increment;

    reference_residual = std::abs(trial_force) - _kx[_qp] * plastic_deformation_increment;

    ++iteration;
    if (iteration > _max_its) // not converging
      throw MooseException("IsotropicSpring: Plasticity model did not converge");
  }
  plastic_deformation_increment *= MathUtils::sign(trial_force);
  _plastic_deformation[_qp](i) += plastic_deformation_increment;
  elastic_deformation_increment = deformation_increment - plastic_deformation_increment;
}
_elastic_deformation[_qp](i) = _deformations[_qp](i) - _plastic_deformation[_qp](i); //Confused here
// std::cout<<" \n elastic deformation " << _elastic_deformation[_qp];

_spring_forces_local(i) = _spring_forces_global_old[_qp](i) + _kx[_qp] * elastic_deformation_increment;
// std::cout<<" \n spring forces local " << _spring_forces_local(0);
}

  // forces
  // _spring_forces_local(0) = _kx[_qp] * _deformations[_qp](0);
  // _spring_forces_local(1) = _ky[_qp] * _deformations[_qp](1);
  // _spring_forces_local(2) = _kz[_qp] * _deformations[_qp](2);
  // convert local forces to global
  _spring_forces_global[_qp] =
      _total_global_to_local_rotation[_qp].transpose() * _spring_forces_local;
// std::cout<<" \n spring forces global" << _spring_forces_global[_qp];

  // moments
  for (unsigned int i = 0; i < 2; ++i)
{
Real rotation_increment = _rotations[_qp](i) - _rotations_old[_qp](i);
    // std::cout<<" \n deformation increment " << deformation_increment ;
Real trial_moment = _spring_moments_global_old[_qp](i) + _krx[_qp] * rotation_increment;

_hardening_variable[_qp](i+3) = _hardening_variable_old[_qp](i+3);
// std::cout<<" \n hardening variable " << _hardening_variable[_qp] ;
_plastic_rotation[_qp](i) = _plastic_rotation_old[_qp](i);
// std::cout<<" \n plastic deformation " << _plastic_deformation[_qp] ;

Real yield_condition = std::abs(trial_moment) - _hardening_variable[_qp](i+3) - _yield_force;
// std::cout<<" \n yield condition " << yield_condition;

Real iteration = 0;
Real plastic_rotation_increment = 0.0;
Real elastic_rotation_increment = rotation_increment;
// std::cout<<" \n elastic deformation increment " << elastic_deformation_increment;

if (yield_condition > 0.0)
{
  Real residual = std::abs(trial_moment) - _hardening_variable[_qp](i+3) - _yield_force -
                  _krx[_qp] * plastic_rotation_increment;

  Real reference_residual =
      std::abs(trial_moment) - _krx[_qp] * plastic_rotation_increment;

  while (std::abs(residual) > _absolute_tolerance ||
         std::abs(residual / reference_residual) > _relative_tolerance)
  {
    _hardening_variable[_qp](i+3) = computeHardeningValue(plastic_rotation_increment,i+3);
    Real hardening_slope = computeHardeningDerivative(plastic_rotation_increment,i+3);

    Real scalar = (std::abs(trial_moment) - _hardening_variable[_qp](i+3) - _yield_force -
                   _krx[_qp] * plastic_rotation_increment) /
                  (_krx[_qp] + hardening_slope);
    plastic_rotation_increment += scalar;

    residual = std::abs(trial_moment) - _hardening_variable[_qp](i+3) - _yield_force -
               _krx[_qp] * plastic_rotation_increment;

    reference_residual = std::abs(trial_moment) - _krx[_qp] * plastic_rotation_increment;

    ++iteration;
    if (iteration > _max_its) // not converging
      throw MooseException("IsotropicSpring: Plasticity model did not converge");
  }
  plastic_rotation_increment *= MathUtils::sign(trial_moment);
  _plastic_rotation[_qp](i) += plastic_rotation_increment;
  elastic_rotation_increment = rotation_increment - plastic_rotation_increment;
}
_elastic_rotation[_qp](i) = _rotations[_qp](i) - _plastic_rotation[_qp](i); //Confused here
// std::cout<<" \n elastic deformation " << _elastic_deformation[_qp];

_spring_moments_local(i) = _spring_moments_global_old[_qp](i) + _krx[_qp] * elastic_rotation_increment;
// std::cout<<" \n spring forces local " << _spring_forces_local(0);
}


  // convert local moments to global
  _spring_moments_global[_qp] =
      _total_global_to_local_rotation[_qp].transpose() * _spring_moments_local;
}

Real
IsotropicSpring::computeHardeningValue(Real scalar, Real j)
{
  if (_hardening_function)
  {
    if(j <= 2)
    (
    const Real disp_old = _plastic_deformation_old[_qp](j);
    const Point p;

    return _hardening_function->value(std::abs(disp_old) + scalar, p) - _yield_force;
  }
    else
    {
      const Real disp_old = _plastic_rotation_old[_qp](j-3);
      const Point p;

      return _hardening_function->value(std::abs(disp_old) + scalar, p) - _yield_force;
    }

  }

  return _hardening_variable_old[_qp](j) + _hardening_constant * scalar;
}

Real IsotropicSpring::computeHardeningDerivative(Real /*scalar*/, Real j)
{
  if (_hardening_function)
  {
    if(j <= 2)
    const Real disp_old = _plastic_deformation_old[_qp](j);
    else
    const Real disp_old = _plastic_rotation_old[_qp](j-3);
    const Point p;

    return _hardening_function->timeDerivative(std::abs(disp_old), p);
  }

  return _hardening_constant;
}


void
IsotropicSpring::computeStiffnessMatrix()
{
  // The stiffness matrix is of the form
  // |  kdd  kdr  |
  // |  krd  krr  |
  // where kdd, krr, krd and kdr are all RankTwoTensors (3x3 matrix) and
  // matrix symmetry is assumed, namely, krd = kdr.transpose()
  // This implementation of the spring element has only diagonal stiffness
  // terms. Therefore, Kdr and Krd are zero.

  // calculating deformation stiffnesses
  _kdd_local(0, 0) = _kx[_qp];
  _kdd_local(1, 1) = _ky[_qp];
  _kdd_local(2, 2) = _kz[_qp];
  // convert stiffness matrix from local to global
  _kdd[_qp] = _total_global_to_local_rotation[_qp].transpose() * _kdd_local *
              _total_global_to_local_rotation[_qp];

  // calculating rotational stiffness
  _krr_local(0, 0) = _krx[_qp];
  _krr_local(1, 1) = _kry[_qp];
  _krr_local(2, 2) = _krz[_qp];
  // convert stiffness matrix from local to global
  _krr[_qp] = _total_global_to_local_rotation[_qp].transpose() * _krr_local *
              _total_global_to_local_rotation[_qp];
}
