// MASTODON includes
#include "KinematicSpring.h"

// MOOSE includes
#include "MooseMesh.h"
#include "Assembly.h"
#include "NonlinearSystem.h"
#include "MooseVariable.h"
#include "RankTwoTensor.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("MastodonApp", KinematicSpring);


template <>
InputParameters
validParams<KinematicSpring>()
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
  params.addRequiredParam<RealVectorValue>("yield_force",
                               "Yield force after which plastic strain starts accumulating");
  params.addRequiredParam<RealVectorValue>("yield_moments",
                              "Yield moments after which plastic strain starts accumulating");
 params.addRequiredParam<RealVectorValue>("hardening_constant_force", "Hardening slope for forces");
 params.addRequiredParam<RealVectorValue>("hardening_constant_moment", "Hardening slope for moments");
 params.addParam<Real>(
     "absolute_tolerance", 1e-8, "Absolute convergence tolerance for Newton iteration");
 return params;
}

KinematicSpring::KinematicSpring(const InputParameters & parameters)
  : Material(parameters),
    _nrot(coupledComponents("rotations")),
    _ndisp(coupledComponents("displacements")),
    _rot_num(3),
    _disp_num(3),
    _deformations(declareProperty<RealVectorValue>("deformations")),
    _deformations_old(getMaterialPropertyOld<RealVectorValue>("deformations")),
    _rotations(declareProperty<RealVectorValue>("rotations")),
    _rotations_old(getMaterialPropertyOld<RealVectorValue>("rotations")),
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
    _total_global_to_local_rotation(declareProperty<RankTwoTensor>("total_global_to_local_rotation")),
    _yield_force(getParam<RealVectorValue>("yield_force")),
    _yield_moments(getParam<RealVectorValue>("yield_moments")),
    _hardening_constant_force(getParam<RealVectorValue>("hardening_constant_force")),
    _hardening_constant_moment(getParam<RealVectorValue>("hardening_constant_moment")),
    _absolute_tolerance(parameters.get<Real>("absolute_tolerance")),
    _plastic_deformation(declareProperty<RealVectorValue>("_plastic_deformation")),
    _plastic_deformation_old(getMaterialPropertyOld<RealVectorValue>("_plastic_deformation")),
    _plastic_rotation(declareProperty<RealVectorValue>("_plastic_rotation")),
    _plastic_rotation_old(getMaterialPropertyOld<RealVectorValue>("_plastic_rotation")),
    _hardening_variable(declareProperty<RealVectorValue>("hardening_variable")),
    _hardening_variable_old(getMaterialPropertyOld<RealVectorValue>("hardening_variable")),
    _hardening_variable_moment(declareProperty<RealVectorValue>("hardening_variable_moment")),
    _hardening_variable_moment_old(getMaterialPropertyOld<RealVectorValue>("hardening_variable_moment")),
    _elastic_deformation(declareProperty<RealVectorValue>("elastic_deformation")),
    _elastic_rotation(declareProperty<RealVectorValue>("elastic_rotation")),
    _max_its(5000)
{
  // Checking for consistency between length of the provided displacements and rotations vector
  if (_ndisp != _nrot)
    mooseError("KinematicSpring: The number of variables supplied in 'displacements' "
               "and 'rotations' input parameters must be equal.");

  // Fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    MooseVariable * disp_variable = getVar("displacements", i);
    _disp_num[i] = disp_variable->number(); // Displacement variable numbers in MOOSE

    MooseVariable * rot_variable = getVar("rotations", i);
    _rot_num[i] = rot_variable->number(); // Rotation variable numbers in MOOSE
  }
}

void
KinematicSpring::initQpStatefulProperties()
{
  _deformations[_qp].zero();
  _spring_forces_global[_qp].zero();
  _elastic_deformation[_qp] = 0.0;
  _rotations[_qp].zero();
  _spring_moments_global[_qp].zero();
  _elastic_rotation[_qp] = 0.0;
}

void
KinematicSpring::computeQpProperties()
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
    mooseError("Error in KinematicSpring: y_orientation should be perpendicular to "
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
KinematicSpring::computeTotalRotation()
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
KinematicSpring::computeDeformations()
{
  // fetch the two end nodes for _current_elem
  std::vector<const Node *> node;
  for (unsigned int i = 0; i < 2; ++i)
    node.push_back(_current_elem->node_ptr(i));

  // Fetch the solution for the two end nodes at time t
  NonlinearSystemBase & nonlinear_sys = _fe_problem.getNonlinearSystemBase();
  const NumericVector<Number> & sol = *nonlinear_sys.currentSolution();


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
  _rotations[_qp] = _local_rot1 - _local_rot0;
}

void
KinematicSpring::computeForces()
{
for (unsigned int i =0; i < _ndisp; i++)
    {
      ///Calculate trial displacment
      Real deformation_increment = _deformations[_qp](i) - _deformations_old[_qp](i);

      if(i==0)
      _k = _kx[_qp];
      if(i==1)
      _k = _ky[_qp];
      if(i==2)
      _k = _kz[_qp];

      ///Calculate trial force
      Real trial_force = _spring_forces_global_old[_qp](i) + _k * deformation_increment;

      _hardening_variable[_qp](i) = _hardening_variable_old[_qp](i);
      _plastic_deformation[_qp](i) = _plastic_deformation_old[_qp](i);

      ///Calculate effective trial force
      Real effective_force = trial_force - _hardening_variable_old[_qp](i);

      /// Calculate yield condition
      Real yield_condition = std::abs(trial_force - _hardening_variable[_qp](i))  - _yield_force(i);

      Real iteration = 0;
      Real plastic_deformation_increment = 0.0;

      /// Assume all deformation as elastic deformation
      Real elastic_deformation_increment = deformation_increment;

      ///Check yield condition i.e. if positive plastic deformation occurs else only elastic deformation
      if (yield_condition > 0.0)
          {
            ///Calculate residual
            Real residual = std::abs(trial_force - _hardening_variable_old[_qp](i))-
                  (_k + _hardening_constant_force(i)) * plastic_deformation_increment  - _yield_force(i);

            /// Determine the plastic deformation required to bring the trial force to yield surface
            while (std::abs(residual) > _absolute_tolerance)
              {

                Real scalar = (std::abs(trial_force - _hardening_variable_old[_qp](i)) -
                              (_k + _hardening_constant_force(i)) * plastic_deformation_increment  - _yield_force(i) ) /
                              (_k + _hardening_constant_force(i));
                plastic_deformation_increment += scalar;
                residual = std::abs(trial_force - _hardening_variable_old[_qp](i))-
                           (_k + _hardening_constant_force(i)) * plastic_deformation_increment  - _yield_force(i);
                ++iteration;
                if (iteration > _max_its) // not converging
                throw MooseException("KinematicModel: Plasticity model did not converge");
              }
            plastic_deformation_increment *= MathUtils::sign(effective_force);
            _hardening_variable[_qp](i) = computeHardeningValue(plastic_deformation_increment,i);
            _plastic_deformation[_qp](i) += plastic_deformation_increment;
            elastic_deformation_increment = deformation_increment - plastic_deformation_increment;
           }
      _elastic_deformation[_qp](i) = _deformations[_qp](i) - _plastic_deformation[_qp](i);

      /// Calculate the spring forces
      _spring_forces_local(i) = _spring_forces_global_old[_qp](i) + _k * elastic_deformation_increment;
    }


  // convert local forces to global
  _spring_forces_global[_qp] =
      _total_global_to_local_rotation[_qp].transpose() * _spring_forces_local;

  // moments
  for (unsigned int i =0; i < 3; i++)
      {

        /// Calculate rotation increment
        Real rotation_increment = _rotations[_qp](i) - _rotations_old[_qp](i);

        if(i==0)
        _k = _krx[_qp];
        if(i==1)
        _k = _kry[_qp];
        if(i==2)
        _k = _krz[_qp];


        Real trial_moment = _spring_moments_global_old[_qp](i) + _k * rotation_increment;

        _hardening_variable_moment[_qp](i) = _hardening_variable_moment_old[_qp](i);
        _plastic_rotation[_qp](i) = _plastic_rotation_old[_qp](i);

        Real effective_moment = trial_moment - _hardening_variable_moment_old[_qp](i);

        Real yield_condition = std::abs(trial_moment - _hardening_variable_moment[_qp](i))  - _yield_moments(i);

        Real iteration = 0;
        Real plastic_rotation_increment = 0.0;
        Real elastic_rotation_increment = rotation_increment;


        if (yield_condition > 0.0)
            {
              Real residual = std::abs(trial_moment - _hardening_variable_moment_old[_qp](i))-
                    (_k + _hardening_constant_moment(i)) * plastic_rotation_increment  - _yield_moments(i);

              while (std::abs(residual) > _absolute_tolerance)
                {

                  Real scalar = (std::abs(trial_moment - _hardening_variable_moment_old[_qp](i)) -
                                (_k + _hardening_constant_moment(i)) * plastic_rotation_increment  - _yield_moments(i) ) /
                                (_k + _hardening_constant_moment(i));
                  plastic_rotation_increment += scalar;
                  residual = std::abs(trial_moment - _hardening_variable_moment_old[_qp](i))-
                             (_k + _hardening_constant_moment(i)) * plastic_rotation_increment  - _yield_moments(i);
                  ++iteration;
                  if (iteration > _max_its) // if not converging
                  throw MooseException("KinematicModel: Plasticity model did not converge");
                }
              plastic_rotation_increment *= MathUtils::sign(effective_moment);
              _hardening_variable_moment[_qp](i) = computeHardeningMomentValue(plastic_rotation_increment,i);
              _plastic_rotation[_qp](i) += plastic_rotation_increment;
              elastic_rotation_increment = rotation_increment - plastic_rotation_increment;
             }
        _elastic_rotation[_qp](i) = _rotations[_qp](i) - _plastic_rotation[_qp](i);

        _spring_moments_local(i) = _spring_moments_global_old[_qp](i) + _k * elastic_rotation_increment;
      }

  // convert local moments to global
  _spring_moments_global[_qp] =
      _total_global_to_local_rotation[_qp].transpose() * _spring_moments_local;
}

Real
KinematicSpring::computeHardeningValue(Real scalar, Real j)
{
  return _hardening_variable_old[_qp](j) + _hardening_constant_force(j) * scalar;
}

Real
KinematicSpring::computeHardeningMomentValue(Real scalar, Real j)
{
  return _hardening_variable_moment_old[_qp](j) + _hardening_constant_moment(j) * scalar;
}




void
KinematicSpring::computeStiffnessMatrix()
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
