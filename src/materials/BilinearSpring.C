// MASTODON includes
#include "BilinearSpring.h"

// MOOSE includes
#include "MooseMesh.h"
#include "Assembly.h"
#include "NonlinearSystem.h"
#include "MooseVariable.h"
#include "RankTwoTensor.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("MastodonApp", BilinearSpring);

template <>
InputParameters
validParams<BilinearSpring>()
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
  params.addRequiredParam<Real>("th_disp", "Threshold displacement after which deterioration starts");
  params.addParam<Real>("deterioration_constant", 0.0, "Deterioration slope");
  params.addParam<Real>(
     "absolute_tolerance", 1e-10, "Absolute convergence tolerance for Newton iteration");
  params.addParam<Real>(
     "relative_tolerance", 1e-8, "Relative convergence tolerance for Newton iteration");
  return params;
}

BilinearSpring::BilinearSpring(const InputParameters & parameters)
  : Material(parameters),
    _nrot(coupledComponents("rotations")),
    _ndisp(coupledComponents("displacements")),
    _rot_num(3),
    _disp_num(3),
    _deformations(declareProperty<RealVectorValue>("deformations")),
    _deformations_old(getMaterialPropertyOld<RealVectorValue>("deformations")),
    _rotations(declareProperty<RealVectorValue>("rotations")),
    _kx(coupledValue("kx")),
    _ky(coupledValue("ky")),
    _kz(coupledValue("kz")),
    _krx(coupledValue("krx")),
    _kry(coupledValue("kry")),
    _krz(coupledValue("krz")),
    _spring_forces_global(declareProperty<RealVectorValue>("global_forces")),
    _spring_forces_global_old(getMaterialPropertyOld<RealVectorValue>("global_forces")),
    _spring_forces_test(declareProperty<RealVectorValue>("test_forces")),
    _spring_forces_test_old(getMaterialPropertyOld<RealVectorValue>("test_forces")),
    _spring_moments_global(declareProperty<RealVectorValue>("global_moments")),
    _kdd(declareProperty<RankTwoTensor>("displacement_stiffness_matrix")),
    _krr(declareProperty<RankTwoTensor>("rotation_stiffness_matrix")),
    _total_global_to_local_rotation(
        declareProperty<RankTwoTensor>("total_global_to_local_rotation")),
    _yield_force(getParam<Real>("yield_force")),
    _hardening_constant(getParam<Real>("hardening_constant")),
    _th_disp(getParam<Real>("th_disp")),
    _deterioration_constant(getParam<Real>("deterioration_constant")),
    _absolute_tolerance(parameters.get<Real>("absolute_tolerance")),
    _relative_tolerance(parameters.get<Real>("relative_tolerance")),
    _plastic_deformation(declareProperty<Real>("_plastic_deformation")),
    _plastic_deformation_old(getMaterialPropertyOld<Real>("_plastic_deformation")),
    _hardening_variable(declareProperty<Real>("hardening_variable")),
    _hardening_variable_old(getMaterialPropertyOld<Real>("hardening_variable")),
    _hardening_variable_pos(declareProperty<Real>("hardening_variable_pos")),
    _hardening_variable_neg(declareProperty<Real>("hardening_variable_neg")),
    _hard_var_carry(declareProperty<Real>("hard_var_carry")),
    _deterioration_variable(declareProperty<Real>("deterioration_variable")),
    _deterioration_variable_old(getMaterialPropertyOld<Real>("deterioration_variable")),
    _hard_var_test(declareProperty<Real>("hard_var_test")),
    _hard_var_test_old(getMaterialPropertyOld<Real>("hard_var_test")),
    _det_var_test(declareProperty<Real>("det_var_test")),
    _det_var_test_old(getMaterialPropertyOld<Real>("det_var_test")),
    _eff_trial_force(declareProperty<Real>("eff_trial_force")),
    _eff_trial_force_old(getMaterialPropertyOld<Real>("eff_trial_force")),
    _elastic_deformation(declareProperty<Real>("elastic_deformation")),
    _max_its(5000),
    _yf(declareProperty<Real>("yf")),
    _yf_pos(declareProperty<Real>("yf_pos")),
    _yf_neg(declareProperty<Real>("yf_neg")),
    _th(declareProperty<Real>("th")),
    _th_disp_pos(declareProperty<Real>("th_disp_pos")),
    _th_disp_neg(declareProperty<Real>("th_disp_neg")),
    _hc(declareProperty<Real>("hc")),
    _hc_pos(declareProperty<Real>("H_pos")),
    _hc_neg(declareProperty<Real>("H_neg")),
    _shift(declareProperty<Real>("shift")),
    _max(declareProperty<Real>("max")),
    _min(declareProperty<Real>("min")),
    _isdamage(declareProperty<bool>("isdamage")),
    _isdamagepos(declareProperty<bool>("isdamagepos")),
    _isdamageneg(declareProperty<bool>("isdamageneg")),
    _isdetpos(declareProperty<bool>("isdetpos")),
    _isdetneg(declareProperty<bool>("isdetneg")),
    _istest(declareProperty<bool>("istest")),
    _istest2(declareProperty<bool>("istest2")),
    _istest3pos(declareProperty<bool>("istest3pos")),
    _istest3neg(declareProperty<bool>("istest3neg")),
    _istemp(declareProperty<bool>("istemp"))
{
  // Checking for consistency between length of the provided displacements and rotations vector
  if (_ndisp != _nrot)
    mooseError("BilinearSpring: The number of variables supplied in 'displacements' "
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
BilinearSpring::initQpStatefulProperties()
{
  _deformations[_qp].zero();
  _spring_forces_global[_qp].zero();
  _spring_forces_test[_qp].zero();
  _elastic_deformation[_qp] = 0.0;
  _yf[_qp]= _yield_force;
  _yf_neg[_qp] = _yield_force;
  _yf_pos[_qp] = _yield_force;
  _th[_qp] = _th_disp;
  _th_disp_pos[_qp] = _th_disp;
  _th_disp_neg[_qp] = _th_disp;
  _hc[_qp] = _hardening_constant;
  _hc_pos[_qp] = _hardening_constant;
  _hc_neg[_qp] = _hardening_constant;
  _shift[_qp] = 0;
  _max[_qp] = 0;
  _min[_qp] = 0;
  _isdamage[_qp] = false;
  _isdamagepos[_qp] = false;
  _isdamageneg[_qp] = false;
  _isdetpos[_qp] = false;
  _isdetneg[_qp] = false;
  _istest[_qp] = false;
  _istest2[_qp] = false;
  _istest3pos[_qp] = false;
  _istest3neg[_qp] = false;
  _istemp[_qp] = false;

}

void
BilinearSpring::computeQpProperties()
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
    mooseError("Error in BilinearSpring: y_orientation should be perpendicular to "
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
BilinearSpring::computeTotalRotation()
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
BilinearSpring::computeDeformations()
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

    std::cout<<" \n _global_disp0 x " << _global_disp0(0);
     std::cout<<" \n _global_disp1 x " << _global_disp1(0);
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
     std::cout<<" \n _local_disp0 x " << _local_disp0(0);
      std::cout<<" \n _local_disp1 x " << _local_disp1(0);
  _rotations[_qp] = _local_rot1 - _local_rot0;
}

void
BilinearSpring::computeForces()
// spring forces = Kdd * deformations
// spring moments = Krr * rotations
{
  std::cout<<" \n deformation " << _deformations[_qp](0);
  std::cout<<" \n deformation(i-1) " << _deformations_old[_qp](0);

  if ( _deformations[_qp](0) >0 && _deformations[_qp](0)< _deformations_old[_qp](0) && (-1* _deformations[_qp](0)) > (-1*_deformations_old[_qp](0)))
  {
    if (_deformations_old[_qp](0) > _max[_qp])
    {
      _max[_qp] = _deformations_old[_qp](0);
    }
  }

  if ( _deformations[_qp](0) <0 && _deformations[_qp](0)> _deformations_old[_qp](0) && -1* _deformations[_qp](0) < -1*_deformations_old[_qp](0))
  {
    if (_deformations_old[_qp](0) < _min[_qp])
    {
      _min[_qp] = _deformations_old[_qp](0);
    }
  }

  if ( _isdetpos[_qp] == false)
  {
      std::cout<<" \n deformation outside if condition " << _deformations[_qp](0);
        std::cout<<" \n threshold disp outside if condition " << _th_disp;
    if ( _deformations[_qp](0) > ( _th_disp ))
    {
      std::cout<<" \n deformation inside if condition " << _deformations[_qp](0);
        std::cout<<" \n threshold disp inside if condition " << _th_disp;
      _isdetpos[_qp] = true;
    }
  }

  if ( _isdetneg[_qp] == false)
  {
    if ( _deformations[_qp](0) < _th_disp && std::abs( _deformations[_qp](0)) > std::abs(_th_disp))
    {
      _isdetneg[_qp] = true;
    }
  }

  Real deformation_increment = _deformations[_qp](0) - _deformations_old[_qp](0);
     // std::cout<<" \n deformation increment " << deformation_increment ;
  Real trial_force = _spring_forces_global_old[_qp](0) + _kx[_qp] * deformation_increment;
     // std::cout<<" \n trial force " << trial_force ;
  _hardening_variable[_qp]= _hardening_variable_old[_qp];
  Real deterioration_var_old = _deterioration_variable_old[_qp];
  _deterioration_variable[_qp] = _deterioration_variable_old[_qp];
  _hard_var_test[_qp] = _hard_var_test_old[_qp];
  _det_var_test[_qp] = _det_var_test_old[_qp];
  _spring_forces_test[_qp](0) = _spring_forces_test_old[_qp](0);

  Real bt = _hardening_variable_old[_qp];


  if(((deterioration_var_old == 0) && (_isdamage[_qp]==false)) ||( (std::abs(_deformations[_qp](0))) < (std::abs(_deformations_old[_qp](0))) ) )
  {
    if(  ((std::abs(_deformations[_qp](0))) < (std::abs(_deformations_old[_qp](0)))) && _istest[_qp]==true  )
    {

      if( _isdetpos[_qp] ==true && _isdetneg[_qp] ==false && _deformations[_qp](0) < _deformations_old[_qp](0))
      {
        _hardening_variable[_qp] = _hard_var_test[_qp];
        deterioration_var_old = _det_var_test[_qp];
        _istest[_qp] = false;
        _hard_var_test[_qp] = 0.0;
        _det_var_test[_qp] = 0.0;
        if( _deformations_old[_qp](0) < _max[_qp])
        {
          _shift[_qp] = std::abs(_spring_forces_test_old[_qp](0) - _spring_forces_global_old[_qp](0));
          _istest2[_qp] = true;
        }
      }

      if( _isdetpos[_qp] == false && _isdetneg[_qp] == true && _deformations[_qp](0) > _deformations_old[_qp](0))
      {
        _hardening_variable[_qp] = _hard_var_test[_qp];
        deterioration_var_old = _det_var_test[_qp];
        _istest[_qp] = false;
        _hard_var_test[_qp] = 0.0;
        _det_var_test[_qp] = 0.0;
        if( _deformations_old[_qp](0) > _min[_qp])
        {
          _shift[_qp] = std::abs(_spring_forces_test[_qp](0) - _spring_forces_global_old[_qp](0));
          _istest2[_qp] = true;
        }
      }
    }

    if ( _isdamage[_qp] == true && _deformations[_qp](0) > _deformations_old[_qp](0))
    {
      _th_disp_neg[_qp] = std::abs(_deformations_old[_qp](0));
      _hc_neg[_qp] = 0.0;
      _hardening_variable_pos[_qp] = _hardening_variable[_qp];
      _hardening_variable_neg[_qp] = 0;
      _isdamageneg[_qp] = true ;
      _yf_neg[_qp] = std::abs(_spring_forces_global_old[_qp](0));
    }

    if ( _isdamage[_qp] == true && _deformations[_qp](0) < _deformations_old[_qp](0))
    {
      _th_disp_pos[_qp] = std::abs(_deformations_old[_qp](0));
      _hc_pos[_qp] = 0.0;
      _hardening_variable_neg[_qp] = _hardening_variable[_qp];
      _hardening_variable_pos[_qp] = 0;
      _isdamagepos[_qp] = true ;
      _yf_pos[_qp] = std::abs(_spring_forces_global_old[_qp](0));
    }

    deterioration_var_old = 0.0;
    bt = _hardening_variable[_qp];
    _eff_trial_force[_qp] = trial_force - _hardening_variable[_qp];

  }
  else
  {
    _eff_trial_force[_qp] = trial_force - deterioration_var_old;
  }

  if(_eff_trial_force[_qp] >=0 && _eff_trial_force_old[_qp] <0 )
  {
    if (_isdetpos[_qp] == true && _isdetneg[_qp] == false && _istemp[_qp] == false)
    {
      _hard_var_carry[_qp] = _hardening_variable[_qp] + _shift[_qp];
      _istemp[_qp] = true;
    }
    _th[_qp] = _th_disp_pos[_qp];
    _hc[_qp] = _hc_pos[_qp];
    if(_istest2[_qp])
    {
      _yf_pos[_qp] = _yf[_qp] - _shift[_qp];
      _istest2[_qp] = false;
      _istest3neg[_qp] = true;
    }
    _yf[_qp] = _yf_pos[_qp];
    if(_isdamage[_qp] ==true || _isdamageneg[_qp]==true || _istest3pos[_qp] == true || _isdetpos[_qp]==true)
    {
      if(_isdetpos[_qp] == true || _istest3pos[_qp] == true)
      {
        _hardening_variable[_qp] = 0;
        _istest3pos[_qp] = false;
      }
      else
      {
        _hardening_variable[_qp] = _hardening_variable_pos[_qp];
      }
    _isdamageneg[_qp] = false;
    }
    if(_hardening_variable[_qp] != bt)
    {
      _eff_trial_force[_qp] = trial_force - _hardening_variable[_qp];
    }
  }

  if(_eff_trial_force[_qp] <0 && _eff_trial_force_old[_qp] >=0 )
  {
    if (_isdetpos[_qp] == false && _isdetneg[_qp] == true && _istemp[_qp] == false)
    {
      _hard_var_carry[_qp] = _hardening_variable[_qp] - _shift[_qp];
      _istemp[_qp] = true;
    }
    _th[_qp] = _th_disp_neg[_qp];
    _hc[_qp] = _hc_neg[_qp];
    if(_istest2[_qp])
    {
      _yf_neg[_qp] = _yf[_qp] - _shift[_qp];
      _istest2[_qp] = false;
      _istest3pos[_qp] = true;
    }
    _yf[_qp] = _yf_neg[_qp];
    if(_isdamage[_qp] ==true || _isdamagepos[_qp]==true || _istest3neg[_qp] == true || _isdetneg[_qp]==true)
    {
      if(_isdetneg[_qp] == true || _istest3neg[_qp] == true)
      {
        _hardening_variable[_qp] = 0;
        _istest3neg[_qp] = false;
      }
      else
      {
        _hardening_variable[_qp] = _hardening_variable_neg[_qp];
      }
    _isdamagepos[_qp] = false;
    }
    if(_hardening_variable[_qp] != bt)
    {
      _eff_trial_force[_qp] = trial_force - _hardening_variable[_qp];
    }

  }

  Real yield_condition = std::abs(_eff_trial_force[_qp]) - _yf[_qp];
  //std::cout<<" \n yield condition " << yield_condition;

  Real iteration = 0;
  Real plastic_deformation_increment = 0.0;
  Real elastic_deformation_increment = deformation_increment;

  if (yield_condition > 0.0)
  {
    if(std::abs(_deformations[_qp](0)) <= (_th[_qp] ))
    {
      if((_eff_trial_force[_qp] >=0 && _isdetneg[_qp] == false && _isdetpos[_qp] ==true ) || (_eff_trial_force[_qp] <=0 && _isdetneg[_qp] == true && _isdetpos[_qp] ==false ))
      {
        Real trial_force_test = 0.0;
        Real eff_trial_test = 0.0;

        if(_hard_var_test_old[_qp] ==0)
        {
          _hard_var_test[_qp] = _hard_var_carry[_qp];
          _hard_var_carry[_qp] = 0;
          _istemp[_qp] = false;
          trial_force_test = trial_force;
          _spring_forces_test[_qp](0) = _spring_forces_global_old[_qp](0);
        }
        else
        {
          _hard_var_test[_qp] = _hard_var_test_old[_qp];
          _spring_forces_test[_qp](0) = _spring_forces_test_old[_qp](0);
          trial_force_test = _spring_forces_test[_qp](0) + _kx[_qp] * deformation_increment ;
        }
        std::cout<<" \n trial_stress_test " << trial_force_test;

        if(std::abs(_deformations[_qp](0)) <= _th_disp)
        {
          eff_trial_test = trial_force_test - _hard_var_test[_qp];
        }
        else
          eff_trial_test = trial_force_test - _det_var_test[_qp];

        Real yield_condn_test = abs(eff_trial_test) - _yield_force;
        Real plastic_inc_test = 0.0;
        Real h_test = _hardening_constant;
        Real elastic_inc_test = deformation_increment;

        if(yield_condn_test >0)
        {
          if(std::abs(_deformations[_qp](0)) <= _th_disp)
          {
            Real residual_test = std::abs(eff_trial_test) - _yield_force -(_kx[_qp] + h_test) * plastic_inc_test;
            Real iteration = 0;
            while (std::abs(residual_test) > _absolute_tolerance)
            {
              Real scalar_test = (std::abs(eff_trial_test) - _yield_force - (_kx[_qp] + h_test) * plastic_inc_test ) / (_kx[_qp] + h_test);
              plastic_inc_test += scalar_test;
              residual_test = std::abs(eff_trial_test) - _yield_force -(_kx[_qp] + h_test) * plastic_inc_test;
              ++iteration;
              if (iteration > _max_its) // not converging
              throw MooseException("BilinearSpring: Plasticity model did not converge");
            }
            plastic_inc_test *= MathUtils::sign(eff_trial_test);
            elastic_inc_test = deformation_increment - plastic_inc_test;
            _hard_var_test[_qp] = _hard_var_test[_qp] + h_test * plastic_inc_test;
            _spring_forces_test[_qp](0) = _spring_forces_test[_qp](0) + _kx[_qp] * elastic_inc_test;
            _istest[_qp] = true;
          }

          if(std::abs(_deformations[_qp](0)) > _th_disp)
          {
            if(_det_var_test[_qp] == 0 && _isdamage[_qp] == false)
            {
              _det_var_test[_qp] = _hard_var_test[_qp];
            }
            else
            {
              _det_var_test[_qp] = _det_var_test[_qp];
            }
            Real residual_test = std::abs(eff_trial_test) - _yield_force -(_kx[_qp] + _deterioration_constant) * plastic_inc_test;
            Real iteration = 0;
            while (std::abs(residual_test) > _absolute_tolerance)
            {
              Real scalar_test = (std::abs(eff_trial_test) - _yield_force - (_kx[_qp] + _deterioration_constant) * plastic_inc_test ) / (_kx[_qp] + _deterioration_constant);
              plastic_inc_test += scalar_test;
              residual_test = std::abs(eff_trial_test) - _yield_force -(_kx[_qp] + _deterioration_constant) * plastic_inc_test;
              ++iteration;
              if (iteration > _max_its) // not converging
              throw MooseException("BilinearSpring: Plasticity model did not converge");
            }
            plastic_inc_test *= MathUtils::sign(eff_trial_test);
            elastic_inc_test = deformation_increment - plastic_inc_test;
            _hard_var_test[_qp] = _hard_var_test[_qp] + h_test * plastic_inc_test;
            _det_var_test[_qp] = _det_var_test[_qp] + _deterioration_constant * plastic_inc_test;
            _spring_forces_test[_qp](0) = _spring_forces_test[_qp](0) + _kx[_qp] * elastic_inc_test;
            _istest[_qp] = true;
          }
        }
        else
        {
          elastic_inc_test = deformation_increment;
          _spring_forces_test[_qp](0) = _spring_forces_test[_qp](0) + _kx[_qp] * elastic_inc_test;
        }
      }

      Real residual = std::abs(trial_force - _hardening_variable[_qp]) - _yf[_qp] - ( _kx[_qp] + _hc[_qp])* plastic_deformation_increment;
      while (std::abs(residual) > _absolute_tolerance)
      {
        Real scalar = (std::abs(trial_force - _hardening_variable[_qp]) -
                      (_kx[_qp] + _hc[_qp]) * plastic_deformation_increment  - _yf[_qp] ) /
                      (_kx[_qp] + _hc[_qp]);
        plastic_deformation_increment += scalar;
        residual = std::abs(trial_force - _hardening_variable[_qp]) - _yf[_qp] - ( _kx[_qp] + _hc[_qp])* plastic_deformation_increment;
        ++iteration;
        if (iteration > _max_its) // not converging
        throw MooseException("BilinearSpring: Plasticity model did not converge");
      }
      plastic_deformation_increment *= MathUtils::sign(_eff_trial_force[_qp]);
      _hardening_variable[_qp] = _hardening_variable[_qp] + _hc[_qp] * plastic_deformation_increment;
      _plastic_deformation[_qp] += plastic_deformation_increment;
      elastic_deformation_increment = deformation_increment - plastic_deformation_increment;
      _deterioration_variable[_qp] = 0;
      std::cout<<" \n back trial_d in pgm1 " << _deterioration_variable[_qp] ;
      _spring_forces_local(0) = _spring_forces_global_old[_qp](0) +  _kx[_qp] * elastic_deformation_increment;
    }

    std::cout<<" \n th in pgm " << _th[_qp] ;
    Real ab= std::abs(_deformations[_qp](0));
    std::cout<<" \n def pgm " << ab ;

    if((std::abs(_deformations[_qp](0))) > ( _th[_qp]  ))
    {
      std::cout<<" \n th in pgm2 " << _th[_qp] ;
      Real ab= std::abs(_deformations[_qp](0));
      std::cout<<" \n def pgm2 " << ab ;
      if( deterioration_var_old == 0 && _isdamage[_qp] == false)
      {
        deterioration_var_old = _hardening_variable[_qp];
      }
      else
      {
        deterioration_var_old = deterioration_var_old;
      }
//start here
      Real residual = std::abs(trial_force - deterioration_var_old) - _yf[_qp] - ( _kx[_qp] + _deterioration_constant)* plastic_deformation_increment;
      while (std::abs(residual) > _absolute_tolerance)
      {
        Real scalar = (std::abs(trial_force - deterioration_var_old) -
                      (_kx[_qp] + _deterioration_constant) * plastic_deformation_increment  - _yf[_qp] ) /
                      (_kx[_qp] + _deterioration_constant);
        plastic_deformation_increment += scalar;
        residual = std::abs(trial_force - deterioration_var_old) - _yf[_qp] - ( _kx[_qp] + _deterioration_constant)* plastic_deformation_increment;
        ++iteration;
        if (iteration > _max_its) // not converging
        throw MooseException("BilinearSpring: Plasticity model did not converge");
      }
      plastic_deformation_increment *= MathUtils::sign(_eff_trial_force[_qp]);
      _plastic_deformation[_qp] += plastic_deformation_increment;
      elastic_deformation_increment = deformation_increment - plastic_deformation_increment;
      _spring_forces_local(0) = _spring_forces_global_old[_qp](0) +  _kx[_qp] * elastic_deformation_increment;

      if(std::abs(_hard_var_test[_qp]) > _absolute_tolerance && (( _isdetneg[_qp] == false && _isdetpos[_qp] ==true ) ||( _isdetneg[_qp] == true && _isdetpos[_qp] ==false ) ))
      {
        _hard_var_test[_qp] = _hard_var_test[_qp] + _hardening_constant * plastic_deformation_increment;
        _det_var_test[_qp] = _det_var_test[_qp] + _deterioration_constant * plastic_deformation_increment;
        _istest[_qp] = true;
      }
      _deterioration_variable[_qp] = deterioration_var_old + _deterioration_constant * plastic_deformation_increment;
      std::cout<<" \n back trial_d in pgm2 " << _deterioration_variable[_qp] ;
      _hardening_variable[_qp] = _hardening_variable[_qp] + _hc[_qp] * plastic_deformation_increment;
    }
  }
  else
  {
    _spring_forces_local(0) = _spring_forces_global_old[_qp](0) +  _kx[_qp] * elastic_deformation_increment;

    _deterioration_variable[_qp] = 0;
  }



  if(_deterioration_variable[_qp] == 0 && deterioration_var_old ==0)
  {
    _isdamage[_qp] = false;
  }
  else
  {
    _isdamage[_qp] =true;
  }


  std::cout<<" \n back trial " << _hardening_variable[_qp] ;
  std::cout<<" \n back trial_d " << _deterioration_variable[_qp] ;
  std::cout<<" \n back trial test " << _hard_var_test[_qp] ;
  std::cout<<" \n back trial_d test " << _det_var_test[_qp] ;
  std::cout<<" \n back trial pos " << _hardening_variable_pos[_qp] ;
  std::cout<<" \n back trial neg" << _hardening_variable_neg[_qp] ;
  std::cout<<" \n bt " << bt;
  std::cout<<" \n dam_neg " << _isdamageneg[_qp];
  std::cout<<" \n dam_pos " << _isdamagepos[_qp];
  std::cout<<" \n damage " << _isdamage[_qp];
  std::cout<<" \n eff trial stress " << _eff_trial_force[_qp];
  std::cout<<" \n elastic_inc " << elastic_deformation_increment;
  std::cout<<" \n H " << _hc[_qp];
  std::cout<<" \n H_neg " << _hc_neg[_qp];
  std::cout<<" \n H_pos " << _hc_pos[_qp];
  std::cout<<" \n max " << _max[_qp];
  std::cout<<" \n min " << _min[_qp];
  std::cout<<" \n neg_det " << _isdetneg[_qp];
  std::cout<<" \n pl_disp " << _plastic_deformation[_qp];
  std::cout<<" \n plastic_inc " << plastic_deformation_increment;
  std::cout<<" \n pos_det " << _isdetpos[_qp];
  std::cout<<" \n shift " << _shift[_qp];
  std::cout<<" \n spring forces local " << _spring_forces_local(0);
  std::cout<<" \n temp " << _istemp[_qp];
  std::cout<<" \n test " << _istest[_qp];
  std::cout<<" \n test2 " << _istest2[_qp];
  std::cout<<" \n test3pos " << _istest3pos[_qp];
  std::cout<<" \n test3neg " << _istest3neg[_qp];
  std::cout<<" \n th_disp " << _th[_qp];
  std::cout<<" \n th_disp_neg " << _th_disp_neg[_qp];
  std::cout<<" \n th_disp_pos " << _th_disp_pos[_qp];
  std::cout<<" \n trial_stress " << trial_force;
  // std::cout<<" \n trial_stress_test " << trial_force_test;
  std::cout<<" \n y0 " << _yf[_qp];
  std::cout<<" \n y0_neg " << _yf_neg[_qp];
  std::cout<<" \n y0_pos " << _yf_pos[_qp];
  std::cout<<" \n yield_condition " << yield_condition;


  _spring_forces_local(1) = _ky[_qp] * _deformations[_qp](1);
  _spring_forces_local(2) = _kz[_qp] * _deformations[_qp](2);
  // convert local forces to global
  _spring_forces_global[_qp] =
      _total_global_to_local_rotation[_qp].transpose() * _spring_forces_local;

  // moments
  _spring_moments_local(0) = _krx[_qp] * _rotations[_qp](0);
  _spring_moments_local(1) = _kry[_qp] * _rotations[_qp](1);
  _spring_moments_local(2) = _krz[_qp] * _rotations[_qp](2);
  // convert local moments to global
  _spring_moments_global[_qp] =
      _total_global_to_local_rotation[_qp].transpose() * _spring_moments_local;
}

// Real
// BilinearSpring::computeHardeningValue(Real scalar)
// {
//
//
//   // return _hardening_variable_old[_qp] + MathUtils::sign(back)* _hardening_constant * scalar;
//   return _hardening_variable_old[_qp] + _hardening_constant * scalar;
// }



void
BilinearSpring::computeStiffnessMatrix()
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
