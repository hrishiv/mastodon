// MASTODON includes
#include "FuelCompressionSpring.h"

// MOOSE includes
#include "MooseMesh.h"
#include "Assembly.h"
#include "NonlinearSystem.h"
#include "MooseVariable.h"
#include "RankTwoTensor.h"
#include "MathUtils.h"
#include "MooseUtils.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("MastodonApp", FuelCompressionSpring);

template <>
InputParameters
validParams<FuelCompressionSpring>()
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
  params.addRequiredCoupledVar("k_trans", "Translational stiffness of the spring");
  params.set<MooseEnum>("constant_on") = "ELEMENT"; // set _qp to 0
  params.addRequiredParam<Real>("yield_force",
                               "Yield force after which plastic strain starts accumulating");
  params.addParam<Real>("hardening_constant" ,0, "Hardening slope");
  params.addParam<Real>("alpha_s", "Hardening factor");
  params.addParam<Real>("recovery_factor",0, "Recovery factor");
  params.addRequiredParam<Real>("th_disp", "Threshold displacement after which deterioration starts");
  params.addParam<Real>("deterioration_constant",0, "Deterioration slope");
  params.addParam<Real>("alpha_c","Post capping stiffness  factor");
  params.addParam<bool>("Mahin_mod", false , "Mahin and Bertero's modification parameter");
  params.addParam<bool>("complete_det", false , "Set to true if you want spring deterioration with"
                        " the spacer grid stiffness and zero recovery after it flattens out");
  params.addRequiredParam<Real>("component", "Translational/Rotational Direction 0,1,2,3,4,5");
  params.addRequiredParam<Real>("k_sg", "Slope after the spring flattens");
  params.addParam<Real>("tension_factor",1e-7, "Factor for tension stiffness");
  params.addParam<Real>(
     "absolute_tolerance", 1e-10, "Absolute convergence tolerance for Newton iteration");
  params.addParam<Real>(
     "relative_tolerance", 1e-8, "Relative convergence tolerance for Newton iteration");
  return params;
}

FuelCompressionSpring::FuelCompressionSpring(const InputParameters & parameters)
  : Material(parameters),
    _nrot(coupledComponents("rotations")),
    _ndisp(coupledComponents("displacements")),
    _rot_num(3),
    _disp_num(3),
    _deformations(declareProperty<RealVectorValue>("deformations")),
    _deformations_old(getMaterialPropertyOld<RealVectorValue>("deformations")),
    _deformations_older(getMaterialPropertyOlder<RealVectorValue>("deformations")),
    _rotations(declareProperty<RealVectorValue>("rotations")),
    _k_trans(coupledValue("k_trans")),
    _spring_forces_local(declareProperty<RealVectorValue>("local_spring_forces")),
    _spring_forces_local_old(getMaterialPropertyOld<RealVectorValue>("local_spring_forces")),
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
    _recovery_factor(getParam<Real>("recovery_factor")),
    _th_disp(getParam<Real>("th_disp")),
    _k_sg(getParam<Real>("k_sg")),
    _Mahin_mod(getParam<bool>("Mahin_mod")),
    _complete_det(getParam<bool>("complete_det")),
    _component(getParam<Real>("component")),
    _tension_factor(getParam<Real>("tension_factor")),
    _absolute_tolerance(parameters.get<Real>("absolute_tolerance")),
    _relative_tolerance(parameters.get<Real>("relative_tolerance")),
    _plastic_deformation(declareProperty<Real>("_plastic_deformation")),
    _plastic_deformation_old(getMaterialPropertyOld<Real>("_plastic_deformation")),
    _hardening_variable(declareProperty<Real>("hardening_variable")),
    _hardening_variable_old(getMaterialPropertyOld<Real>("hardening_variable")),
    _deterioration_variable(declareProperty<Real>("deterioration_variable")),
    _deterioration_variable_old(getMaterialPropertyOld<Real>("deterioration_variable")),
    _net_def(declareProperty<Real>("net_def")),
    _net_def_old(getMaterialPropertyOld<Real>("net_def")),
    _elastic_deformation(declareProperty<Real>("elastic_deformation")),
    _max_its(5000),
    _min_c(0.0),
    _k(declareProperty<Real>("k")),
    _k_old(getMaterialPropertyOld<Real>("k")),
    _yf(declareProperty<Real>("yf")),
    _yf_old(getMaterialPropertyOld<Real>("yf")),
    _yf_pos(declareProperty<Real>("yf_pos")),
    _yf_pos_old(getMaterialPropertyOld<Real>("yf_pos")),
    _yf_neg(declareProperty<Real>("yf_neg")),
    _yf_neg_old(getMaterialPropertyOld<Real>("yf_neg")),
    _th(declareProperty<Real>("th")),
    _th_old(getMaterialPropertyOld<Real>("th")),
    _th_disp_pos(declareProperty<Real>("th_disp_pos")),
    _th_disp_pos_old(getMaterialPropertyOld<Real>("th_disp_pos")),
    _th_disp_neg(declareProperty<Real>("th_disp_neg")),
    _th_disp_neg_old(getMaterialPropertyOld<Real>("th_disp_neg")),
    _hc(declareProperty<Real>("hc")),
    _hc_old(getMaterialPropertyOld<Real>("hc")),
    _hc_mod(declareProperty<Real>("hc_mod")),
    _hc_mod_old(getMaterialPropertyOld<Real>("hc_mod")),
    // _hc_pos(declareProperty<Real>("hc_pos")),
    // _hc_neg(declareProperty<Real>("hc_neg")),
    _max_d(declareProperty<Real>("max_d")),
    _max_d_old(getMaterialPropertyOld<Real>("max_d")),
    _min_d(declareProperty<Real>("min_d")),
    _min_d_old(getMaterialPropertyOld<Real>("min_d")),
    _max_f(declareProperty<Real>("max_f")),
    _max_f_old(getMaterialPropertyOld<Real>("max_f")),
    _min_f(declareProperty<Real>("min_f")),
    _min_f_old(getMaterialPropertyOld<Real>("min_f")),
    _prev_pos_d(declareProperty<Real>("prev_pos_d")),
    _prev_pos_d_old(getMaterialPropertyOld<Real>("prev_pos_d")),
    _prev_neg_d(declareProperty<Real>("prev_neg_d")),
    _prev_neg_d_old(getMaterialPropertyOld<Real>("prev_neg_d")),
    _prev_pos_f(declareProperty<Real>("prev_pos_f")),
    _prev_pos_f_old(getMaterialPropertyOld<Real>("prev_pos_f")),
    _prev_neg_f(declareProperty<Real>("prev_neg_f")),
    _prev_neg_f_old(getMaterialPropertyOld<Real>("prev_neg_f")),
    _isdet_start(declareProperty<bool>("isdet_start")),
    _isdet_start_old(getMaterialPropertyOld<bool>("isdet_start")),
    _isplastic(declareProperty<bool>("isplastic")),
    _isplastic_old(getMaterialPropertyOld<bool>("isplastic")),
    _length_org(declareProperty<Real>("length_org")),
    _length_new(declareProperty<Real>("length_new"))
{
  // Checking for consistency between length of the provided displacements and rotations vector
  if (_ndisp != _nrot)
    mooseError("FuelCompressionSpring: The number of variables supplied in 'displacements' "
               "and 'rotations' input parameters must be equal.");

  // Fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    MooseVariable * disp_variable = getVar("displacements", i);
    _disp_num[i] = disp_variable->number(); // Displacement variable numbers in MOOSE

    MooseVariable * rot_variable = getVar("rotations", i);
    _rot_num[i] = rot_variable->number(); // Rotation variable numbers in MOOSE
  }

  /// Make sure that only hardening_constant or alpha_s and deterioration_constant or alpha_c
  /// is only used as a input by the user

  if (!parameters.isParamSetByUser("hardening_constant") && !isParamValid("alpha_s"))
   mooseError("FuelCompressionSpring: Either hardening_constant or alpha_s must be defined");

  if (parameters.isParamSetByUser("hardening_constant") && isParamValid("alpha_s"))
   mooseError("FuelCompressionSpring: Only the hardening_constant or only the alpha_s can be "
              "defined but not both");



///If alpha_s and alpha_c is provided calculate hardening_constant and deterioration_constant
   if (isParamValid("alpha_s") && !parameters.isParamSetByUser("hardening_constant"))
  {
    _alpha_s = getParam<Real>("alpha_s");
    _hardening_constant = ( _k_trans[0] * _alpha_s ) / (1.0 - _alpha_s);
  }



}

void
FuelCompressionSpring::initQpStatefulProperties()
{

  _deformations[_qp].zero();
  _spring_forces_global[_qp].zero();
  _spring_forces_test[_qp].zero();
  _elastic_deformation[_qp] = 0.0;
  _net_def[_qp]= 0.0;
  _yf[_qp]= _yield_force;
  _yf_neg[_qp] = _yield_force;
  _yf_pos[_qp] = _tension_factor * _yield_force ; // most likely unnecessary
  _th[_qp] = _th_disp;
  _th_disp_pos[_qp] = 1e16;
  _th_disp_neg[_qp] = _th_disp;
  _hc[_qp] = _hardening_constant;
  _hc_mod[_qp] = _hardening_constant;
  _max_d[_qp] = _tension_factor * _yield_force / _k_trans[_qp];  // most likely unnecessary
  _min_d[_qp] = _yield_force / _k_trans[_qp];
  _max_f[_qp] = _yield_force * _tension_factor;         //most likely unnecessary
  _min_f[_qp] = _yield_force;
  _k[_qp] = _k_trans[_qp];
  _isdet_start[_qp] = false;
  _isplastic[_qp] = false;
  _prev_pos_d[_qp] = 0;   // most likely unnecessary
  _prev_neg_d[_qp] = 0;
  _prev_pos_f[_qp] = _tension_factor * _yield_force ; // most likely unnecessary
  _prev_neg_f[_qp] = _yield_force;
  // std::cout<<" \n k0 i " << _k[_qp];
    // std::cout<<" \n k_trans i " << _k_trans[_qp];
    // std::cout<<" \n rf i " << _recovery_factor;

}

void
FuelCompressionSpring::computeQpProperties()
{
  // std::cout<<" \n k0 i " << _k[_qp];
  // Compute initial orientation and length of the spring in global coordinate system
  // Fetch the two nodes of the link element
  std::vector<const Node *> node;
  for (unsigned int i = 0; i < 2; ++i)
    node.push_back(_current_elem->node_ptr(i));
  RealGradient x_orientation;
  for (unsigned int i = 0; i < _ndisp; ++i)
    x_orientation(i) = (*node[1])(i) - (*node[0])(i);
  _length_org[_qp] = x_orientation.norm();
  x_orientation /= x_orientation.norm(); // Normalizing with length to get orientation

    // std::cout << " x = " << x_orientation << "\n";

  // Get y orientation of the spring in global coordinate system
  RealGradient y_orientation = getParam<RealGradient>("y_orientation");
  Real dot = x_orientation * y_orientation;
  // std::cout << " y = " << y_orientation << "\n";

  // Check if x and y orientations are perpendicular
  if (abs(dot) > 1e-4)
    mooseError("Error in FuelCompressionSpring: y_orientation should be perpendicular to "
               "the axis of the beam.");

  // Calculate z orientation in the global coordinate system as a cross product of the x and y
  // orientations
  RealGradient z_orientation = x_orientation.cross(y_orientation);
    // std::cout << " z = " << z_orientation << "\n";
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
FuelCompressionSpring::computeTotalRotation()
{
  _qp = 0;

  // Currently this forumation is limited to small deformations in the spring,
  // namely, it is assumed that there are no rigid body rotations in the spring,
  // and that the total rotation matrix from global to local coordinates
  // (calculated below) remains the same as the one at t = 0 throughout the
  // duration of the simulation.
  _total_global_to_local_rotation[_qp] = _original_global_to_local_rotation;


   // std::cout << " rotation matrix = " << _total_global_to_local_rotation[0] << "\n";
}

void
FuelCompressionSpring::computeDeformations()
{
  // std::cout<<" \n k0 i " << _k[_qp];
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
    // if (_d == 1)
    // {
      _global_disp0(i) = sol(node[0]->dof_number(nonlinear_sys.number(), _disp_num[i], 0));
      _global_disp1(i) = sol(node[1]->dof_number(nonlinear_sys.number(), _disp_num[i], 0));
      _global_rot0(i) = sol(node[0]->dof_number(nonlinear_sys.number(), _rot_num[i], 0));
      _global_rot1(i) = sol(node[1]->dof_number(nonlinear_sys.number(), _rot_num[i], 0));
    // }
    // else
    // {
    //   _global_disp0(i) = sol(node[1]->dof_number(nonlinear_sys.number(), _disp_num[i], 0));
    //   _global_disp1(i) = sol(node[0]->dof_number(nonlinear_sys.number(), _disp_num[i], 0));
    //   _global_rot0(i) = sol(node[1]->dof_number(nonlinear_sys.number(), _rot_num[i], 0));
    //   _global_rot1(i) = sol(node[0]->dof_number(nonlinear_sys.number(), _rot_num[i], 0));
    // }




  }
  // std::cout << "disp1 = " << _global_disp1<< "\n";
  // std::cout << "disp0 = " << _global_disp0<< "\n";
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
  _deformations[_qp] =  (_local_disp1 - _local_disp0);
  _rotations[_qp] = _local_rot1 - _local_rot0;
  // std::cout << " deformation " << _deformations[0] << "\n";

  RealGradient x_orientation_nw;
  for (unsigned int i = 0; i < _ndisp; ++i)
    x_orientation_nw(i) = ((*node[1])(i) + _global_disp1(i)) - ((*node[0])(i) + _global_disp0(i));
  _length_new[_qp] = x_orientation_nw.norm();
  // std::cout << "x_orientation_nw = " << x_orientation_nw << "\n";
  // std::cout << "final length = " << _length_new[0] << "\n";

  // _deformations[_qp](0) = _length_org[_qp] - _length_new[_qp];
}

void
FuelCompressionSpring::computeForces()
// spring forces = Kdd * deformations
// spring moments = Krr * rotations
{
if( _component <3)
{
  Real i = 0;
  Real disp_miss;
  bool isdispmiss = false;
  bool iselastic = false;
  Real def = MathUtils::round(_deformations[_qp](i) * 100000000)/100000000; //To counteract error in linear interpolation
  Real def_old = MathUtils::round(_deformations_old[_qp](i) * 100000000)/100000000;
// std::cout<<" \n spgf old " << _spring_forces_local_old[_qp](i);
  // //test if the spring is fully compressed and goes in other direction
  // if(_isdet_start[_qp] == false && def < -1 * _th_disp)
  // _isdet_start[_qp] = true;



  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_deformation[_qp] = MathUtils::round(_plastic_deformation_old[_qp] * 1e8) / 1e8;
//Added now
  _deterioration_variable[_qp] = _deterioration_variable_old[_qp];
  _k[_qp] = _k_old[_qp];
  _yf[_qp] = _yf_old[_qp];
  _yf_pos[_qp] = _yf_pos_old[_qp];
  _yf_neg[_qp] = _yf_neg_old[_qp];
  _th[_qp] = _th_old[_qp];
  _th_disp_pos[_qp] = _th_disp_pos_old[_qp];
  _th_disp_neg[_qp] = _th_disp_neg_old[_qp];
  _hc[_qp] = _hc_old[_qp];
  _hc_mod[_qp] = _hc_mod_old[_qp];
  _max_d[_qp] = _max_d_old[_qp];
  _min_d[_qp] = _min_d_old[_qp];
  _max_f[_qp] = _max_f_old[_qp];
  _min_f[_qp] = _min_f_old[_qp];
  _prev_pos_d[_qp] = _prev_pos_d_old[_qp];
  _prev_neg_d[_qp] = _prev_neg_d_old[_qp];
  _prev_pos_f[_qp] = _prev_pos_f_old[_qp];
  _prev_neg_f[_qp] = _prev_neg_f_old[_qp];
  _isdet_start[_qp] = _isdet_start_old[_qp];
  _isplastic[_qp] = _isplastic_old[_qp];

 // std::cout<<" \n k start  " << _k_trans[_qp];
 // std::cout<<" \n rf  " << _recovery_factor;
  // std::cout<<" \n is det changed  " << _isdet_start[_qp];
  // std::cout<<" \n k start changed " << _k[_qp];
  //
  // std::cout<<" \n plastic def start " << _plastic_deformation[_qp];
  // std::cout<<" \n prev_neg_d " << _prev_neg_d[_qp];
  // std::cout<<" \n prev neg f " << _prev_neg_f[_qp];
  // std::cout<<" \n min f " << _min_f[_qp];
  // std::cout<<" \n yf neg " << _yf_neg[_qp];

//Change of loading in negative region
// _k[_qp] = _k_trans[_qp];

if( def > def_old && std::abs(def) < std::abs( def_old) && _deformations_older[_qp](i) >= def_old && std::abs(_deformations_older[_qp](i)) <= std::abs(def_old))
{
  if( _isdet_start[_qp] == false)
  {
    _k[_qp] = _k_trans[_qp];
 // std::cout<<" \n k00  " << _k[_qp];
    if(std::abs(def_old) < _th_disp)
    {
      _prev_neg_d[_qp] = std::abs(def_old);
      _prev_neg_f[_qp] = std::abs(_spring_forces_local_old[_qp](i));
      if(std::abs(def_old) > _min_d[_qp])
      {
        _min_d[_qp] = std::abs(def_old);
        _min_f[_qp] = std::abs(_spring_forces_local_old[_qp](i));
        _yf_neg[_qp] = _min_f[_qp];
      }
    }
    else
    {
      _prev_neg_d[_qp] = _th_disp;
      _prev_neg_f[_qp] = _yield_force + (((_th_disp - ( _yield_force / _k_trans[_qp])) * _k_trans[_qp] * _hardening_constant ) / (_k_trans[_qp] + _hardening_constant));
      // if(std::abs(def_old) > _min_d[_qp])
      {
        _min_d[_qp] = _th_disp;
        _min_f[_qp] = _yield_force + (((_th_disp - ( _yield_force / _k_trans[_qp])) * _k_trans[_qp] * _hardening_constant ) / (_k_trans[_qp] + _hardening_constant));
        _yf_neg[_qp] = _min_f[_qp];
      }
    }
  }
  else
  {
    _k[_qp] = _k_sg;
  }


}

// std::cout<<" \n prev_neg_d 2 " << _prev_neg_d[_qp];
// std::cout<<" \n prev neg f 2 " << _prev_neg_f[_qp];
// std::cout<<" \n min f 2 " << _min_f[_qp];
// std::cout<<" \n yf neg 2 " << _yf_neg[_qp];
//
// std::cout<<" \n hardening variable 1 " << _hardening_variable[_qp];
//
//  std::cout<<" \n deformation " << def ;
//  std::cout<<" \n def old " << def_old;
// std::cout<<" \n def older " << _deformations_older[_qp](i);


Real deformation_increment = def - _deformations_old[_qp](i);
deformation_increment = MathUtils::round(deformation_increment * 100000000)/100000000;
 // std::cout<<" \n def increment " << deformation_increment;

if(_tension_factor != 0 && _min_c != 0)
{
  _net_def[_qp] =  MathUtils::round((1 - _recovery_factor) *  _min_c / _tension_factor / _k_trans[_qp] * 1e9) / 1e9;
}
else
{
  _net_def[_qp] =  MathUtils::round((1 - _recovery_factor) * _plastic_deformation[_qp] * 1e9) / 1e9;
}

// std::cout<<" \n net def  " << _net_def[_qp];
// std::cout<<" \n min c  " << _min_c;
//
 // std::cout<<" \n k0  " << _k[_qp];
 //  std::cout<<" \n k_trans  " << _k_trans[_qp];
//
// std::cout<<" \n det start beginning  " << _isdet_start[_qp];
//Changed GreaterEqual to Equal
if(MooseUtils::absoluteFuzzyLessThan(def,  _net_def[_qp],0) && MooseUtils::absoluteFuzzyEqual(def_old,  _net_def_old[_qp],1e-10) )
{
    if( _isdet_start[_qp] == false)
    {
      if((_Mahin_mod == true ) && MooseUtils::absoluteFuzzyGreaterThan( _min_f[_qp], _prev_neg_f[_qp] ,1e-8) && _isplastic[_qp] == true)
      {
        _k[_qp] = (_prev_neg_f[_qp] + (_spring_forces_local_old[_qp](i))) / (_prev_neg_d[_qp] + (def_old)) ;
         // std::cout<<" \n k01  " << _k[_qp];
        _yf[_qp] = _prev_neg_f[_qp];
        Real final_slope = ( _min_f[_qp] - _prev_neg_f[_qp]) / ( _min_d[_qp] - _prev_neg_d[_qp]);
        _hc_mod[_qp] = ( final_slope * _k_trans[_qp] ) /( _k_trans[_qp] - final_slope);
      }

      else
      {
        _k[_qp] = (_min_f[_qp] + (_spring_forces_local_old[_qp](i))) / (_min_d[_qp] + (def_old)) ;
        // std::cout<<" \n k02  " << _k[_qp];
        _yf[_qp] = _yf_neg[_qp];
        _hc_mod[_qp] = _hardening_constant;
      }
      _th[_qp] = _th_disp_neg[_qp];
      _hardening_variable[_qp] = 0;
    }
    else
    {
      _k[_qp] = _k_sg;
    }

}
// std::cout<<" \n prev_neg_d 3 " << _prev_neg_d[_qp];
// std::cout<<" \n prev neg f 3 " << _prev_neg_f[_qp];
// std::cout<<" \n min f 3 " << _min_f[_qp];
// std::cout<<" \n yf neg 3 " << _yf_neg[_qp];
// std::cout<<" \n yf 3 " << _yf[_qp];
//
// std::cout<<" \n hardening variable 2" << _hardening_variable[_qp];

//Modification for attaing slope in negative direction when it jumps the transition point
if(MooseUtils::absoluteFuzzyLessThan(def,  _net_def[_qp],0) && MooseUtils::absoluteFuzzyGreaterThan(def_old,  _net_def_old[_qp],1e-10) )
{
  if( _isdet_start[_qp] == false)
  {
    if((_Mahin_mod == true ) && MooseUtils::absoluteFuzzyGreaterThan( _min_f[_qp], _prev_neg_f[_qp] ,1e-8) && _isplastic[_qp] == true)
    {
      _k[_qp] = (_prev_neg_f[_qp] + (_spring_forces_local_old[_qp](i))) / (_prev_neg_d[_qp] + (_net_def_old[_qp])) ;
      // std::cout<<" \n k03  " << _k[_qp];
      _yf[_qp] = _prev_neg_f[_qp];
      Real final_slope = ( _min_f[_qp] - _prev_neg_f[_qp]) / ( _min_d[_qp] - _prev_neg_d[_qp]);
      _hc_mod[_qp] = ( final_slope * _k_trans[_qp] ) /( _k_trans[_qp] - final_slope);
    }

    else
    {
      _k[_qp] = (_min_f[_qp] + (_spring_forces_local_old[_qp](i))) / (_min_d[_qp] + (_net_def_old[_qp])) ;
      // std::cout<<" \n k04  " << _k[_qp];
      _yf[_qp] = _yf_neg[_qp];
      _hc_mod[_qp] = _hardening_constant;
    }
    _th[_qp] = _th_disp_neg[_qp];
    _hardening_variable[_qp] = 0;
  }
  else
  {
    _k[_qp] = _k_sg;
  }

}

// std::cout<<" \n prev_neg_d 4 " << _prev_neg_d[_qp];
// std::cout<<" \n prev neg f 4 " << _prev_neg_f[_qp];
// std::cout<<" \n yf neg 4 " << _yf_neg[_qp];
// std::cout<<" \n min f 4 " << _min_f[_qp];
// std::cout<<" \n yf 4 " << _yf[_qp];
// std::cout<<" \n hardening variable 3 " << _hardening_variable[_qp];

//Added for negative reversal  //think one more condition with regard to net deformation
if(def < 0  && def < def_old && std::abs(def) > std::abs( def_old) && _deformations_older[_qp](i) <= def_old && MooseUtils::absoluteFuzzyLessThan(def,  _net_def[_qp],0))
{
  if( _isdet_start[_qp] == false)
  {
    if((_Mahin_mod == true ) && MooseUtils::absoluteFuzzyGreaterThan( _min_f[_qp], _prev_neg_f[_qp] ,1e-8) && _isplastic[_qp] == true)
    {
      _k[_qp] = (_prev_neg_f[_qp] + (_spring_forces_local_old[_qp](i))) / (_prev_neg_d[_qp] + (_deformations_old[_qp](i))) ;
      // std::cout<<" \n k05  " << _k[_qp];
      _yf[_qp] = _prev_neg_f[_qp];
      Real final_slope = ( _min_f[_qp] - _prev_neg_f[_qp]) / ( _min_d[_qp] - _prev_neg_d[_qp]);
      _hc_mod[_qp] = ( final_slope * _k_trans[_qp] ) /( _k_trans[_qp] - final_slope);
    }

    else
    {
      _k[_qp] = (_min_f[_qp] + (_spring_forces_local_old[_qp](i))) / (_min_d[_qp] + (_deformations_old[_qp](i))) ;
      // std::cout<<" \n k06  " << _k[_qp];
      _yf[_qp] = _yf_neg[_qp];
      _hc_mod[_qp] = _hardening_constant;
    }
    _th[_qp] = _th_disp_neg[_qp];
    _hardening_variable[_qp] = 0;
  }
  else
  {
        _k[_qp] = _k_sg;
  }


}
//
// std::cout<<" \n hardening variable 4 " << _hardening_variable[_qp];
// std::cout<<" \n prev_neg_d 5 " << _prev_neg_d[_qp];
// std::cout<<" \n prev neg f 5 " << _prev_neg_f[_qp];
// std::cout<<" \n yf neg 5 " << _yf_neg[_qp];
// std::cout<<" \n min f 5 " << _min_f[_qp];
// std::cout<<" \n yf 5 " << _yf[_qp];



Real trial_force= 0;
 // std::cout<<" \n k  " << _k[_qp];
//   std::cout<<" \n sf old " << _spring_forces_local_old[_qp](i);
// std::cout<<" \n def inc " << deformation_increment;

trial_force = MathUtils::round((_spring_forces_local_old[_qp](i) + _k[_qp] * deformation_increment) * 1e9) / 1e9;
// std::cout<<" \n spgf " << _spring_forces_local_old[_qp](i);
// std::cout<<" \n def inc " << deformation_increment;
// std::cout<<" \n trial force0 " << trial_force;


//if(MooseUtils::absoluteFuzzyGreaterEqual(_length_new[_qp], _length_org[_qp] ,1e-10))
//if(MooseUtils::absoluteFuzzyGreaterThan(def, 0 ,1e-10))
if(MooseUtils::absoluteFuzzyGreaterEqual(def, 0 ,1e-10))
{
  _k[_qp] = _k_trans[_qp] * _tension_factor;
//  trial_force = (_spring_forces_local_old[_qp](i) + _k[_qp] * deformation_increment);
trial_force =  _k[_qp] * def;
 iselastic = true;
 // std::cout<<" \n tf " << trial_force;

}
 // std::cout<<" \n trial force " << trial_force;
 // std::cout<<" \n k2  " << _k[_qp];
if(MooseUtils::absoluteFuzzyLessEqual(_spring_forces_local_old[_qp](i), _k_trans[_qp] * def_old * _tension_factor ,1e-10) && MooseUtils::absoluteFuzzyGreaterThan(trial_force,  _k_trans[_qp] * _deformations[_qp](i) * _tension_factor ,1e-15) && MooseUtils::absoluteFuzzyGreaterThan(def, def_old , 1e-10) && def<0) //netdf changed to df
{
  if (_k[_qp] != ( _k_trans[_qp] * _tension_factor) )
  {
    if( _tension_factor !=0)
    {
      disp_miss = ((_spring_forces_local_old[_qp](i) - def_old * _k[_qp]  ) /  (  _tension_factor * _k_trans[_qp] - _k[_qp]  ));
    }
    else
    {
      disp_miss = _net_def[_qp];
    }
// std::cout<<" \n disp_miss  " << disp_miss;
    _min_c = disp_miss * _tension_factor * _k_trans[_qp];
    // std::cout<<" \n min cc  " << _min_c;
    isdispmiss = true ;
  }


  _k[_qp] = _k_trans[_qp] * _tension_factor;

  if(disp_miss != def_old )
  {
    trial_force = def * _k[_qp];
  }
  else
  {
    trial_force = (_spring_forces_local_old[_qp](i) + _k[_qp] * deformation_increment);
  }
}




// std::cout<<" \n hardening variable 5 " << _hardening_variable[_qp];
 // std::cout<<" \n k3  " << _k[_qp];
 // std::cout<<" \n trial force f" << trial_force;
Real elastic_deformation_increment = deformation_increment;
// std::cout << "\n is det start = " << _isdet_start[_qp];

if(_isdet_start[_qp] == false)
{
Real effective_force = trial_force - _hardening_variable[_qp];
// std::cout << "\n eff force = " << effective_force;
// std::cout << "\n hardening_variable 6 = " << _hardening_variable[_qp];
// std::cout << "\n yf = " << _yf[_qp];
Real yield_condn = std::abs(effective_force) - _yf[_qp];
Real iteration = 0;
Real plastic_deformation_increment = 0.0;

//  std::cout<<" \n yield condition " << yield_condn;
// std::cout<<" \n th disp " << _th_disp;
// //Real elastic_deformation_increment = deformation_increment;
// std::cout << "\n is det start 1 = " << _isdet_start[_qp];
 if (MooseUtils::absoluteFuzzyGreaterThan(yield_condn, 0.0, 1e-12)  && MooseUtils::absoluteFuzzyLessThan(def, 0, 1e-16) && (std::abs(def)<= _th_disp))//MooseUtils::absoluteFuzzyLessEqual(std::abs(def),  _th_disp, 1e-9))

{
  // std::cout << "\n is det start a = " << _isdet_start[_qp];
  if(_isplastic[_qp]==false)
  {
    _isplastic[_qp] = true;
  }
// std::cout << "\n is det start b = " << _isdet_start[_qp];
  _k[_qp] = _k_trans[_qp];
  trial_force = _spring_forces_local_old[_qp](i) + _k[_qp] * deformation_increment;
  effective_force = trial_force - _hardening_variable[_qp];

    _hc[_qp] = _hc_mod[_qp];

  Real residual = std::abs(trial_force  - _hardening_variable[_qp]) - _yf[_qp] -
                  (_k[_qp] + _hc[_qp]) * plastic_deformation_increment;

// std::cout << "\n is det start c = " << _isdet_start[_qp];
  while (std::abs(residual) > _absolute_tolerance)
  {

    Real scalar = (std::abs(trial_force - _hardening_variable[_qp]) - _yf[_qp] -
                   (_k[_qp] + _hc[_qp] )* plastic_deformation_increment) /
                  (_k[_qp] + _hc[_qp]);
    plastic_deformation_increment += scalar;

    residual = std::abs(trial_force  - _hardening_variable[_qp]) - _yf[_qp] -
                    (_k[_qp] + _hc[_qp]) * plastic_deformation_increment;


    ++iteration;
    if (iteration > _max_its) // not converging
      throw MooseException("FuelCompressionSpring: Plasticity model did not converge");
  }
  // std::cout << "\n is det start d = " << _isdet_start[_qp];
  plastic_deformation_increment *= MathUtils::sign(trial_force - _hardening_variable[_qp]);
  _plastic_deformation[_qp] += plastic_deformation_increment;
  elastic_deformation_increment = deformation_increment - plastic_deformation_increment;
  _hardening_variable[_qp] += _hc[_qp] * plastic_deformation_increment;
  _spring_forces_local[_qp](i) = _spring_forces_local_old[_qp](i) + _k[_qp] * elastic_deformation_increment;
// std::cout<<" \n spring forces local plastic " << _spring_forces_local[_qp](i);
// std::cout<<" \n hardening variable 7" << _hardening_variable[_qp];
// std::cout << "\n is det start 2 = " << _isdet_start[_qp];
}
else
{
_elastic_deformation[_qp] = def - _plastic_deformation[_qp];
if(MooseUtils::absoluteFuzzyGreaterEqual(def, -1 * _th_disp, 1e-10))
 {

if(MooseUtils::absoluteFuzzyGreaterThan(def, _plastic_deformation[_qp], 1e-10 ) && MooseUtils::absoluteFuzzyLessThan(def_old, _plastic_deformation_old [_qp],1e-10) && def<0)
   {
     _spring_forces_local[_qp](i) = _spring_forces_local_old[_qp](i) + _k[_qp] * (_plastic_deformation[_qp] - def_old) + _k_trans[_qp] * (def - _plastic_deformation[_qp]);
   // std::cout<<" \n spring forces local transition " << _spring_forces_local[_qp](i);
   }
   else if (MooseUtils::absoluteFuzzyLessThan(def,  _net_def[_qp],0) && MooseUtils::absoluteFuzzyGreaterThan(def_old,  _net_def_old[_qp],1e-10) )
//added else if block for the time skipping
   {
     _spring_forces_local[_qp](i) = _spring_forces_local_old[_qp](i) + _tension_factor * _k_trans[_qp] * (_net_def[_qp] - def_old) + _k[_qp] * (def - _net_def[_qp]);
   }
   // else if (isdispmiss == true)
   // {
   //  _spring_forces_local[_qp](i) = trial_force;
   // }
   else
   {
     _spring_forces_local[_qp](i) = _spring_forces_local_old[_qp](i) + _k[_qp] * elastic_deformation_increment;
    // std::cout<<" \n spring forces local0 " << _spring_forces_local[_qp](i);

   }

   // std::cout<<" \n hardening variable 8 " << _hardening_variable[_qp];

   if (iselastic == true || isdispmiss == true)
   {
       _spring_forces_local[_qp](i) = trial_force;
       // std::cout<<" \n spring forces local0 " << _spring_forces_local[_qp](i);
   }

}
else
{
  if(MooseUtils::absoluteFuzzyGreaterThan(def_old, -1 * _th_disp, 1e-10) && def <0)
  {
    Real elpl_def_inc = - _th_disp - def_old;
    Real el_def_inc = def + _th_disp;
    Real plastic_deformation_increment = 0.0;
//   std::cout<<" \n elpl_def_inc " << elpl_def_inc;
// std::cout<<" \n el_def_inc " << el_def_inc;

    trial_force = _spring_forces_local_old[_qp](i) + _k_trans[_qp] * elpl_def_inc;
    effective_force = trial_force - _hardening_variable[_qp];

      _hc[_qp] = _hc_mod[_qp];

    Real residual = std::abs(trial_force  - _hardening_variable[_qp]) - _yf[_qp] -
                    (_k_trans[_qp] + _hc[_qp]) * plastic_deformation_increment;


    while (std::abs(residual) > _absolute_tolerance)
    {

      Real scalar = (std::abs(trial_force - _hardening_variable[_qp]) - _yf[_qp] -
                     (_k_trans[_qp] + _hc[_qp] )* plastic_deformation_increment) /
                    (_k_trans[_qp] + _hc[_qp]);
      plastic_deformation_increment += scalar;

      residual = std::abs(trial_force  - _hardening_variable[_qp]) - _yf[_qp] -
                      (_k_trans[_qp] + _hc[_qp]) * plastic_deformation_increment;


      ++iteration;
      if (iteration > _max_its) // not converging
        throw MooseException("FuelCompressionSpring: Plasticity model did not converge");
    }
    plastic_deformation_increment *= MathUtils::sign(trial_force - _hardening_variable[_qp]);
    _plastic_deformation[_qp] += plastic_deformation_increment;
    elastic_deformation_increment = elpl_def_inc - plastic_deformation_increment;
    _hardening_variable[_qp] += _hc[_qp] * plastic_deformation_increment;
//     std::cout<<" \n hardening variable 9 " << _hardening_variable[_qp];
// std::cout<<" \n elastic deformation increment " << elastic_deformation_increment;
// std::cout<<" \n spring forces old " << _spring_forces_local_old[_qp](i);
    _spring_forces_local[_qp](i) = _spring_forces_local_old[_qp](i) + _k_trans[_qp] * elastic_deformation_increment ;

// std::cout<<" \n spring forces local before add " << _spring_forces_local[_qp](i);
    _spring_forces_local[_qp](i) += _k_sg * el_def_inc;
// std::cout<<" \n spring forces local after add " << _spring_forces_local[_qp](i);
if(_complete_det == true && _isdet_start[_qp] == false)
{
  _recovery_factor = 0;
  _isdet_start[_qp] = true;
  _plastic_deformation[_qp] = -(_th_disp  - (_yield_force + (((_th_disp - ( _yield_force / _k_trans[_qp])) * _k_trans[_qp] * _hardening_constant ) / (_k_trans[_qp] + _hardening_constant))) / (_k_sg));
  _k[_qp] = _k_sg;
  // std::cout << "\n det starts 1" ;
  // std::cout << "\n is det start 4 = " << _isdet_start[_qp];
}

    //added till here
  }
  else
  {
    _spring_forces_local[_qp](i) = _spring_forces_local_old[_qp](i) + _k_sg * elastic_deformation_increment;
    // std::cout<<" \n sf sg " << _spring_forces_local[_qp](i);
    if(_complete_det == true && _isdet_start[_qp] == false)
    {
      _recovery_factor = 0;
      _isdet_start[_qp] = true;
      _plastic_deformation[_qp] = -(_th_disp - (_yield_force + (((_th_disp - ( _yield_force / _k_trans[_qp])) * _k_trans[_qp] * _hardening_constant ) / (_k_trans[_qp] + _hardening_constant))) / (_k_sg));
      _k[_qp] = _k_sg;
      // std::cout << "\n is det start 5 = " << _isdet_start[_qp];
      // std::cout << "\n det starts 2" ;
    }
  }

}
}
}
else
{
  if(MooseUtils::absoluteFuzzyGreaterThan(def, _plastic_deformation[_qp], 1e-10 ) && MooseUtils::absoluteFuzzyLessThan(def_old, _plastic_deformation_old [_qp],1e-10) && def<0)
     {
       _spring_forces_local[_qp](i) = _spring_forces_local_old[_qp](i) + _k[_qp] * (_plastic_deformation[_qp] - def_old) + _k_trans[_qp] * (def - _plastic_deformation[_qp]);
     // std::cout<<" \n spring forces local transition " << _spring_forces_local[_qp](i);
     }
     else if (MooseUtils::absoluteFuzzyLessThan(def,  _net_def[_qp],0) && MooseUtils::absoluteFuzzyGreaterThan(def_old,  _net_def_old[_qp],1e-10) )
  //added else if block for the time skipping
     {
       _spring_forces_local[_qp](i) = _spring_forces_local_old[_qp](i) + _tension_factor * _k_trans[_qp] * (_net_def[_qp] - def_old) + _k[_qp] * (def - _net_def[_qp]);
     }
     // else if ( isdispmiss == true)
     // {
     //   _spring_forces_local[_qp](i) = trial_force;
     // }
     else
     {
       _spring_forces_local[_qp](i) = _spring_forces_local_old[_qp](i) + _k[_qp] * elastic_deformation_increment;
     // std::cout<<" \n spring forces local0 " << _spring_forces_local[_qp](i);

     }

     if (iselastic == true || isdispmiss == true)
     {
         _spring_forces_local[_qp](i) = trial_force;
          // std::cout<<" \n spring force =  " << _spring_forces_local(0);
     }
}

}
// std::cout<<" \n plastic def " << _plastic_deformation[_qp];
// std::cout<<" \n hardening variable end " << _hardening_variable[_qp];
//
//
// std::cout<<" \n kf  " << _k[_qp];
// std::cout<<" \n det start end  " << _isdet_start[_qp];


for(unsigned int i=0; i<3; i++)
{

  if(i != 0)
  {
    _spring_forces_local[_qp](i) = 1e-6 * _deformations[_qp](i);
  }

}
// std::cout << "spring forces local " << _spring_forces_local[_qp];

  // _spring_forces_local(1) = 0;
  // _spring_forces_local(2) = 0;
 // std::cout<<" \n spring forces local=  " << _spring_forces_local;
// std::cout<<" \n tr " << _total_global_to_local_rotation[_qp].transpose();

  // convert local forces to global
  _spring_forces_global[_qp] =
      _total_global_to_local_rotation[_qp].transpose() * _spring_forces_local[_qp];

    // std::cout<<" \n spring forces global= " << _spring_forces_global[_qp];


for(unsigned int i=3; i<6; i++)
{
  // if(i=!_component)              //Uncomment for rotational
    _spring_moments_local(i-3) = 1e-6 * _rotations[_qp](i-3);
  }


  // convert local moments to global
  _spring_moments_global[_qp] =
      _total_global_to_local_rotation[_qp].transpose() * _spring_moments_local;


}





void
FuelCompressionSpring::computeStiffnessMatrix()
{
  // std::cout<<" \n k0 i " << _k[_qp];
  // The stiffness matrix is of the form
  // |  kdd  kdr  |
  // |  krd  krr  |
  // where kdd, krr, krd and kdr are all RankTwoTensors (3x3 matrix) and
  // matrix symmetry is assumed, namely, krd = kdr.transpose()
  // This implementation of the spring element has only diagonal stiffness
  // terms. Therefore, Kdr and Krd are zero.

  // calculating deformation stiffnesses
  for(unsigned int i=0; i<3; i++)
  {
    if(i==0)
    {
      _kdd_local(i, i) = _k[_qp];
    }

    else
    {
      _kdd_local(i, i) = 1e-6; ///Actually 0
    }

  }


  // convert stiffness matrix from local to global
  _kdd[_qp] = _total_global_to_local_rotation[_qp].transpose() * _kdd_local *
              _total_global_to_local_rotation[_qp];

  for(unsigned int i=3; i<6; i++)
  {
    // if(i==_component)
    //   _krr_local(i-3, i-3) = _k_trans[_qp];  //Uncomment for rotational
    // else
      _krr_local(i-3, i-3) = 1e-6;
  }

  // convert stiffness matrix from local to global
  _krr[_qp] = _total_global_to_local_rotation[_qp].transpose() * _krr_local *
              _total_global_to_local_rotation[_qp];

               // std::cout<<" \n kdd " << _kdd[_qp];
}
