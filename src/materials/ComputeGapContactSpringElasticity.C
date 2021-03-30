// MASTODON includes
#include "ComputeGapContactSpringElasticity.h"

// MOOSE includes
#include "MooseMesh.h"
#include "Assembly.h"
#include "NonlinearSystem.h"
#include "MooseVariable.h"
#include "ColumnMajorMatrix.h"
#include "MathUtils.h"
#include "MooseUtils.h"

// libmesh includes
#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

registerMooseObject("MastodonApp", ComputeGapContactSpringElasticity);

template <>
InputParameters
validParams<ComputeGapContactSpringElasticity>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Node to node contact element modeled to represent the normal"
                             "contact forces and tangential Coulomb frictional forces after the"
                             "gap between the nodes is closed");
  params.addRequiredCoupledVar("displacements",
                               "The displacement variables appropriate for the simulation of "
                               "geometry and coordinate system.");
  params.addRequiredParam<RealGradient>("y_orientation",
                                        "Orientation of the y direction along with tangential "
                                        "contact stiffness is provided. This should be "
                                        "perpendicular to the axis of the spring.");

  // Input parameters
  // params.addRequiredParam<Real>("gap", "Gap between the nodes");
  params.addRequiredParam<Real>(
      "Kn", "Contact normal stiffness");
  params.addRequiredParam<Real>(
      "Kt", "Tangential contact stiffness");
  params.addParam<Real>(
      "mu", 0.0, "Coefficient of friction.");
  // params.addParam<Real>("rel_tol", 1e-6, "Relative tolerance for error in adaptive algorithm");
  // params.addParam<Real>("abs_tol", 1e-6, "Absolute tolerance for error in adaptive algorithm");
  params.set<MooseEnum>("constant_on") = "ELEMENT";
  return params;
}

ComputeGapContactSpringElasticity::ComputeGapContactSpringElasticity(const InputParameters & parameters)
  : Material(parameters),
    _ndisp(coupledComponents("displacements")),
    _disp_num(3),
    _length(declareProperty<Real>("initial_gap")),
    _Kn(getParam<Real>("Kn")),
    _Kt(getParam<Real>("Kt")),
    _mu(getParam<Real>("mu")),
    // _rel_tol(getParam<Real>("rel_tol")),
    // _abs_tol(getParam<Real>("abs_tol")),
    _local_disp(declareProperty<ColumnMajorMatrix>("local_deformations")),
    _basic_def(declareProperty<ColumnMajorMatrix>("basic_deformations")),
    _basic_def_old(getMaterialPropertyOld<ColumnMajorMatrix>("basic_deformations")),
    _Fb(declareProperty<ColumnMajorMatrix>("basic_force")),
    _Fb_old(getMaterialPropertyOld<ColumnMajorMatrix>("basic_force")),
    _Fl(declareProperty<ColumnMajorMatrix>("local_forces")),
    _Fg(declareProperty<ColumnMajorMatrix>("contact_forces")),
    _Kb(declareProperty<ColumnMajorMatrix>("basic_stiffness")),
    _Kl(declareProperty<ColumnMajorMatrix>("local_stiffness_matrix")),
    _Kg(declareProperty<ColumnMajorMatrix>("global_stiffness_matrix")),
    _transform_gl(declareProperty<ColumnMajorMatrix>("total_global_to_local_transformation")),
    _transform_lb(declareProperty<ColumnMajorMatrix>("total_local_to_basic_transformation"))

{
  // Fetch coupled variables (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    MooseVariable * disp_variable = getVar("displacements", i);
    _disp_num[i] = disp_variable->number(); // Displacement variable numbers in MOOSE
  }
}

void
ComputeGapContactSpringElasticity::initQpStatefulProperties()
{
}

void
ComputeGapContactSpringElasticity::computeQpProperties()
{
  // Compute transformation matrices
  computeTransformationMatrix();

  // Compute the global displacement values
  computeDeformation();

  // Initialize the stiffness matrices and force vectors
  initializeContactSpring();

  // Compute axial forces and stiffness terms
  computeForces();

  // Finalize forces and stiffness matrix and convert them into global co-ordinate system
  finalize();
}

void
ComputeGapContactSpringElasticity::computeTransformationMatrix()
{
  // Fetch the two nodes of the link element
  std::vector<const Node *> node;
  for (unsigned int i = 0; i < 2; ++i)
    node.push_back(_current_elem->node_ptr(i));

  // Defining orientation of contact spring (direction cosines)
  RealGradient x_orientation;
  for (unsigned int i = 0; i < _ndisp; ++i)
    x_orientation(i) = (*node[1])(i) - (*node[0])(i);

  // Length of the contact spring element
  //Also, the initial gap between the nodes
  _length[_qp] = x_orientation.norm();

  if (_length[_qp] == 0.0)
    mooseError("Error in ComputeGapContactSpringElasticity material block, ",
               name(),
               ". contact spring element cannot be of zero length.");

  // Normalizing the orientation vector
  x_orientation /= _length[_qp];

  // Get y orientation of the contact spring in global coordinate system
  RealGradient y_orientation = getParam<RealGradient>("y_orientation");
  Real dot = x_orientation * y_orientation;

  // Check if x and y orientations are perpendicular
  if (abs(dot) > 1e-4)
    mooseError("Error in ComputeGapContactSpringElasticity material block, ",
               name(),
               ". y_orientation should be perpendicular to the axis of the contact spring.");

  //Calculate z orientation of the contact spring as cross product of the
  // x and y orientations

  RealGradient z_orientation = x_orientation.cross(y_orientation);

  // Transformation matrix from global to local coordinate system
  _transform_gl[_qp].reshape(6, 6);
  _transform_gl[_qp].zero();
  _transform_gl[_qp](0, 0) = _transform_gl[_qp](3, 3) = x_orientation(0);
  _transform_gl[_qp](0, 1) = _transform_gl[_qp](3, 4) = x_orientation(1);
  _transform_gl[_qp](0, 2) = _transform_gl[_qp](3, 5) = x_orientation(2);
  _transform_gl[_qp](1, 0) = _transform_gl[_qp](4, 3) = y_orientation(0);
  _transform_gl[_qp](1, 1) = _transform_gl[_qp](4, 4) = y_orientation(1);
  _transform_gl[_qp](1, 2) = _transform_gl[_qp](4, 5) = y_orientation(2);
  _transform_gl[_qp](2, 0) = _transform_gl[_qp](5, 3) = z_orientation(0);
  _transform_gl[_qp](2, 1) = _transform_gl[_qp](5, 4) = z_orientation(1);
  _transform_gl[_qp](2, 2) = _transform_gl[_qp](5, 5) = z_orientation(2);




  // Transformation matrix from local to basic system
  _transform_lb[_qp].reshape(3, 6);
  _transform_lb[_qp].zero();
  _transform_lb[_qp](0, 0) = _transform_lb[_qp](1, 1) = _transform_lb[_qp](2, 2) = -1;
  _transform_lb[_qp](0, 3) = _transform_lb[_qp](1, 4) = _transform_lb[_qp](2, 5) = 1;


}

void
ComputeGapContactSpringElasticity::computeDeformation()
{
  // Fetch the two end nodes for _current_elem
  std::vector<const Node *> node;
  for (unsigned int i = 0; i < 2; ++i)
    node.push_back(_current_elem->node_ptr(i));

  // Fetch the solution for the two end nodes at current time
  NonlinearSystemBase & nonlinear_sys = _fe_problem.getNonlinearSystemBase();
  const NumericVector<Number> & sol = *nonlinear_sys.currentSolution();
  const NumericVector<Number> & sol_old = nonlinear_sys.solutionOld();

  // Calculating global displacements 6 x 1 matrix with
  // first three rows corresponding to node 0 dofs and next three to node 1 dofs
  ColumnMajorMatrix global_disp(6, 1);
  ColumnMajorMatrix global_disp_old(6, 1);

  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    global_disp(i) =
        sol(node[0]->dof_number(nonlinear_sys.number(), _disp_num[i], 0)); // node 0 displacements
    global_disp(i + 3) =
        sol(node[1]->dof_number(nonlinear_sys.number(), _disp_num[i], 0)); // node 1 displacements
    global_disp_old(i) = sol_old(
        node[0]->dof_number(nonlinear_sys.number(), _disp_num[i], 0)); // node 0 displacements
    global_disp_old(i + 3) = sol_old(
        node[1]->dof_number(nonlinear_sys.number(), _disp_num[i], 0)); // node 1 displacements
  }

  //Convert the global deformations to basic basic deformations
  //Resize the local and basic deformation matrix
  _local_disp[_qp].reshape(6, 1);
  _basic_def[_qp].reshape(3, 1);
  // _basic_def_old[_qp].reshape(3, 1);

  _local_disp[_qp] = _transform_gl[_qp] * global_disp;
  _basic_def[_qp] = _transform_lb[_qp] * _local_disp[_qp];


}

void
ComputeGapContactSpringElasticity::initializeContactSpring()
{
  // Initialize stiffness matrices
  _Kb[_qp].reshape(3, 3);
  if(MooseUtils::absoluteFuzzyLessThan( _basic_def[_qp](0, 0), 0 , 1e-10) &&
   MooseUtils::absoluteFuzzyGreaterEqual(std::abs( _basic_def[_qp](0, 0)), _length[_qp], 1e-10))
   {
     _Kb[_qp](0, 0) = _Kn;     //normal stiffness
     _Kb[_qp](1, 1) = 0;     //tangential stiffness
     _Kb[_qp](2, 2) = 0;     //tangential stiffness
   }
   else
   {
     _Kb[_qp](0, 0) = 0;     //normal stiffness
     _Kb[_qp](1, 1) = 0;     //tangential stiffness
     _Kb[_qp](2, 2) = 0;     //tangential stiffness
   }


    // std::cout<< "basic stiffness " << _Kb << "\n";

  // Initialize force matrices
  _Fb[_qp].reshape(3, 1);   //force matrix in basic system
  _Fb[_qp].zero();
  _Fl[_qp].reshape(6, 1);   //forces in local system
  _Fl[_qp].zero();
  _Fg[_qp].reshape(6, 1);   //forces in global system
  _Fg[_qp].zero();
}



void
ComputeGapContactSpringElasticity::computeForces()
{
//  Real gap = basic_def[_qp](0, 0) - _length[_qp];
   // std::cout<< "normal def " << _basic_def[_qp](0, 0) << "\n";
   // std::cout<< "initial gap " << _length[_qp] << "\n";
  if(MooseUtils::absoluteFuzzyLessThan( _basic_def[_qp](0, 0), 0 , 1e-10) &&
   MooseUtils::absoluteFuzzyGreaterEqual(std::abs( _basic_def[_qp](0, 0)), _length[_qp], 1e-10))
  {
    _Fb[_qp](0, 0) =   _Kn * (_basic_def[_qp](0, 0) +  _length[_qp] );


    Real frictional_force = std::abs(_mu * _Fb[_qp](0, 0));
    Real old_res_tan_force = sqrt(pow(_Fb_old[_qp](1, 0),2) + pow(_Fb_old[_qp](2, 0),2));
    if(MooseUtils::absoluteFuzzyGreaterEqual(old_res_tan_force, frictional_force, 1e-10))
    {
      if(MooseUtils::absoluteFuzzyEqual(_basic_def[_qp](1,0), 0.0, 1e-10 ))
      _Fb[_qp](1, 0) =  0;
      else
      _Fb[_qp](1, 0) =  _mu * std::abs(_Fb[_qp](0, 0));
      if(MooseUtils::absoluteFuzzyEqual(_basic_def[_qp](2,0), 0.0, 1e-10 ))
      _Fb[_qp](2, 0) =  0;
      else
      _Fb[_qp](2, 0) = _mu * std::abs(_Fb[_qp](0, 0));
    }
    else
    {
       _Fb[_qp](1, 0) = _Kt * _basic_def[_qp](1, 0);
       _Fb[_qp](2, 0) = _Kt * _basic_def[_qp](2, 0);

      Real res_tan_force = sqrt(pow(_Fb[_qp](1, 0),2) + pow(_Fb[_qp](2, 0),2));

      if(MooseUtils::absoluteFuzzyGreaterEqual( res_tan_force, frictional_force, 1e-10))
      {
        if(MooseUtils::absoluteFuzzyEqual(_basic_def[_qp](1,0), 0.0, 1e-10 ))
        _Fb[_qp](1, 0) =  0;
        else
        _Fb[_qp](1, 0) =  _mu * std::abs(_Fb[_qp](0, 0));
        if(MooseUtils::absoluteFuzzyEqual(_basic_def[_qp](2,0), 0.0, 1e-10 ))
        _Fb[_qp](2, 0) =  0;
        else
        _Fb[_qp](2, 0) = _mu * std::abs(_Fb[_qp](0, 0));
      }

    }
  }
  else
  {
    _Fb[_qp](0, 0) = 0;
   _Fb[_qp](1, 0) =  0;
   _Fb[_qp](2, 0) =  0;
  }
}

void
ComputeGapContactSpringElasticity::finalize()
{

  // Convert forces from basic to local to global coordinate system
  _Fl[_qp] = _transform_lb[_qp].transpose() * _Fb[_qp]; // local forces
  _Fg[_qp] = _transform_gl[_qp].transpose() * _Fl[_qp]; // global forces



  // Convert stiffness matrix from basic to local coordinate system
  _Kl[_qp] = _transform_lb[_qp].transpose() * _Kb[_qp] * _transform_lb[_qp];


  // Converting stiffness matrix from local to global coordinate system
  _Kg[_qp] = _transform_gl[_qp].transpose() * _Kl[_qp] * _transform_gl[_qp];


}
