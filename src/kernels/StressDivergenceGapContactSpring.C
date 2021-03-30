/*************************************************/
/*           DO NOT MODIFY THIS HEADER           */
/*                                               */
/*                     MASTODON                  */
/*                                               */
/*    (c) 2015 Battelle Energy Alliance, LLC     */
/*            ALL RIGHTS RESERVED                */
/*                                               */
/*   Prepared by Battelle Energy Alliance, LLC   */
/*     With the U. S. Department of Energy       */
/*                                               */
/*     See COPYRIGHT for full restrictions       */
/*************************************************/

// MASTODON includes
#include "StressDivergenceGapContactSpring.h"

// MOOSE includes
#include "Assembly.h"
#include "Material.h"
#include "MooseVariable.h"
#include "SystemBase.h"
#include "NonlinearSystem.h"
#include "MooseMesh.h"

// libmesh includes
#include "libmesh/quadrature.h"

registerMooseObject("MastodonApp", StressDivergenceGapContactSpring);


InputParameters
StressDivergenceGapContactSpring::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Kernel for gap contact spring element");
  params.addRequiredParam<unsigned int>(
      "component",
      "An integer corresponding to the direction "
      "the variable this kernel acts in. (0 for x, "
      "1 for y, 2 for z, 3 for rot_x, 4 for rot_y and 5 for rot_z).");
  params.addRequiredCoupledVar("displacements", "The displacement variables for gap contact spring.");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

StressDivergenceGapContactSpring::StressDivergenceGapContactSpring(const InputParameters & parameters)
  : Kernel(parameters),
    _component(getParam<unsigned int>("component")),
    _ndisp(coupledComponents("displacements")),
    _disp_var(_ndisp),
    _Fg(getMaterialPropertyByName<ColumnMajorMatrix>("contact_forces")),
    _Kg(getMaterialPropertyByName<ColumnMajorMatrix>("global_stiffness_matrix"))
{
  if (_component > 3)
    mooseError("Error in StressDivergenceGapContactSpring block ",
               name(),
               ". Please enter an integer value between 0 and 2 for the 'component' parameter.");


  for (unsigned int i = 0; i < _ndisp; ++i)
    _disp_var[i] = coupled("displacements", i);

}

void
StressDivergenceGapContactSpring::computeResidual()
{
  // Accessing residual vector, re, from MOOSE assembly
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  mooseAssert(re.size() == 2, "Spring element has two nodes only.");
  _local_re.resize(re.size());
  _local_re.zero();

  // Calculating residual for node 0 (external forces on node 0)
  _local_re(0) = _Fg[0](_component);

  // Calculating residual for node 1 (external forces on node 1)
  _local_re(1) = _Fg[0](_component + 3);

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _save_in.size(); ++i)
      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
  }
}

void
StressDivergenceGapContactSpring::computeJacobian()
{
  // Access Jacobian; size is n x n (n is number of nodes)
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  // i and j are looping over nodes
  for (unsigned int i = 0; i < _test.size(); ++i)
    for (unsigned int j = 0; j < _phi.size(); ++j)
      _local_ke(i, j) += _Kg[0](i * 3 + _component, j * 3 + _component);

  ke += _local_ke;

  if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Number> diag(rows);
    for (unsigned int i = 0; i < rows; ++i)
      diag(i) = _local_ke(i, i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
      _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
  }
}

void
StressDivergenceGapContactSpring::computeOffDiagJacobian(const unsigned int jvar_num)
// coupling one variable to another (disp x to disp y, etc)
{
  if (jvar_num == _var.number())
    // jacobian calculation if jvar is the same as the current variable i.e.,
    // diagonal elements
    computeJacobian();

  else
  // jacobian calculation for off-diagonal elements
  {
    unsigned int coupled_component = 0;
    bool coupled = false;
    // finding which variable jvar is
    for (unsigned int i = 0; i < _ndisp; ++i)
    {
      if (jvar_num == _disp_var[i])
      {
        coupled_component = i;
        coupled = true;
        break;
      }
      // else if (jvar_num == _rot_var[i])
      // {
      //   coupled_component = i + 3;
      //   coupled = true;
      //   break;
      // }
    }
    // getting the jacobian from assembly
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), jvar_num);
    if (coupled)
    {
      for (unsigned int i = 0; i < _test.size(); ++i)
        for (unsigned int j = 0; j < _phi.size(); ++j)
          ke(i, j) += _Kg[0](i * 3 + _component, j * 3 + coupled_component);
    }
  }
}
