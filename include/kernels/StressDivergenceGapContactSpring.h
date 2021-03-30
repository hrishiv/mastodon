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

#ifndef STRESSDIVERGENCEGAPCONTACTSPRING_H
#define STRESSDIVERGENCEGAPCONTACTSPRING_H

#include "Kernel.h"
#include "ColumnMajorMatrix.h"

// Forward Declarations
// class StressDivergenceGapContactSpring;
//
// template <>
// InputParameters validParams<StressDivergenceGapContactSpring>();

class StressDivergenceGapContactSpring : public Kernel
{
public:
  static InputParameters validParams();
  StressDivergenceGapContactSpring(const InputParameters & parameters);
  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(unsigned int jvar) override;
  // using Kernel::computeOffDiagJacobian;

protected:
  virtual Real computeQpResidual() override { return 0.0; }

  /// Direction along which force/moment is calculated
  const unsigned int _component;

  /// Number of coupled displacement variables
  unsigned int _ndisp;

  /// Variable numbers corresponding to displacement variables
  std::vector<unsigned int> _disp_var;

  /// Number of coupled rotational variables
  //unsigned int _nrot;

  /// Variable numbers corresponding to rotational variables
  //std::vector<unsigned int> _rot_var;

  /// Global isolator forces
  const MaterialProperty<ColumnMajorMatrix> & _Fg;

  /// GLobal stiffness matrix
  const MaterialProperty<ColumnMajorMatrix> & _Kg;
};

#endif // StressDivergenceGapContactSpring_H
