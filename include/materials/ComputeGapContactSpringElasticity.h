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

#ifndef ComputeGapContactSpringElasticity_H
#define ComputeGapContactSpringElasticity_H

#include "Material.h"
#include "ColumnMajorMatrix.h"

/**
 * ComputeGapContactSpringElasticity calculates the deformations,forces and stiffness matrix
 * of a contact spring element.
 **/

// Forward Declarations
class ComputeGapContactSpringElasticity;

template <>
InputParameters validParams<ComputeGapContactSpringElasticity>();

class ComputeGapContactSpringElasticity : public Material
{
public:
  ComputeGapContactSpringElasticity(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  void initQpStatefulProperties() override;

  /// compute transformation matrix for the contact spring element
  void computeTransformationMatrix();

  /// compute global deformations of the contact spring element
  void computeDeformation();

  /// initilize properties
  void initializeContactSpring();

  /// computes the basic forces and stiffness (along axis direction)
  void computeForces();

  /// convert basic forces and stiffness into global co-ordinate system
  void finalize();

  /// number of coupled displacement variables
  unsigned int _ndisp;

  /// variable numbers corresponding to the displacement variables
  std::vector<unsigned int> _disp_num;

  /// initial gap between the nodes
  MaterialProperty<Real> & _length;


  /// Gap between the nodes
  // const Real & _gap;

  /// contact normal stiffness
  const Real & _Kn;

  /// tangential contact stiffness
  const Real & _Kt;

  /// Coefficient of frction
  const Real & _mu;


  // /// relative tolerance for error in adaptive algorithm
  // const Real & _rel_tol;
  //
  // /// absolute tolerance for error in adaptive algorithm
  // const Real & _abs_tol;

  /// displacements in the contact spring local system, namely, deformations
  MaterialProperty<ColumnMajorMatrix> & _local_disp;

  /// displacements in the contact spring basic system, namely, deformations
  MaterialProperty<ColumnMajorMatrix> & _basic_def;

  /// displacements in the contact spring basic system, namely, old deformations
  const MaterialProperty<ColumnMajorMatrix> & _basic_def_old;

  // /// velocity in current step
  // MaterialProperty<Real> & _vel;
  //
  // /// velocity in previous step
  // const MaterialProperty<Real> & _vel_old;

  /// Spring forces in the basic coordinate system
  MaterialProperty<ColumnMajorMatrix> & _Fb;

  /// Spring forces of the previous step in the basic coordinate system
  const MaterialProperty<ColumnMajorMatrix> & _Fb_old;

  /// contact spring forces in the local coordinate system
  MaterialProperty<ColumnMajorMatrix> & _Fl;

  /// contact spring forces in the global coordinate system
  MaterialProperty<ColumnMajorMatrix> & _Fg;

  /// contact spring stiffness in the basic coordinate system
  MaterialProperty<ColumnMajorMatrix> & _Kb;

  /// contact spring stiffness in the local coordinate system
  MaterialProperty<ColumnMajorMatrix> & _Kl;

  /// contact spring stiffness in the global coordinate system
  MaterialProperty<ColumnMajorMatrix> & _Kg;

  /// transformation matrix from global coordinate system to contact spring local configuration at any time
  MaterialProperty<ColumnMajorMatrix> & _transform_gl;

  /// transformation matrix from contact spring local configuration to basic coordinate system at any time
  MaterialProperty<ColumnMajorMatrix> & _transform_lb;

  // const MaterialProperty<RealVectorValue> & _force;
};

#endif // ComputeGapContactSpringElasticity_H
