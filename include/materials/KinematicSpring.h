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


#pragma once


// Forward Declarations
class KinematicSpring;

template <>
InputParameters validParams<KinematicSpring>();

/**
 * KinematicSpring material simulates a plastic spring with linear kinematic hardening
 */
class KinematicSpring : public Material
{
public:

  KinematicSpring(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() ;
  /// Compute properties for each qp
  virtual void computeQpProperties() override;

  /// Computes the total rotation matrix accounting for rigid body rotations.
  /// Relevant for large rotation problems. Currently approximated as the
  /// initial rotation (t=0) matrix and the exact calculation of total rotation
  /// is not implemented.
  void computeTotalRotation();

  /// Computes the displacement and rotation deformations of the spring element
  void computeDeformations();

  /// Computes the spring forces in global coordinate system
  void computeForces();

  /// Computes stiffness matrix relating the forces and displacements in the
  /// global coordinate system
  void computeStiffnessMatrix();

  /// Computes the linear hardening variable based on Prager rule
  virtual Real computeHardeningValue(Real scalar, Real j);
  virtual Real computeHardeningMomentValue(Real scalar, Real j);

  /// Number of coupled rotational variables
  unsigned int _nrot;

  /// Number of coupled displacement variables
  unsigned int _ndisp;

  /// Variable numbers corresponding to the rotational variables
  std::vector<unsigned int> _rot_num;

  /// Variable numbers corresponding to the displacement variables
  std::vector<unsigned int> _disp_num;

  /// Deformations in the spring calculated at the currrent and previous time step in the spring local coordinates
  MaterialProperty<RealVectorValue> & _deformations;
  const MaterialProperty<RealVectorValue> & _deformations_old;

  /// Rotations in the spring calculated at current and previous timestep in the spring local coordinates
  MaterialProperty<RealVectorValue> & _rotations;
  const MaterialProperty<RealVectorValue> & _rotations_old;

  /// Axial stiffness of the spring (local x direction)
  const VariableValue & _kx;

  /// Shear stiffness of the spring in local y direction
  const VariableValue & _ky;

  /// Shear stiffness of the spring in local z direction
  const VariableValue & _kz;

  /// Torsional stiffness of the spring (local x direction)
  const VariableValue & _krx;

  /// Rotational stiffness of the spring in the local y direction
  const VariableValue & _kry;

  /// Rotational stiffness of the spring in the local z direction
  const VariableValue & _krz;

  /// Global displacement at node 0 of the spring element
  RealVectorValue _global_disp0; // node 0

  /// Global displacement at node 1 of the spring element
  RealVectorValue _global_disp1; // node 1

  /// Global rotation at node 0 of the spring element
  RealVectorValue _global_rot0; // node 0

  /// Global rotation at node 1 of the spring element
  RealVectorValue _global_rot1; // node 1

  /// Local displacement at node 0 of the spring element
  RealVectorValue _local_disp0; // node 0

  /// Local displacement at node 1 of the spring element
  RealVectorValue _local_disp1; // node 1

  /// Local rotation at node 0 of the spring element
  RealVectorValue _local_rot0; // node 0

  /// Local rotation at node 1 of the spring element
  RealVectorValue _local_rot1; // node 1

  /// Spring forces in the local coordinate system
  RealVectorValue _spring_forces_local;

  // Spring moments in the local coordinate system
  RealVectorValue _spring_moments_local;

  /// Spring displacement stiffness matrix in the local coordinate system
  RankTwoTensor _kdd_local;

  /// Spring rotational stiffness matrix in the local coordinate system
  RankTwoTensor _krr_local;

  /// Spring forces in the global coordinate system at current and previous timestep
  MaterialProperty<RealVectorValue> & _spring_forces_global;
  const MaterialProperty<RealVectorValue> & _spring_forces_global_old;

  /// Spring moments in the global coordinate system at current and previous timestep
  MaterialProperty<RealVectorValue> & _spring_moments_global;
  const MaterialProperty<RealVectorValue> & _spring_moments_global_old;

  /// Spring displacement stiffness matrix in the global coordinate system
  MaterialProperty<RankTwoTensor> & _kdd;

  /// Spring rotational stiffness matrix in the global coordinate system
  MaterialProperty<RankTwoTensor> & _krr;

  /// Rotational transformation from global coordinate system to spring local configuration at t = 0
  RankTwoTensor _original_global_to_local_rotation;

  /// Rotational transformation from global coordinate system to spring local configuration at any time
  MaterialProperty<RankTwoTensor> & _total_global_to_local_rotation;



 //  yield and hardening property input
 RealVectorValue _yield_force;
 RealVectorValue _yield_moments;
 const RealVectorValue & _hardening_constant_force;
 const RealVectorValue & _hardening_constant_moment;


 /// convergence tolerance
 Real _absolute_tolerance;

 /// Variable to store stiffness during iteration
 Real _k;

 /// Plastic displacements and rotation at current and previous time step
 MaterialProperty<RealVectorValue> & _plastic_deformation;
 const MaterialProperty<RealVectorValue> & _plastic_deformation_old;
 MaterialProperty<RealVectorValue> & _plastic_rotation;
 const MaterialProperty<RealVectorValue> & _plastic_rotation_old;

/// Hardening variables at current and previous time step
 MaterialProperty<RealVectorValue> & _hardening_variable;
 const MaterialProperty<RealVectorValue> & _hardening_variable_old;
 MaterialProperty<RealVectorValue> & _hardening_variable_moment;
 const MaterialProperty<RealVectorValue> & _hardening_variable_moment_old;


/// Elastic displacements and rotation at current and previous time step
 MaterialProperty<RealVectorValue> & _elastic_deformation;
 MaterialProperty<RealVectorValue> & _elastic_rotation;


 /// maximum no. of iterations
 const unsigned int _max_its;

};
