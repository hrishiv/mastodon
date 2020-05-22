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
class BilinearSpring;

template <>
InputParameters validParams<BilinearSpring>();

/**
 * BilinearSpring material simulates a linear spring with a diagonal stiffness
 * matrix, including rotational stiffnesses.
 */
class BilinearSpring : public Material
{
public:

  BilinearSpring(const InputParameters & parameters);

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

  // virtual Real computeHardeningValue(Real scalar);

  /// Number of coupled rotational variables
  unsigned int _nrot;

  /// Number of coupled displacement variables
  unsigned int _ndisp;

  /// Variable numbers corresponding to the rotational variables
  std::vector<unsigned int> _rot_num;

  /// Variable numbers corresponding to the displacement variables
  std::vector<unsigned int> _disp_num;

  /// Deformations in the spring calculated at a time t in the spring local coordinates
  MaterialProperty<RealVectorValue> & _deformations;
  const MaterialProperty<RealVectorValue> & _deformations_old;

  /// Rotations in the spring calculated at a time t in the spring local coordinates
  MaterialProperty<RealVectorValue> & _rotations;

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

  /// Spring forces in the global coordinate system
  MaterialProperty<RealVectorValue> & _spring_forces_global;
  const MaterialProperty<RealVectorValue> & _spring_forces_global_old;

  MaterialProperty<RealVectorValue> & _spring_forces_test;
  const MaterialProperty<RealVectorValue> & _spring_forces_test_old;

  /// Spring moments in the global coordinate system
  MaterialProperty<RealVectorValue> & _spring_moments_global;

  /// Spring displacement stiffness matrix in the global coordinate system
  MaterialProperty<RankTwoTensor> & _kdd;

  /// Spring rotational stiffness matrix in the global coordinate system
  MaterialProperty<RankTwoTensor> & _krr;

  /// Rotational transformation from global coordinate system to spring local configuration at t = 0
  RankTwoTensor _original_global_to_local_rotation;

  /// Rotational transformation from global coordinate system to spring local configuration at any time
  MaterialProperty<RankTwoTensor> & _total_global_to_local_rotation;



  //  yield stress and hardening property input
  const Real _yield_force;
  const Real _hardening_constant;
  const Real _th_disp;
  const Real _deterioration_constant;


  /// convergence tolerance
  Real _absolute_tolerance;
  Real _relative_tolerance;


  MaterialProperty<Real> & _plastic_deformation;
  const MaterialProperty<Real> & _plastic_deformation_old;


  MaterialProperty<Real> & _hardening_variable; //back_trial
  const MaterialProperty<Real> & _hardening_variable_old; //back_trial_old
  MaterialProperty<Real> & _hardening_variable_pos; //back_trial_pos
  MaterialProperty<Real> & _hardening_variable_neg; //back_trial_neg
  MaterialProperty<Real> & _hard_var_carry; //back_trial_carry
  MaterialProperty<Real> & _deterioration_variable; //back_trial_d
  const MaterialProperty<Real> & _deterioration_variable_old; //back_trial_old_d old
  MaterialProperty<Real> & _hard_var_test; //back_trial_test
  const MaterialProperty<Real> & _hard_var_test_old; //back_trial_test_old
  MaterialProperty<Real> & _det_var_test; //back_trial_d_test
  const MaterialProperty<Real> & _det_var_test_old; //back_trial__d_test_old
  MaterialProperty<Real> & _eff_trial_force;
  const MaterialProperty<Real> & _eff_trial_force_old;

  MaterialProperty<Real> & _elastic_deformation;

  /// maximum no. of iterations
  const unsigned int _max_its;

  MaterialProperty<Real> & _yf; //updated yield force to be used for each time step
  MaterialProperty<Real> & _yf_pos; //updated yield force for positive direction
  MaterialProperty<Real> & _yf_neg; //updated yield force for negative direction
  MaterialProperty<Real> & _th; //updated threshold displacement to be used for each time step
  MaterialProperty<Real> & _th_disp_pos; //updated threshold displacement for positive direction
  MaterialProperty<Real> & _th_disp_neg; //updated threshold displacement for negative direction
  MaterialProperty<Real> & _hc; //updated hardening constant to be used for each time step
  MaterialProperty<Real> & _hc_pos; //updated hardening constant for positive direction
  MaterialProperty<Real> & _hc_neg; //updated hardening constant for negative direction
  MaterialProperty<Real> & _shift;  //vertical shift of yield surface
  MaterialProperty<Real> & _max;    //max positive displacement reached
  MaterialProperty<Real> & _min;    //max negative displacement reached

  /// Variable initialized to check if different conditions exist or not
  MaterialProperty<bool> & _isdamage; //damage wrt previous time step calculated in every time step
//
  MaterialProperty<bool> & _isdamagepos; //damage pos
//  const MaterialProperty<bool> & _isdamagepos_old;
  MaterialProperty<bool> & _isdamageneg; //damage neg
//  const MaterialProperty<bool> & _isdamageneg_old;
  MaterialProperty<bool> & _isdetpos; //has positive side k deteriorated
//  const MaterialProperty<bool> & _isdetpos_old;
  MaterialProperty<bool> & _isdetneg; //has negative side k deteriorated
//  const MaterialProperty<bool> & _isdetneg_old;
  MaterialProperty<bool> & _istest; //only one side deteriorated
//  const MaterialProperty<bool> & _istest_old;
  MaterialProperty<bool> & _istest2; //shift the yield surface as only one side deteriorated and the displacement
  ///doesnot reach the previous max/min displacment
//  const MaterialProperty<bool> & _istest2_old;
  MaterialProperty<bool> & _istest3pos;
//  const MaterialProperty<bool> & _istest3pos_old;
  MaterialProperty<bool> & _istest3neg;
//  const MaterialProperty<bool> & _istest3neg_old;
  MaterialProperty<bool> & _istemp; //initialize to store the back trial for only one side deteriorated
//  const MaterialProperty<bool> & _temp_old;


// Real _yf; //updated yield force to be used for each time step
// Real _yf_pos; //updated yield force for positive direction
// Real _yf_neg; //updated yield force for negative direction
// Real _shift;  //vertical shift of yield surface
// Real _max;
// Real _min;
//
//   /// Variable initialized to check if different conditions exist or not
// bool _isdamage; //damage wrt previous time step calculated in every time step
// bool _isdamagepos; //damage pos
// bool _isdamageneg; //damage neg
// bool _isdetpos; //has positive side k deteriorated
// bool _isdetneg; //has negative side k deteriorated
// bool _istest;
// bool _istest2;
// bool _istest3pos;
// bool _istest3neg;
// bool _temp;


};
