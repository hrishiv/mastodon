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
class FuelCompressionSpring;

template <>
InputParameters validParams<FuelCompressionSpring>();

/**
 * FuelCompressionSpring material simulates a linear spring with a diagonal stiffness
 * matrix, including rotational stiffnesses.
 */
class FuelCompressionSpring : public Material
{
public:

  FuelCompressionSpring(const InputParameters & parameters);

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
  const MaterialProperty<RealVectorValue> & _deformations_older;

  /// Rotations in the spring calculated at a time t in the spring local coordinates
  MaterialProperty<RealVectorValue> & _rotations;

  /// Axial stiffness of the spring (local x direction)
  const VariableValue & _k_trans;


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
  MaterialProperty<RealVectorValue> &_spring_forces_local;
  const MaterialProperty<RealVectorValue> &_spring_forces_local_old;
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
  Real _hardening_constant;
  Real _alpha_s;
  Real _recovery_factor;
  const Real _th_disp;
  // Real _deterioration_constant;
  // Real _alpha_c;
  Real _k_sg;
  bool _Mahin_mod;
  bool _complete_det;  //To account for the complete permanent deformation after exceeding
                      //the threshold where the spring flattens
  const Real _component;
  const Real _tension_factor;


  /// convergence tolerance
  Real _absolute_tolerance;
  Real _relative_tolerance;


  MaterialProperty<Real> & _plastic_deformation;
  const MaterialProperty<Real> & _plastic_deformation_old;


  MaterialProperty<Real> & _hardening_variable; //back_trial
  const MaterialProperty<Real> & _hardening_variable_old; //back_trial_old
  MaterialProperty<Real> & _deterioration_variable; //back_trial_d
  const MaterialProperty<Real> & _deterioration_variable_old; //back_trial_old_d old
  MaterialProperty<Real> & _net_def;
  const MaterialProperty<Real> & _net_def_old;

  // MaterialProperty<Real> & _temp;
  // const MaterialProperty<Real> & _temp_old;

  MaterialProperty<Real> & _elastic_deformation;

  /// maximum no. of iterations
  const unsigned int _max_its;
  Real _min_c;

  MaterialProperty<Real> & _k;
  const MaterialProperty<Real> & _k_old;
  //Real  _k;
  MaterialProperty<Real> & _yf; //updated yield force to be used for each time step
  const MaterialProperty<Real> & _yf_old;
  MaterialProperty<Real> & _yf_pos; //updated yield force for positive direction
  const MaterialProperty<Real> & _yf_pos_old;
  MaterialProperty<Real> & _yf_neg; //updated yield force for negative direction
  const MaterialProperty<Real> & _yf_neg_old;
  MaterialProperty<Real> & _th; //updated threshold displacement to be used for each time step
  const MaterialProperty<Real> & _th_old;
  MaterialProperty<Real> & _th_disp_pos; //updated threshold displacement for positive direction
  const MaterialProperty<Real> & _th_disp_pos_old;
  MaterialProperty<Real> & _th_disp_neg; //updated threshold displacement for negative direction
  const MaterialProperty<Real> & _th_disp_neg_old;
  MaterialProperty<Real> & _hc; //updated hardening constant to be used for each time step
  const MaterialProperty<Real> & _hc_old;
  MaterialProperty<Real> & _hc_mod; //updated hardening constant to adjust the mahin modification
    const MaterialProperty<Real> & _hc_mod_old;
  // MaterialProperty<Real> & _hc_pos; //updated hardening constant for positive direction
  // MaterialProperty<Real> & _hc_neg; //updated hardening constant for negative direction
  MaterialProperty<Real> & _max_d;    //max positive displacement
    const MaterialProperty<Real> & _max_d_old;
  MaterialProperty<Real> & _min_d;    //max negative displacement
  const MaterialProperty<Real> & _min_d_old;
  MaterialProperty<Real> & _max_f;    //max positive force
    const MaterialProperty<Real> & _max_f_old;
  MaterialProperty<Real> & _min_f;    //max negative force
  const MaterialProperty<Real> & _min_f_old;
  MaterialProperty<Real> & _prev_pos_d;    //prev cycle positive displacement
  const MaterialProperty<Real> & _prev_pos_d_old;
  MaterialProperty<Real> & _prev_neg_d;    //prev cycle negative displacement
  const MaterialProperty<Real> & _prev_neg_d_old;
  MaterialProperty<Real> & _prev_pos_f;    //prev cycle positive force
  const MaterialProperty<Real> & _prev_pos_f_old;
  MaterialProperty<Real> & _prev_neg_f;    //prev cycle negative force
  const MaterialProperty<Real> & _prev_neg_f_old;

  MaterialProperty<bool> & _isdet_start; //set to true at the first step of the deterioration curve
  const MaterialProperty<bool> & _isdet_start_old;
  MaterialProperty<bool> & _isplastic;
  const MaterialProperty<bool> & _isplastic_old;
  // /// Property that stores the time step limit
  // MaterialProperty<Real> & _matl_timestep_limit;

  MaterialProperty<Real> & _length_org;

  MaterialProperty<Real> & _length_new;

};
