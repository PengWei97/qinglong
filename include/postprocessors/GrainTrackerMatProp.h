//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// 目的：基于GrainTracker管理grain id -> material properties
// material properties: euler angles, elastic tensor

#include "GrainTracker.h"
#include "RankFourTensor.h"
#include "EulerAngles.h"

class EulerAngleProvider;

class GrainTrackerMatProp : public GrainTracker
{
public:
  static InputParameters validParams();

  GrainTrackerMatProp(const InputParameters & parameters);

  const EulerAngles & getEulerAngles(unsigned int grain_id) const;
  const RankFourTensor & getElasticTensor(unsigned int grain_id) const;

protected:
  virtual void newGrainCreated(unsigned int new_grain_id);

  virtual EulerAngles newGrainEuler(unsigned int new_grain_id);

  virtual RankFourTensor newGrainElastic(unsigned int new_grain_id);

  std::vector<EulerAngles> & _euler_angles_vec;

  std::vector<RankFourTensor> & _C_ijkl_vec;

  /// generate random rotations when the Euler Angle provider runs out of data (otherwise error out)
  const bool _random_rotations;

  /// unrotated elasticity tensor
  RankFourTensor _C_ijkl;

  /// object providing the Euler angles
  const EulerAngleProvider & _euler;
};
