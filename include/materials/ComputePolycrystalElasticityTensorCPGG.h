// 目的：基于晶体塑性有限元模型 + grainTracker 自定义类来管理弹性模量
  // 问题1：弹性模量在整个过程中是否保持不变
// 参考1：poly_grain_growth_2D_eldrforce.i + ComputePolycrystalElasticityTensor + GrainTrackerElasticity
// 参考2：GrainPropertyReadFileCP + ComputeElasticityTensorCP

// step 1: 基于 ComputePolycrystalElasticityTensor 创建

#pragma once

#include "ComputeElasticityTensorBase.h"
#include "GrainDataTracker.h"
#include "RotationTensor.h"

class EulerAngleProvider;

class ComputePolycrystalElasticityTensorCPGG : public ComputeElasticityTensorBase
{
public:
  static InputParameters validParams();

  ComputePolycrystalElasticityTensorCPGG(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor();

  const GrainDataTracker<RankFourTensor> & _grain_tracker; /// Grain tracker object
  const unsigned int _op_num; /// Number of order parameters
  const std::vector<const VariableValue *> _vals; /// Order parameters
};