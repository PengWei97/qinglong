// 目的：基于晶体塑性有限元模型 + grainTracker 自定义类来管理弹性模量
  // 问题1：弹性模量在整个过程中是否保持不变
// 参考1：poly_grain_growth_2D_eldrforce.i + ComputePolycrystalElasticityTensor + GrainTrackerElasticity
// 参考2：GrainPropertyReadFileCP + ComputeElasticityTensorCP

// step 1: 基于 ComputePolycrystalElasticityTensor 创建
// 

#pragma once

// ComputePolyElasticTensorCpl ~ Compute Polycrystal Elastic Tensor in Coupled

#include "ComputeElasticityTensorBase.h"
#include "GrainTrackerMatProp.h"

class EulerAngleProvider;

class ComputePolyElasticTensorCpl : public ComputeElasticityTensorBase
{
public:
  static InputParameters validParams();

  ComputePolyElasticTensorCpl(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor();

  /// vector of elasticity tensor material properties
  std::vector<MaterialProperty<RankFourTensor> *> _D_elastic_tensor;

  const GrainTrackerMatProp & _grain_tracker;

  const unsigned int _op_num; /// Number of order parameters
  const std::vector<const VariableValue *> _vals; /// Order parameters

  /// Crystal Rotation Matrix used to rotate the slip system direction and normal
  MaterialProperty<RankTwoTensor> & _crysrot;

  Real _length_scale;
  Real _pressure_scale;  
  const Real _JtoeV; /// Conversion factor from J to eV3
};