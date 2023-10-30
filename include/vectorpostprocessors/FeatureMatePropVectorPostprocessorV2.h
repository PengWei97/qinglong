/*
  Pengwei ECUST 2023.10.25 #25807: 
  In order to obtain average material parameters for each grain ID during parallel solving, including parameters of types Real, RankTwoTensor, and RankFourTensor. These parameters will subsequently be used for calculations in the Materials class.
*/ 

#pragma once

#include "FeatureVolumeVectorPostprocessorCopy.h"
#include "RankTwoTensor.h"

class FeatureMatePropVectorPostprocessorV2 : public FeatureVolumeVectorPostprocessorCopy
{
public:
  static InputParameters validParams();

  FeatureMatePropVectorPostprocessorV2(const InputParameters & parameters);

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;    

  // Returns the average material property for the given grain ID.
  Real getMatPropAves(unsigned int feature_id) const;

  std::vector<Real> getGrainIdtoSlipResistances(unsigned int feature_id);

protected:  
  // Get the material property of each [qp] in each element
  const MaterialProperty<Real> & _mat_prop_qp;
  // Calculate the average material parameters corresponding to each grain
  VectorPostprocessorValue & _mat_prop_aves;

  const MaterialProperty<Real> & _test_real;
  const MaterialProperty<std::vector<Real>> & _slip_resistance;

  std::map<unsigned int, std::vector<Real>> _grain_id_to_slip_resistances_map;

  /// Add material property contributions to one or entries in the feature volume vector
  void accumulateMatProps(const Elem * elem,
                         const std::vector<unsigned int> & var_to_features,
                         std::size_t num_features);
  
  /// Calculate the integral value of material property
  Real computeMatPropIntegral(std::size_t var_index) const;
};