/*
  Pengwei ECUST 2023.10.25 #25807: 
  In order to obtain average material parameters for each grain ID during parallel solving, including parameters of types Real, RankTwoTensor, and RankFourTensor. These parameters will subsequently be used for calculations in the Materials class.
*/ 

#pragma once

#include "FeatureVolumeVectorPostprocessor.h"

class FeatureMatePropVectorPostprocessor : public FeatureVolumeVectorPostprocessor
{
public:
  static InputParameters validParams();

  FeatureMatePropVectorPostprocessor(const InputParameters & parameters);

  virtual void execute() override;
  virtual void finalize() override;    

  // Returns the average material property for the given grain ID.
  Real getMatPropAves(unsigned int feature_id) const;

protected:  
  // Get the material property of each [qp] in each element
  const MaterialProperty<Real> & _mat_prop_qp;

  // Calculate the average material parameters corresponding to each grain
  VectorPostprocessorValue & _mat_prop_aves;

  /// Add material property contributions to one or entries in the feature volume vector
  void accumulateMatProps(const Elem * elem,
                         const std::vector<unsigned int> & var_to_features,
                         std::size_t num_features);
  
  /// Calculate the integral value of material property
  Real computeMatPropIntegral(std::size_t var_index) const;
};