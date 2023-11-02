//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FeatureDataVectorPostprocessor.h"

enum class state_variable {slip_resistance, backstress};

class FeatureMatPropCplVectorPostprocessor : public FeatureDataVectorPostprocessor
{
public:
  static InputParameters validParams();

  FeatureMatPropCplVectorPostprocessor(const InputParameters & parameters);  

  virtual void execute() override;
  virtual void finalize() override;

  // _slip_resistance, _backstress
  std::vector<Real> getStateVariable(unsigned int feature_id, const state_variable & state_variable_name  = state_variable::slip_resistance) const; //

protected:
  virtual void accumulateVolumes(const Elem * elem,
                         const std::vector<unsigned int> & var_to_features,
                         std::size_t num_features) override;

  std::vector<Real> computeStateVarivalesIntegral(std::size_t var_index, const MaterialProperty<std::vector<Real>> & state_variable_copy) const;

  void sum_state_variables(std::vector<std::vector<Real>> & stat_variables, const unsigned int & num_features);

  // Get from material class CPKalindindiCplUpdate
  const MaterialProperty<std::vector<Real>> & _slip_resistance_copy;
  const MaterialProperty<std::vector<Real>> & _backstress_copy;

  const std::string _base_name;
  const unsigned int _number_slip_systems;

  std::vector<std::vector<Real>> & _slip_resistances;
  std::vector<std::vector<Real>> & _backstresses;
};

