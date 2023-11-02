//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FeatureMatPropCplVectorPostprocessor.h"

registerMooseObject("qinglongApp", FeatureMatPropCplVectorPostprocessor);

InputParameters
FeatureMatPropCplVectorPostprocessor::validParams()
{
  InputParameters params = FeatureDataVectorPostprocessor::validParams();

  params.addParam<std::string>("base_name", "Name to append to reporters.");    
  params.addParam<unsigned int>("number_slip_systems", 12,
                                "The total number of possible active slip systems for the crystalline material");

  params.addClassDescription("This object is designed to pull information from the data structures "
                             "of a \"FeatureFloodCount\" or derived object (e.g. individual "
                             "feature volumes)"); // TODO   
  return params;
}

FeatureMatPropCplVectorPostprocessor::FeatureMatPropCplVectorPostprocessor(const InputParameters & parameters)
  : FeatureDataVectorPostprocessor(parameters),
    _number_slip_systems(getParam<unsigned int>("number_slip_systems")),
    _slip_resistances(declareRestartableData<std::vector<std::vector<Real>>>("slip_resistances")),    
    _backstresses(declareRestartableData<std::vector<std::vector<Real>>>("backstresses")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _slip_resistance_copy(getMaterialProperty<std::vector<Real>>(_base_name + "slip_resistance_copy")),    
    _backstress_copy(getMaterialProperty<std::vector<Real>>(_base_name + "backstress_copy"))
{
}

void
FeatureMatPropCplVectorPostprocessor::execute()
{
  const auto num_features = _feature_counter.getTotalFeatureCount();

  // Reset the state variables vector
  _slip_resistances.resize(num_features);
  _backstresses.resize(num_features);
  for (auto grain_index : make_range(num_features))
  {
    _slip_resistances[grain_index].assign(_number_slip_systems, 0.0);
    _backstresses[grain_index].assign(_number_slip_systems, 0.0);
  }

  FeatureDataVectorPostprocessor::execute();
}

void
FeatureMatPropCplVectorPostprocessor::accumulateVolumes(
    const Elem * elem,
    const std::vector<unsigned int> & var_to_features,
    std::size_t libmesh_dbg_var(num_features))
{
  unsigned int dominant_feature_id = FeatureFloodCount::invalid_id;
  Real max_var_value = std::numeric_limits<Real>::lowest();

  for (MooseIndex(var_to_features) var_index = 0; var_index < var_to_features.size(); ++var_index)
  {
    // Only sample "active" variables
    if (var_to_features[var_index] != FeatureFloodCount::invalid_id)
    {
      auto feature_id = var_to_features[var_index];
      mooseAssert(feature_id < num_features, "Feature ID out of range");
      auto integral_value = computeIntegral(var_index);

      // Compute volumes in a simplistic but domain conservative fashion
      if (_single_feature_per_elem)
      {
        if (integral_value > max_var_value)
        {
          // Update the current dominant feature and associated value
          max_var_value = integral_value;
          dominant_feature_id = feature_id;
        }
      }
      // Solution based volume calculation (integral value)
      else
        _feature_volumes[feature_id] += integral_value;

      std::vector<Real> slip_resistance_integral_value = computeStateVarivalesIntegral(var_index, _slip_resistance_copy);
      std::vector<Real> backstress_integral_value = computeStateVarivalesIntegral(var_index, _backstress_copy);

      for (auto sr_index : make_range(_number_slip_systems))
      {
        _slip_resistances[feature_id][sr_index] += slip_resistance_integral_value[sr_index];
        _backstresses[feature_id][sr_index] += backstress_integral_value[sr_index];
      }
    }
  }

  // Accumulate the entire element volume into the dominant feature. Do not use the integral value
  if (_single_feature_per_elem && dominant_feature_id != FeatureFloodCount::invalid_id)
    _feature_volumes[dominant_feature_id] += _assembly.elementVolume(elem);
}

std::vector<Real> 
FeatureMatPropCplVectorPostprocessor::computeStateVarivalesIntegral(std::size_t var_index, const MaterialProperty<std::vector<Real>> & state_variable_copy) const
{
  std::vector<Real> sum(_number_slip_systems, 0.0);

  for (auto sr_index : make_range(_number_slip_systems))
    for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
      sum[sr_index] += _JxW[qp] * _coord[qp] * (*_coupled_sln[var_index])[qp] * state_variable_copy[qp][sr_index];

  return sum;
}

void
FeatureMatPropCplVectorPostprocessor::finalize()
{
  // Do the parallel sum
  _communicator.sum(_feature_volumes);

  auto num_features = _feature_volumes.size();
  sum_state_variables(_slip_resistances, num_features);
  sum_state_variables(_backstresses, num_features);
}

void 
FeatureMatPropCplVectorPostprocessor::sum_state_variables(std::vector<std::vector<Real>> & stat_variables, const unsigned int & num_features)
{
  for (auto grain_index : make_range(num_features))
    _communicator.sum(stat_variables[grain_index]);

  for (auto grain_index : make_range(num_features))
  {
    if (_feature_volumes[grain_index] > 0.0)
      for (auto sr_index : make_range(_number_slip_systems))
        stat_variables[grain_index][sr_index] /= _feature_volumes[grain_index];
    else
      stat_variables[grain_index].assign(_number_slip_systems, 0.0);
  }
}

std::vector<Real> 
FeatureMatPropCplVectorPostprocessor::getStateVariable(unsigned int feature_id, const state_variable & state_variable_name) const //
{
  mooseAssert(feature_id < _slip_resistances.size(), "feature_id is out of range");

  switch (state_variable_name) 
  {
    case state_variable::slip_resistance:
        return _slip_resistances[feature_id];
    case state_variable::backstress:
        return _backstresses[feature_id];
  }

  return _backstresses[feature_id];
}