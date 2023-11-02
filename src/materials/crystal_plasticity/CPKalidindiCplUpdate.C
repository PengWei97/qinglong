//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CPKalidindiCplUpdate.h"
#include "libmesh/int_range.h"
#include "FEProblem.h"

registerMooseObject("TensorMechanicsApp", CPKalidindiCplUpdate);

InputParameters
CPKalidindiCplUpdate::validParams()
{
  InputParameters params = CrystalPlasticityStressUpdateBase::validParams();
  params.addClassDescription("Kalidindi version of homogeneous crystal plasticity.");
  params.addParam<Real>("r", 1.0, "Latent hardening coefficient");
  params.addParam<Real>("h", 541.5, "hardening constants");
  params.addParam<Real>("t_sat", 109.8, "saturated slip system strength");
  params.addParam<Real>("gss_a", 2.5, "coefficient for hardening");
  params.addParam<Real>("ao", 0.001, "slip rate coefficient");
  params.addParam<Real>("xm", 0.1, "exponent for slip rate");
  params.addParam<Real>("gss_initial", 60.8, "initial lattice friction strength of the material");

  params.addParam<MaterialPropertyName>(
      "total_twin_volume_fraction",
      "Total twin volume fraction, if twinning is considered in the simulation");

  params.addRequiredParam<VectorPostprocessorName>(
      "vpp_name", "The name of the VectorPostprocessor that you want to use");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");  
  return params;
}

CPKalidindiCplUpdate::CPKalidindiCplUpdate(const InputParameters & parameters)
  : CrystalPlasticityStressUpdateBase(parameters),
    // Constitutive values
    _r(getParam<Real>("r")),
    _h(getParam<Real>("h")),
    _tau_sat(getParam<Real>("t_sat")),
    _gss_a(getParam<Real>("gss_a")),
    _ao(getParam<Real>("ao")),
    _xm(getParam<Real>("xm")),
    _gss_initial(getParam<Real>("gss_initial")),

    // resize vectors used in the consititutive slip hardening
    _hb(_number_slip_systems, 0.0),
    _slip_resistance_increment(_number_slip_systems, 0.0),

    // resize local caching vectors used for substepping
    _previous_substep_slip_resistance(_number_slip_systems, 0.0),
    _slip_resistance_before_update(_number_slip_systems, 0.0),

    // Twinning contributions, if used
    _include_twinning_in_Lp(parameters.isParamValid("total_twin_volume_fraction")),
    _twin_volume_fraction_total(_include_twinning_in_Lp
                                    ? &getMaterialPropertyOld<Real>("total_twin_volume_fraction")
                                    : nullptr),

    _slip_resistance_copy(declareProperty<std::vector<Real>>(_base_name + "slip_resistance_copy")),
    _grain_tracker(getUserObject<GrainTrackerMatProp>("grain_tracker")),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _fe_problem(*parameters.get<FEProblem *>("_fe_problem")),
    _t_step(_fe_problem.timeStep()),
    _vpp_name(getParam<VectorPostprocessorName>("vpp_name")),
    _number_active_grains(declareProperty<unsigned int>("number_active_grains")),
    _number_active_grains_old(getMaterialPropertyOld<unsigned int>("number_active_grains"))
{
}

void
CPKalidindiCplUpdate::initQpStatefulProperties()
{
  CrystalPlasticityStressUpdateBase::initQpStatefulProperties();

  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance[_qp][i] = _gss_initial;
    _slip_increment[_qp][i] = 0.0;
  }

  _slip_resistance_copy[_qp] = _slip_resistance[_qp];
  _number_active_grains[_qp] = 1;
}

void
CPKalidindiCplUpdate::setInitialConstitutiveVariableValues()
{
  // Would also set old dislocation densities here if included in this model
  _slip_resistance[_qp] = _slip_resistance_old[_qp];
  _previous_substep_slip_resistance = _slip_resistance_old[_qp];

  if (_t_step > 0)
    convertStateVariablesFromPFtoCP();
}

// v2
void
CPKalidindiCplUpdate::convertStateVariablesFromPFtoCP()
{
  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  std::vector<unsigned int> active_op_indexs;
  unsigned int max_op_index = 0;
  Real max_op_value = (*_vals[max_op_index])[_qp];
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue; 

    active_op_indexs.push_back(op_index);
    if (max_op_value < (*_vals[op_index])[_qp])
    {
      max_op_index = op_index;
      max_op_value = (*_vals[max_op_index])[_qp];
    }
  }

  _number_active_grains[_qp] = active_op_indexs.size();

  if (_number_active_grains[_qp] == 1) // inside grain
  {
    if (_number_active_grains[_qp] == _number_active_grains_old[_qp])
      _slip_resistance[_qp] = _slip_resistance_old[_qp];
    else
      _slip_resistance[_qp].resize(_number_slip_systems, _gss_initial);
  }
  else if (_number_active_grains[_qp] > 1) // at GB or junction
  {
    std::vector<std::vector<Real>> grain_to_resistance_slips(active_op_indexs.size());

    // TODO： want need set it to be _vpp_object_ptr
    const FeatureMatPropVectorPostprocessor * vpp_object_ptr = dynamic_cast<FeatureMatPropVectorPostprocessor *>(const_cast<VectorPostprocessor *>(&_fe_problem.getVectorPostprocessorObjectByName(_vpp_name)));

    if (!vpp_object_ptr)
      mooseError("Pointer cast failed! Object is not of expected type");    
    
    if (_number_active_grains[_qp] > _number_active_grains_old[_qp]) // have some active grains
    {
      std::vector<Real> vector_active_op_val;
      for (auto op_index : active_op_indexs)
        vector_active_op_val.push_back((*_vals[op_index])[_qp]);

      std::sort(vector_active_op_val.rbegin(), vector_active_op_val.rend()); // descending sort

      unsigned int new_active_number = _number_active_grains[_qp] - _number_active_grains_old[_qp];
      Real critical_op_val = vector_active_op_val[new_active_number];

      for (unsigned int i = 0; i < active_op_indexs.size(); i++)
      {
        auto & op_index = active_op_indexs[i];
        Real op_val = (*_vals[op_index])[_qp];
        auto grain_id = op_to_grains[op_index];

        if (op_val >= critical_op_val && op_val == max_op_index)
          grain_to_resistance_slips[i] = _slip_resistance_old[_qp];
        else if (op_val >= critical_op_val)
          grain_to_resistance_slips[i] = vpp_object_ptr->getSlipResistance(grain_id);
        else
          grain_to_resistance_slips[i].resize(_number_slip_systems, _gss_initial);
      }
    }
    else // No new activated grains or reduce the number of activated grains
    {
      for (unsigned int i = 0; i < active_op_indexs.size(); i++)
      {
        auto & op_index = active_op_indexs[i];
        if (op_index != max_op_index) // second: _slip_resistance_ave
        {
          auto grain_id = op_to_grains[op_index];
          grain_to_resistance_slips[i] = vpp_object_ptr->getSlipResistance(grain_id);
        }
        else // primary
          grain_to_resistance_slips[i] = _slip_resistance_old[_qp];
      }
    }

    // Interpolate _slip_resistance[_qp]
    Real sum_h = 0.0;
    _slip_resistance[_qp].clear();
    std::vector<Real> sum_slip_resistance(_number_slip_systems, 0.0);

    for (unsigned int i = 0; i < active_op_indexs.size(); i++)
    {
      auto & op_index = active_op_indexs[i];
      
      Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;
      sum_h += h;

      for (unsigned int sr_index = 0; sr_index < _number_slip_systems; sr_index++)
        sum_slip_resistance[sr_index] += grain_to_resistance_slips[i][sr_index] * h;
    }

    const Real tol = 1.0e-10;
    sum_h = std::max(sum_h, tol);

    for (unsigned int sr_index = 0; sr_index < _number_slip_systems; sr_index++)
      _slip_resistance[_qp][sr_index] = sum_slip_resistance[sr_index] / sum_h;
  }
}

// v1
// void
// CPKalidindiCplUpdate::convertStateVariablesFromPFtoPF()
// {
//   // Get list of active order parameters from grain tracker
//   const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

//   unsigned int max_id = 0;
//   Real max_var = (*_vals[max_id])[_qp];
//   for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
//   {
//     auto grain_id = op_to_grains[op_index];
//     if (grain_id == FeatureFloodCount::invalid_id)
//       continue;  
    
//     if ((*_vals[op_index])[_qp] > max_var)
//     {
//       max_id = op_index;
//       max_var = (*_vals[max_id])[_qp];
//     }    
//   }

//   auto max_grain_id = op_to_grains[max_id];

//   // TODO want need set it to be _vpp_object_ptr
//   const FeatureMatPropVectorPostprocessor * vpp_object_ptr = dynamic_cast<FeatureMatPropVectorPostprocessor *>(const_cast<VectorPostprocessor *>(&_fe_problem.getVectorPostprocessorObjectByName(_vpp_name)));

//   if (!vpp_object_ptr)
//     mooseError("Pointer cast failed! Object is not of expected type");

//   std::vector<Real> max_slip_resistance = vpp_object_ptr->getSlipResistance(max_grain_id);

//   _slip_resistance[_qp] = max_slip_resistance;
//   _previous_substep_slip_resistance = max_slip_resistance;
// }

void
CPKalidindiCplUpdate::setSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _slip_resistance[_qp] = _previous_substep_slip_resistance;
}

bool
CPKalidindiCplUpdate::calculateSlipRate()
{
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_increment[_qp][i] =
        _ao * std::pow(std::abs(_tau[_qp][i] / _slip_resistance[_qp][i]), 1.0 / _xm);
    if (_tau[_qp][i] < 0.0)
      _slip_increment[_qp][i] *= -1.0;

    if (std::abs(_slip_increment[_qp][i]) * _substep_dt > _slip_incr_tol)
    {
      if (_print_convergence_message)
        mooseWarning("Maximum allowable slip increment exceeded ",
                     std::abs(_slip_increment[_qp][i]) * _substep_dt);

      return false;
    }
  }
  return true;
}

void
CPKalidindiCplUpdate::calculateEquivalentSlipIncrement(
    RankTwoTensor & equivalent_slip_increment)
{
  if (_include_twinning_in_Lp)
  {
    for (const auto i : make_range(_number_slip_systems))
      equivalent_slip_increment += (1.0 - (*_twin_volume_fraction_total)[_qp]) *
                                   _flow_direction[_qp][i] * _slip_increment[_qp][i] * _substep_dt;
  }
  else // if no twinning volume fraction material property supplied, use base class
    CrystalPlasticityStressUpdateBase::calculateEquivalentSlipIncrement(equivalent_slip_increment);
}

void
CPKalidindiCplUpdate::calculateConstitutiveSlipDerivative(
    std::vector<Real> & dslip_dtau)
{
  for (const auto i : make_range(_number_slip_systems))
  {
    if (MooseUtils::absoluteFuzzyEqual(_tau[_qp][i], 0.0))
      dslip_dtau[i] = 0.0;
    else
      dslip_dtau[i] = _ao / _xm *
                      std::pow(std::abs(_tau[_qp][i] / _slip_resistance[_qp][i]), 1.0 / _xm - 1.0) /
                      _slip_resistance[_qp][i];
  }
}

bool
CPKalidindiCplUpdate::areConstitutiveStateVariablesConverged()
{
  return isConstitutiveStateVariableConverged(_slip_resistance[_qp],
                                              _slip_resistance_before_update,
                                              _previous_substep_slip_resistance,
                                              _resistance_tol);
}

void
CPKalidindiCplUpdate::updateSubstepConstitutiveVariableValues()
{
  // Would also set substepped dislocation densities here if included in this model
  _previous_substep_slip_resistance = _slip_resistance[_qp];
}

void
CPKalidindiCplUpdate::cacheStateVariablesBeforeUpdate()
{
  _slip_resistance_before_update = _slip_resistance[_qp];
}

void
CPKalidindiCplUpdate::calculateStateVariableEvolutionRateComponent()
{
  for (const auto i : make_range(_number_slip_systems))
  {
    // Clear out increment from the previous iteration
    _slip_resistance_increment[i] = 0.0;

    _hb[i] = _h * std::pow(std::abs(1.0 - _slip_resistance[_qp][i] / _tau_sat), _gss_a);
    const Real hsign = 1.0 - _slip_resistance[_qp][i] / _tau_sat;
    if (hsign < 0.0)
      _hb[i] *= -1.0;
  }

  for (const auto i : make_range(_number_slip_systems))
  {
    for (const auto j : make_range(_number_slip_systems))
    {
      unsigned int iplane, jplane;
      iplane = i / 3;
      jplane = j / 3;

      if (iplane == jplane) // self vs. latent hardening
        _slip_resistance_increment[i] +=
            std::abs(_slip_increment[_qp][j]) * _hb[j]; // q_{ab} = 1.0 for self hardening
      else
        _slip_resistance_increment[i] +=
            std::abs(_slip_increment[_qp][j]) * _r * _hb[j]; // latent hardenign
    }
  }
}

bool
CPKalidindiCplUpdate::updateStateVariables()
{
  // Now perform the check to see if the slip system should be updated
  for (const auto i : make_range(_number_slip_systems))
  {
    _slip_resistance_increment[i] *= _substep_dt;
    if (_previous_substep_slip_resistance[i] < _zero_tol && _slip_resistance_increment[i] < 0.0)
      _slip_resistance[_qp][i] = _previous_substep_slip_resistance[i];
    else
      _slip_resistance[_qp][i] =
          _previous_substep_slip_resistance[i] + _slip_resistance_increment[i];

    if (_slip_resistance[_qp][i] < 0.0)
      return false;
  }
  return true;
}
