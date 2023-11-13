//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DeformedGrainEBSDMaterial.h"
#include "GrainTrackerInterface.h"

registerMooseObject("qinglongApp", DeformedGrainEBSDMaterial);

InputParameters
DeformedGrainEBSDMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale in m, where default is nm");
  params.addParam<Real>("time_scale", 1.0e-6, "Time scale in sec, where default is micro sec");
  params.addParam<Real>("Elas_Mod", 2.50e10, "Elastic Modulus in J/m^3");
  params.addParam<Real>("Burg_vec", 3.0e-10, "Length of Burger Vector in m");
  params.addParam<Real>("stored_factor", 0.5, "the scaling factor in stored energy function"); 
  params.addParam<bool>("concurrent_recovery", false, "The concurrent recovery would be considered if true");

  params.addParam<Real>("rho_end_l2", 2.10e12, "Dislocation density of grains after long-term concurrent recovery in level 2");
  params.addParam<Real>("rho_end_l3", 1.8e13, "Dislocation density of grains after long-term concurrent recovery in level 3");
  params.addParam<Real>("rho_critical", 3.9e13, "Critual dislocation density of grains");
  params.addParam<Real>("a_rho_l2", 4.6e-4, "Evolution coefficient of dislocation during medium time recovery at level 2"); 
  params.addParam<Real>("a_rho_l3", 6.0e-4, "Evolution coefficient of dislocation during medium time recovery at level 3");

  params.addRequiredParam<UserObjectName>("grain_tracker",
                                          "The GrainTracker UserObject to get values from.");
  params.addRequiredParam<UserObjectName>("GNDs_provider", "GNDs provider for EBSD reader");
  return params;
}

DeformedGrainEBSDMaterial::DeformedGrainEBSDMaterial(const InputParameters & parameters)
  : Material(parameters),
    _op_num(coupledComponents("v")),
    _vals(coupledValues("v")),
    _length_scale(getParam<Real>("length_scale")),
    _time_scale(getParam<Real>("time_scale")),
    _Elas_Mod(getParam<Real>("Elas_Mod")),
    _Burg_vec(getParam<Real>("Burg_vec")),
    _stored_factor(getParam<Real>("stored_factor")),
    _JtoeV(6.24150974e18),
    _concurrent_recovery(getParam<bool>("concurrent_recovery")),
    _beta(declareProperty<Real>("beta")),
    _rho_eff(declareProperty<Real>("rho_eff")),

    _rho_end_l2(getParam<Real>("rho_end_l2")),
    _rho_end_l3(getParam<Real>("rho_end_l3")),
    _rho_critical(getParam<Real>("rho_critical")),
    _a_rho_l2(getParam<Real>("a_rho_l2")),
    _a_rho_l3(getParam<Real>("a_rho_l3")),

    _grain_tracker(getUserObject<GrainTrackerInterface>("grain_tracker")),
    _GNDs_provider(getUserObject<EBSDReader>("GNDs_provider"))
{
  if (_op_num == 0)
    paramError("op_num", "Model requires op_num > 0");
}

void
DeformedGrainEBSDMaterial::computeQpProperties()
{

  Real SumEtai2 = 0.0;
  for (unsigned int i = 0; i < _op_num; ++i)
    SumEtai2 += (*_vals[i])[_qp] * (*_vals[i])[_qp];

  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  Real rho_i;

  _rho_eff[_qp] = 0.0;
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (op_to_grains[op_index] == FeatureFloodCount::invalid_id)
      continue;

    rho_i = getGNDsFromEBSD(grain_id);

    _rho_eff[_qp] += rho_i * (*_vals[op_index])[_qp] * (*_vals[op_index])[_qp];
  }

  _rho_eff[_qp] /= SumEtai2;

  _beta[_qp] = _stored_factor * _Elas_Mod * _Burg_vec * _Burg_vec * _JtoeV * _length_scale;
}

Real
DeformedGrainEBSDMaterial::getGNDsFromEBSD(const unsigned int & grain_id)
{  
  auto & time_current = _fe_problem.time(); // current simulation time s

  Real rho_init = 2.0e15 * (_length_scale * _length_scale); // TODO - need to be more physically based
  Real rho_i = rho_init;

  if (grain_id < _GNDs_provider.getGrainNum())
    rho_init = _GNDs_provider.getAvgData(grain_id)._custom[0] * (_length_scale * _length_scale); // GNDs for each grain, 1/m^2 

  // // 550du
  // const Real rho_critical = 4.0e13 * (_length_scale * _length_scale); // Real rho_yellow_up = 90;4.0e13
  // const Real rho_end_2_550 = 2.3e13 * (_length_scale * _length_scale); // the minimum GND due to recovery 1/m^2 3.0e7
  // const Real a_rho_2_550 = 1.0e-4;
  // const Real rho_end_3_550 = 8.5e13 * (_length_scale * _length_scale); // the minimum GND due to recovery 1/m^2 3.0e7
  // const Real a_rho_3_550 = 1.0e-4;
    
  // if (_concurrent_recovery) //  && (time_current > 50.0)
  // {
  //   if ((rho_init > rho_end_2_550) && (rho_init < rho_critical)) // level 2
  //     rho_i = (rho_init - rho_end_2_550) * std::exp(-a_rho_2_550 * time_current) + rho_end_2_550; // func2 - level 2
  //   else if (rho_init >= rho_end_3_550) // level 3
  //     rho_i = (rho_init - rho_end_3_550) * std::exp(-a_rho_3_550 * time_current) + rho_end_3_550; // func1 - level 3
  //   else // internal + lower
  //     rho_i = rho_init;
  // }

  // // // 700du V5
  // // const Real rho_critical_blue = 1.8e13 * (_length_scale * _length_scale); // Real rho_yellow_up = 90;
  // // const Real rho_critical_yellow = 3.9e13 * (_length_scale * _length_scale); // Real rho_yellow_up = 90;
  // // const Real rho_end_2 = 2.10e12 * (_length_scale * _length_scale); // the minimum GND due to recovery 1/m^2 3.0e7
  // // const Real a_rho_2 = 4.6e-4;
  // // const Real rho_end_3 = 3.9e13 * (_length_scale * _length_scale); // the minimum GND due to recovery 1/m^2 3.0e7
  // // const Real a_rho_3 = 6.0e-4;
  
  // // if (time_current > 21600)
  // // {
  // //   Real rho_init_700du = rho_i;
  // //   auto time_current_700 = time_current - 21600;
  // //   if (_concurrent_recovery) // && (time_current > 50.0)
  // //   {
  // //     if ((rho_init_700du > rho_end_2) && (rho_init_700du < rho_critical_blue)) // blue grain
  // //       rho_i = (rho_init_700du - rho_end_2) * std::exp(-a_rho_2 * time_current_700) + rho_end_2; // func2 - level 2
  // //     else if (rho_init_700du >= rho_critical_yellow) // yellow grain
  // //       rho_i = (rho_init_700du - rho_end_3) * std::exp(-a_rho_3 * time_current_700) + rho_end_3; // func1 - level 3
  // //     else
  // //       rho_i = rho_init_700du;
  // //   }
  // // }

  // // 700du V6
  // const Real rho_critical_700 = 4.0e13 * (_length_scale * _length_scale); // Real rho_yellow_up = 90;4.0e13
  // const Real rho_end_2_700 = 5.30e12 * (_length_scale * _length_scale); // the minimum GND due to recovery 1/m^2 3.0e7
  // const Real a_rho_2_700 = 4.6e-4;
  // const Real rho_end_3_700 = 4.0e13 * (_length_scale * _length_scale); // the minimum GND due to recovery 1/m^2 3.0e7
  // const Real a_rho_3_700 = 6.0e-4;
    
  // Real rho_init_700du = rho_i;
  // auto time_current_700 = time_current - 21600;  
  // if (_concurrent_recovery && (time_current_700 > 0.0)) //  && (time_current > 50.0)
  // {
  //   if ((rho_init_700du > rho_end_2_700) && (rho_init_700du < rho_critical_700)) // level 2
  //     rho_i = (rho_init_700du - rho_end_2_700) * std::exp(-a_rho_2_700 * time_current_700) + rho_end_2_700; // func2 - level 2
  //   else if (rho_init_700du >= rho_end_3_700) // level 3
  //     rho_i = (rho_init_700du - rho_end_3_700) * std::exp(-a_rho_3_550 * time_current_700) + rho_end_3_700; // func1 - level 3
  //   else // internal + lower
  //     rho_i = rho_init_700du;
  // }

  // 700du
  const Real rho_end_l2 = _rho_end_l2 * (_length_scale * _length_scale);
  const Real rho_critical = _rho_critical * (_length_scale * _length_scale);
  const Real rho_end_l3 = _rho_end_l3 * (_length_scale * _length_scale);
  
  if (_concurrent_recovery) // && (time_current > 50.0)
  {
    if ((rho_init > rho_end_l2) && (rho_init < _rho_critical)) // blue grain
      rho_i = (rho_init - rho_end_l2) * std::exp(-_a_rho_l2 * time_current) + rho_end_l2; // func2 - level 2
    else if (rho_init >= rho_end_l3) // yellow grain
      rho_i = (rho_init - rho_end_l3) * std::exp(-_a_rho_l3 * time_current) + rho_end_l3; // func1 - level 3
    else
      rho_i = rho_init;
  }

  // // 700du
  // const Real rho_critical_blue = 1.8e13 * (_length_scale * _length_scale); // Real rho_yellow_up = 90;
  // const Real rho_critical_yellow = 3.9e13 * (_length_scale * _length_scale); // Real rho_yellow_up = 90;
  // const Real rho_end_2 = 2.10e12 * (_length_scale * _length_scale); // the minimum GND due to recovery 1/m^2 3.0e7
  // const Real a_rho_2 = 4.6e-4;
  // const Real rho_end_3 = 3.9e13 * (_length_scale * _length_scale); // the minimum GND due to recovery 1/m^2 3.0e7
  // const Real a_rho_3 = 6.0e-4;
  
  // if (_concurrent_recovery) // && (time_current > 50.0)
  // {
  //   if ((rho_init > rho_end_2) && (rho_init < rho_critical_blue)) // blue grain
  //     rho_i = (rho_init - rho_end_2) * std::exp(-a_rho_2 * time_current) + rho_end_2; // func2 - level 2
  //   else if (rho_init >= rho_critical_yellow) // yellow grain
  //     rho_i = (rho_init - rho_end_3) * std::exp(-a_rho_3 * time_current) + rho_end_3; // func1 - level 3
  //   else
  //     rho_i = rho_init;
  // }

  return rho_i;
}