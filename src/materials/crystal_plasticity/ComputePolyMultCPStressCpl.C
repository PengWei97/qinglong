#include "ComputePolyMultCPStressCpl.h"

registerMooseObject("TensorMechanicsApp", ComputePolyMultCPStressCpl);

InputParameters
ComputePolyMultCPStressCpl::validParams()
{
  InputParameters params = ComputePolyMultCPStressCopy::validParams();
  
  params.addClassDescription(
      "Crystal Plasticity and phase field coupled model");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors");
  params.addParam<Real>("length_scale", 1.0e-9, "Length scale of the problem, in meters");
  params.addParam<Real>("pressure_scale", 1.0e6, "Pressure scale of the problem, in pa");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}
ComputePolyMultCPStressCpl::ComputePolyMultCPStressCpl(
    const InputParameters & parameters)
  : ComputePolyMultCPStressCopy(parameters),
  _elastic_energy_name(_base_name + "elastic_energy"),
  _elastic_energy(declareProperty<Real>(_elastic_energy_name)),    
  _grain_tracker(getUserObject<GrainTrackerMatProp>("grain_tracker")),
  _op_num(coupledComponents("v")),
  _vals(coupledValues("v")),
  _D_elastic_energy(_op_num),
  _length_scale(getParam<Real>("length_scale")),
  _pressure_scale(getParam<Real>("pressure_scale")),
  _JtoeV(6.24150974e18)
{
  _D_elastic_energy.resize(_op_num);
  // Loop over variables (ops)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    // declare elasticity tensor derivative properties
    _D_elastic_energy[op_index] = &declarePropertyDerivative<Real>(
        _elastic_energy_name, coupledName("v", op_index));
  }
}

void
ComputePolyMultCPStressCpl::initialSetup()
{
  ComputePolyMultCPStressCopy::initialSetup();
}

void
ComputePolyMultCPStressCpl::initQpStatefulProperties()
{
  ComputePolyMultCPStressCopy::initQpStatefulProperties();
  _elastic_energy[_qp] = 0.0;
}

void
ComputePolyMultCPStressCpl::postSolveQp(RankTwoTensor & cauchy_stress, RankFourTensor & jacobian_mult)
{
  ComputePolyMultCPStressCopy::postSolveQp(cauchy_stress, jacobian_mult);

  computeMechanicalEnergy(); 
}

void
ComputePolyMultCPStressCpl::computeMechanicalEnergy()
{
  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  // Calculate elastic energy
  _elastic_energy[_qp] = 0.5 * _total_lagrangian_strain[_qp].doubleContraction(_pk2[_qp]);

  std::cout << "_elastic_energy[_qp] " << _elastic_energy[_qp] << std::endl;
  Real sum_h = 0.0;
  unsigned int max_id = 0;
  Real max_vaule = (*_vals[max_id])[_qp];
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
  {
    (*_D_elastic_energy[op_index])[_qp] = 0.0;  
    
    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;
    sum_h += h;  

    if ((*_vals[op_index])[_qp] > max_vaule) 
    {
      max_id = op_index; 
      max_vaule = (*_vals[max_id])[_qp];
    }
  }

  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);

  // Calculate elastic energy derivative: Cderiv = dhdopi/sum_h * (Cop - _Cijkl)
  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    Real dhdopi = libMesh::pi * std::cos(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5)) / 2.0;
    Real & C_deriv = (*_D_elastic_energy[op_index])[_qp];
    
    Real elastic_energy_grain_id = _elastic_energy[_qp];
    // if (op_index != max_id)
    //   elastic_energy_grain_id = _grain_tracker.getElasticEnergy(grain_id);

    C_deriv = (elastic_energy_grain_id - _elastic_energy[_qp]) * dhdopi / sum_h;

    // Convert from XPa to eV/(xm)^3, where X is pressure scale and x is length scale;
    C_deriv *= _JtoeV * (_length_scale * _length_scale * _length_scale) * _pressure_scale;
  }  
}
