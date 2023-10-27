#include "ComputePolyElasticTensorCpl.h"
#include "RotationTensor.h"
#include "EulerAngles.h"

registerMooseObject("qinglongApp", ComputePolyElasticTensorCpl);

InputParameters
ComputePolyElasticTensorCpl::validParams()
{
  InputParameters params = ComputeElasticityTensorBase::validParams();
  params.addClassDescription(
      "Compute an evolving elasticity tensor coupled to a grain growth phase field and cp model.");
  params.addRequiredParam<UserObjectName>(
      "grain_tracker", "Name of GrainTracker user object that provides RankFourTensors");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

ComputePolyElasticTensorCpl::ComputePolyElasticTensorCpl(
    const InputParameters & parameters)
  : ComputeElasticityTensorBase(parameters),
  _grain_tracker(getUserObject<GrainTrackerMatProp>("grain_tracker")),
  _op_num(coupledComponents("v")),
  _vals(coupledValues("v")),
  _crysrot(declareProperty<RankTwoTensor>(_base_name + "crysrot"))
{ 
}

void
ComputePolyElasticTensorCpl::computeQpElasticityTensor()
{
  // Get list of active order parameters from grain tracker
  const auto & op_to_grains = _grain_tracker.getVarToFeatureVector(_current_elem->id());

  // Calculate elasticity tensor
  _elasticity_tensor[_qp].zero();
  Real sum_h = 0.0;
  unsigned int max_id = 0;
  Real max_var = (*_vals[max_id])[_qp];

  for (MooseIndex(op_to_grains) op_index = 0; op_index < op_to_grains.size(); ++op_index)
  {
    auto grain_id = op_to_grains[op_index];
    if (grain_id == FeatureFloodCount::invalid_id)
      continue;

    // Interpolation factor for elasticity tensors
    Real h = (1.0 + std::sin(libMesh::pi * ((*_vals[op_index])[_qp] - 0.5))) / 2.0;

    // Sum all rotated elasticity tensors
    _elasticity_tensor[_qp] += _grain_tracker.getElasticTensor(grain_id) * h;
    sum_h += h;

    if ((*_vals[op_index])[_qp] > max_var)
    {
      max_id = op_index;
      max_var = (*_vals[max_id])[_qp];
    }
  }

  const Real tol = 1.0e-10;
  sum_h = std::max(sum_h, tol);
  _elasticity_tensor[_qp] /= sum_h;

  auto grain_id = op_to_grains[max_id];
  auto & angles = _grain_tracker.getEulerAngles(max_id);

  RotationTensor R(angles);
  _crysrot[_qp] = R.transpose();
}