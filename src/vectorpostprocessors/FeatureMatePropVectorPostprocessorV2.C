/*
  Pengwei ECUST 2023.10.25 #25807: 
  In order to obtain average material parameters for each grain ID during parallel solving, including parameters of types Real, RankTwoTensor, and RankFourTensor. These parameters will subsequently be used for calculations in the Materials class.
*/ 

#include "FeatureMatePropVectorPostprocessorV2.h"

registerMooseObject("PhaseFieldApp", FeatureMatePropVectorPostprocessorV2);

InputParameters
FeatureMatePropVectorPostprocessorV2::validParams()
{
  InputParameters params = FeatureVolumeVectorPostprocessorCopy::validParams();

  params.addClassDescription("This object is designed to obtain the average value of the material parameters within each grain region.");
  params.addRequiredParam<MaterialPropertyName>("mat_prop", "The name of the material property");
  return params;
}

FeatureMatePropVectorPostprocessorV2::FeatureMatePropVectorPostprocessorV2(
    const InputParameters & parameters)
  : FeatureVolumeVectorPostprocessorCopy(parameters),
    _mat_prop_qp(getMaterialProperty<Real>("mat_prop")),
    _slip_resistance(getMaterialProperty<std::vector<Real>>("slip_resistance")),
    _test_real(getMaterialProperty<Real>("test_real")),
    _mat_prop_aves(declareVector("mat_prop_aves")),
    _grain_id_to_slip_resistances_map()
{
  std::cout << "the done of FeatureMatePropVectorPostprocessorV2" << std::endl;
}

void
FeatureMatePropVectorPostprocessorV2::initialize()
{
  std::cout << "FeatureMatePropVectorPostprocessorV2::initialize()" << std::endl;
  std::cout << "_test_real: " << _test_real[0] << std::endl;
  std::cout << "_slip_resistance["<< 0 << "].size() " << _slip_resistance[0].size() << std::endl;
}

void
FeatureMatePropVectorPostprocessorV2::execute()
{
  std::cout << "FeatureMatePropVectorPostprocessorV2::execute()" << std::endl;

  FeatureVolumeVectorPostprocessorCopy::execute();

  const auto num_features = _feature_counter.getTotalFeatureCount();
  _mat_prop_aves.assign(num_features, 0);

  // loop each element
  for (const auto & elem : _mesh.getMesh().active_local_element_ptr_range()) 
  {
    _fe_problem.setCurrentSubdomainID(elem, 0); // const Elem * elem
    _fe_problem.prepare(elem, 0);
    _fe_problem.reinitElem(elem, 0);    

    const auto & var_to_features = _feature_counter.getVarToFeatureVector(elem->id());

    accumulateMatProps(elem, var_to_features, num_features);
  }
}

void
FeatureMatePropVectorPostprocessorV2::accumulateMatProps(
    const Elem * elem,
    const std::vector<unsigned int> & var_to_features,
    std::size_t libmesh_dbg_var(num_features))
{
  unsigned int dominant_feature_id = FeatureFloodCount::invalid_id;

  for (MooseIndex(var_to_features) var_index = 0; var_index < var_to_features.size(); ++var_index)
  {
    if (var_to_features[var_index] != FeatureFloodCount::invalid_id)
    {
      auto feature_id = var_to_features[var_index];
      mooseAssert(feature_id < num_features, "Feature ID out of range");

      auto mat_integral_value = computeMatPropIntegral(var_index);
      _mat_prop_aves[feature_id] += mat_integral_value;
    }
  }  
}

void
FeatureMatePropVectorPostprocessorV2::finalize()
{
  FeatureVolumeVectorPostprocessorCopy::finalize();

  _communicator.sum(_mat_prop_aves);

  auto num_features = _feature_volumes.size();
  for (MooseIndex(num_features) feature_num = 0; feature_num < num_features; ++feature_num)
  {
    if (_feature_volumes[feature_num] > 0.0)
      _mat_prop_aves[feature_num] = _mat_prop_aves[feature_num]/_feature_volumes[feature_num];
    else
      _mat_prop_aves[feature_num] = 0.0;
  }
}

Real
FeatureMatePropVectorPostprocessorV2::getMatPropAves(unsigned int feature_id) const
{
  mooseAssert(feature_id < _mat_prop_aves.size(), "feature_id is out of range");
  return _mat_prop_aves[feature_id];
}

Real
FeatureMatePropVectorPostprocessorV2::computeMatPropIntegral(std::size_t var_index) const
{
  Real sum = 0;

  for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
    sum += _JxW[qp] * _coord[qp] * MetaPhysicL::raw_value(_mat_prop_qp[qp]) * (*_coupled_sln[var_index])[qp];
  
  return sum;
}


std::vector<Real>
FeatureMatePropVectorPostprocessorV2::getGrainIdtoSlipResistances(unsigned int feature_id)
{
  std::vector<Real> temp_vec(12, 0.0);

  return temp_vec;
}