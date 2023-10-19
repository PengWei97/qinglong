#include "ComputePolycrystalCrystalPlasticityStress.h"


registerMooseObject("TensorMechanicsApp", ComputePolycrystalCrystalPlasticityStress);

InputParameters
ComputePolycrystalCrystalPlasticityStress::validParams()
{
  InputParameters params = ComputeMultipleCrystalPlasticityStress::validParams();

  params.addClassDescription("xxx");
  return params;
}

ComputePolycrystalCrystalPlasticityStress::ComputePolycrystalCrystalPlasticityStress(
    const InputParameters & parameters)
  : ComputeMultipleCrystalPlasticityStress(parameters)
{
}

void
ComputePolycrystalCrystalPlasticityStress::computeQpStress()
{
  ComputeMultipleCrystalPlasticityStress::computeQpStress();
}