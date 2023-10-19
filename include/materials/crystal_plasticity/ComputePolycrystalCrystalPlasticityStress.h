#pragma once

#include "ComputeMultipleCrystalPlasticityStress.h"

class ComputePolycrystalCrystalPlasticityStress : public ComputeMultipleCrystalPlasticityStress
{
public:
  static InputParameters validParams();

  ComputePolycrystalCrystalPlasticityStress(const InputParameters & parameters);

protected:
  virtual void computeQpStress() override;
};
