//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ACGGElasticEnergyCpl.h"

#include "Material.h"

registerMooseObject("PhaseFieldApp", ACGGElasticEnergyCpl);

InputParameters
ACGGElasticEnergyCpl::validParams()
{
  InputParameters params = ACBulk<Real>::validParams();
  params.addClassDescription("Adds elastic energy contribution to the Allen-Cahn equation");
  params.addRequiredParam<MaterialPropertyName>(
      "D_elastic_energy_name", "The elastic energy derivative for the specific order parameter");
  return params;
}

ACGGElasticEnergyCpl::ACGGElasticEnergyCpl(const InputParameters & parameters)
  : ACBulk<Real>(parameters),
    _D_elastic_energy(getMaterialProperty<Real>("D_elastic_energy_name"))
{
}

Real
ACGGElasticEnergyCpl::computeDFDOP(PFFunctionType type)
{
  switch (type)
  {
    case Residual:
      return _D_elastic_energy[_qp]; // Compute the elastic energy driving force

    case Jacobian:
      return 0.0;
  }

  mooseError("Invalid type passed in");
}
