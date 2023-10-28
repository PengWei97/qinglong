#pragma once

#include "ComputePolyMultCPStressCopy.h"
#include "GrainTrackerMatProp.h"

class ComputePolyMultCPStressCpl : public ComputePolyMultCPStressCopy
{
public:
  static InputParameters validParams();

  ComputePolyMultCPStressCpl(const InputParameters & parameters);

  virtual void initialSetup() override;

protected:
  virtual void initQpStatefulProperties() override;

  virtual void postSolveQp(RankTwoTensor & cauchy_stress, RankFourTensor & jacobian_mult) override;  

  virtual void computeMechanicalEnergy();

  std::string _elastic_energy_name;
  MaterialProperty<Real> & _elastic_energy;
  std::vector<MaterialProperty<Real> *> _D_elastic_energy;
  
  const GrainTrackerMatProp & _grain_tracker; /// Grain tracker object 
  const unsigned int _op_num; /// Number of order parameters
  const std::vector<const VariableValue *> _vals; /// Order parameters  

  Real _length_scale;
  Real _pressure_scale;
  const Real _JtoeV; /// Conversion factor from J to eV

  const VectorPostprocessorValue & _mat_prop_vect; 
};