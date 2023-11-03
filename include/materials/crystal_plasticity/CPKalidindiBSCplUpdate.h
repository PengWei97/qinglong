#pragma once

#include "CPKalidindiBackstressUpdate.h"
#include "GrainTrackerMatProp.h"
#include "FeatureMatPropCplVectorPostprocessor.h"
#include "FEProblem.h"

typedef MaterialProperty<std::vector<Real>> MaterialPropVec;
typedef std::vector<std::vector<Real>> Matrix2D;

class CPKalidindiBSCplUpdate;

class CPKalidindiBSCplUpdate : public CPKalidindiBackstressUpdate
{
public:
  static InputParameters validParams();

  CPKalidindiBSCplUpdate(const InputParameters & parameters);

protected:

  virtual void initQpStatefulProperties() override;

  virtual void setInitialConstitutiveVariableValues() override;

  void convertStateVariablesFromPFtoCP();

  MaterialProperty<std::vector<Real>> & _slip_resistance_copy;
  MaterialProperty<std::vector<Real>> & _backstress_copy;

  MaterialProperty<unsigned int> & _number_active_grains;
  const MaterialProperty<unsigned int> & _number_active_grains_old;

  const GrainTrackerMatProp & _grain_tracker;
  const unsigned int _op_num;
  const std::vector<const VariableValue *> _vals;

  FEProblemBase & _fe_problem;
  int & _t_step;
  const VectorPostprocessorName & _vpp_name;
};