#pragma once

#include "CPKalidindiBackstressUpdate.h"
#include "GrainTrackerMatProp.h"
#include "FeatureMatPropVectorPostprocessor.h"
#include "FEProblem.h"

/*
   Constitutive model form:
   Reference 1: L. Zhao, P. Chakraborty, M.R. Tonks, I. Szlufarska, On the plastic driving force of grain boundary migration: A fully coupled phase field and crystal plasticity model, Computational Materials Science. 128 (2017) 320–330.
   Reference 2: L. Chen, J. Chen, R.A. Lebensohn, Y.Z. Ji, T.W. Heo, S. Bhattacharyya, K. Chang, S. Mathaudhu, Z.K. Liu, L.-Q. Chen, An integrated fast Fourier transform-based phase-field and crystal plasticity approach to model recrystallization of three dimensional polycrystals, Computer Methods in Applied Mechanics and Engineering. 285 (2015) 829–848.

   CP to PF: // TODO

   PF to CP: // TODO

   Elastic energy and Plasticity: // TODO
*/

class CPKalidindiBackstressCplUpdate;

class CPKalidindiBackstressCplUpdate : public CPKalidindiBackstressUpdate
{
public:
  static InputParameters validParams();

  CPKalidindiBackstressCplUpdate(const InputParameters & parameters);

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
