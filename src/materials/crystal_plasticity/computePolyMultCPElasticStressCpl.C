//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "computePolyMultCPElasticStressCpl.h"

#include "CrystalPlasticityStressUpdateBase.h"
#include "libmesh/utility.h"
#include "Conversion.h"
#include "MooseException.h"

registerMooseObject("TensorMechanicsApp", computePolyMultCPElasticStressCpl);

InputParameters
computePolyMultCPElasticStressCpl::validParams()
{
  InputParameters params = ComputeFiniteStrainElasticStress::validParams();

  params.addClassDescription(
      "Crystal Plasticity base class: handles the Newton iteration over the stress residual and "
      "calculates the Jacobian based on constitutive laws from multiple material classes "
      "that are inherited from CrystalPlasticityStressUpdateBase");
  params.addParam<std::string>(
      "base_name",
      "Optional parameter that allows the user to define multiple mechanics material systems on "
      "the same block, i.e. for multiple phases");

  params.addRequiredParam<std::vector<MaterialName>>(
      "crystal_plasticity_models",
      "The material objects to use to calculate crystal plasticity stress and strains.");
  params.addParam<std::vector<MaterialName>>("eigenstrain_names",
                                             "The material objects to calculate eigenstrains.");
  params.addParam<MooseEnum>("tan_mod_type",
                             MooseEnum("exact none", "none"),
                             "Type of tangent moduli for preconditioner: default elastic");
  params.addParam<Real>("rtol", 1e-6, "Constitutive stress residual relative tolerance");
  params.addParam<Real>("abs_tol", 1e-6, "Constitutive stress residual absolute tolerance");
  params.addParam<unsigned int>("maxiter", 100, "Maximum number of iterations for stress update");
  params.addParam<unsigned int>(
      "maxiter_state_variable", 100, "Maximum number of iterations for state variable update");
  params.addParam<unsigned int>(
      "maximum_substep_iteration", 1, "Maximum number of substep iteration");
  params.addParam<bool>("use_line_search", false, "Use line search in constitutive update");
  params.addParam<Real>("min_line_search_step_size", 0.01, "Minimum line search step size");
  params.addParam<Real>("line_search_tol", 0.5, "Line search bisection method tolerance");
  params.addParam<unsigned int>(
      "line_search_maxiter", 20, "Line search bisection method maximum number of iteration");
  params.addParam<MooseEnum>("line_search_method",
                             MooseEnum("CUT_HALF BISECTION", "CUT_HALF"),
                             "The method used in line search");
  params.addParam<bool>(
      "print_state_variable_convergence_error_messages",
      false,
      "Whether or not to print warning messages from the crystal plasticity specific convergence "
      "checks on the stress measure and general constitutive model quantinties.");
  return params;
}

computePolyMultCPElasticStressCpl::computePolyMultCPElasticStressCpl(
    const InputParameters & parameters)
  : ComputeFiniteStrainElasticStress(parameters),
    _num_models(getParam<std::vector<MaterialName>>("crystal_plasticity_models").size()),
    _num_eigenstrains(getParam<std::vector<MaterialName>>("eigenstrain_names").size()),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _rtol(getParam<Real>("rtol")),
    _abs_tol(getParam<Real>("abs_tol")),
    _maxiter(getParam<unsigned int>("maxiter")),
    _maxiterg(getParam<unsigned int>("maxiter_state_variable")),
    _tan_mod_type(getParam<MooseEnum>("tan_mod_type").getEnum<TangentModuliType>()),
    _max_substep_iter(getParam<unsigned int>("maximum_substep_iteration")),
    _use_line_search(getParam<bool>("use_line_search")),
    _min_line_search_step_size(getParam<Real>("min_line_search_step_size")),
    _line_search_tolerance(getParam<Real>("line_search_tol")),
    _line_search_max_iterations(getParam<unsigned int>("line_search_maxiter")),
    _line_search_method(getParam<MooseEnum>("line_search_method").getEnum<LineSearchMethod>()),
    _plastic_deformation_gradient(declareProperty<RankTwoTensor>("plastic_deformation_gradient")),
    _plastic_deformation_gradient_old(
        getMaterialPropertyOld<RankTwoTensor>("plastic_deformation_gradient")),
    _eigenstrain_deformation_gradient(
        _num_eigenstrains ? &declareProperty<RankTwoTensor>("eigenstrain_deformation_gradient")
                          : nullptr),
    _eigenstrain_deformation_gradient_old(
        _num_eigenstrains
            ? &getMaterialPropertyOld<RankTwoTensor>("eigenstrain_deformation_gradient")
            : nullptr),
    _deformation_gradient(getMaterialProperty<RankTwoTensor>(_base_name + "deformation_gradient")),
    _deformation_gradient_old(
        getMaterialPropertyOld<RankTwoTensor>(_base_name + "deformation_gradient")),
    _pk2(declareProperty<RankTwoTensor>("second_piola_kirchhoff_stress")),
    _pk2_old(getMaterialPropertyOld<RankTwoTensor>("second_piola_kirchhoff_stress")),
    _total_lagrangian_strain(
        declareProperty<RankTwoTensor>("total_lagrangian_strain")), // Lagrangian strain
    _updated_rotation(declareProperty<RankTwoTensor>("updated_rotation")),
    _crysrot(getMaterialProperty<RankTwoTensor>(
        _base_name + "crysrot")), // defined in the elasticity tensor classes for crystal plasticity
    _print_convergence_message(getParam<bool>("print_state_variable_convergence_error_messages"))
{
  _convergence_failed = false;
}

void
computePolyMultCPElasticStressCpl::initQpStatefulProperties()
{
  _plastic_deformation_gradient[_qp].zero();
  _plastic_deformation_gradient[_qp].addIa(1.0);

  _pk2[_qp].zero();

  _total_lagrangian_strain[_qp].zero();

  _updated_rotation[_qp].zero();
  _updated_rotation[_qp].addIa(1.0);
}

void
computePolyMultCPElasticStressCpl::computeQpStress()
{
  // Does not support face/boundary material property calculation
  if (isBoundaryMaterial())
    return;

  _temporary_deformation_gradient_old = _deformation_gradient_old[_qp];
  if (_temporary_deformation_gradient_old.det() == 0)
    _temporary_deformation_gradient_old.addIa(1.0);

  _delta_deformation_gradient = _deformation_gradient[_qp] - _temporary_deformation_gradient_old;

  preSolveQp();

  _substep_dt = _dt;

  _temporary_deformation_gradient = _delta_deformation_gradient;
  _temporary_deformation_gradient += _temporary_deformation_gradient_old;

  solveQp();

  postSolveQp(_stress[_qp], _Jacobian_mult[_qp]);
}

void
computePolyMultCPElasticStressCpl::preSolveQp()
{
  _pk2[_qp] = _pk2_old[_qp];
  _inverse_plastic_deformation_grad_old = _plastic_deformation_gradient_old[_qp].inverse();
}

void
computePolyMultCPElasticStressCpl::solveQp()
{
  _inverse_plastic_deformation_grad = _inverse_plastic_deformation_grad_old;

  solveStress();

  _plastic_deformation_gradient[_qp] =
      _inverse_plastic_deformation_grad.inverse(); // the postSoveStress

  // save off the old F^{p} inverse now that have converged on the stress and state variables
  _inverse_plastic_deformation_grad_old = _inverse_plastic_deformation_grad;
}

void
computePolyMultCPElasticStressCpl::postSolveQp(RankTwoTensor & cauchy_stress,
                                                    RankFourTensor & jacobian_mult)
{
  cauchy_stress = _elastic_deformation_gradient * _pk2[_qp] *
                  _elastic_deformation_gradient.transpose() / _elastic_deformation_gradient.det();

  calcTangentModuli(jacobian_mult);

  _total_lagrangian_strain[_qp] =
      _deformation_gradient[_qp].transpose() * _deformation_gradient[_qp] -
      RankTwoTensor::Identity();
  _total_lagrangian_strain[_qp] = _total_lagrangian_strain[_qp] * 0.5;

  // Calculate crystal rotation to track separately
  RankTwoTensor rot;
  _elastic_deformation_gradient.getRUDecompositionRotation(rot);
  _updated_rotation[_qp] = rot * _crysrot[_qp];
}

void
computePolyMultCPElasticStressCpl::solveStress()
{
  unsigned int iteration = 0;
  RankTwoTensor dpk2;
  Real rnorm, rnorm0, rnorm_prev;

  // Calculate stress residual
  calculateResidual();
  calculateJacobian();

  rnorm = _residual_tensor.L2norm();
  rnorm0 = rnorm;

  // Calculate stress increment
  dpk2 = -_jacobian.invSymm() * _residual_tensor;
  _pk2[_qp] = _pk2[_qp] + dpk2;

  calculateResidual(); // _residual_tensor
  calculateJacobian(); // _jacobian

  rnorm_prev = rnorm;
  rnorm = _residual_tensor.L2norm();
}

void
computePolyMultCPElasticStressCpl::calculateResidual()
{
  RankTwoTensor elastic_strain, pk2_new;

  _inverse_plastic_deformation_grad = _inverse_plastic_deformation_grad_old;

  _elastic_deformation_gradient = _temporary_deformation_gradient *
                                  _inverse_plastic_deformation_grad;

  elastic_strain = (_elastic_deformation_gradient.transpose() * _elastic_deformation_gradient - RankTwoTensor::Identity()) * 0.5;

  pk2_new = _elasticity_tensor[_qp] * elastic_strain;
  _residual_tensor = _pk2[_qp] - pk2_new;
}

void
computePolyMultCPElasticStressCpl::calculateJacobian()
{
  // may not need to cache the dfpinvdpk2 here. need to double check
  RankFourTensor dfedfpinv, deedfe, dfpinvdpk2, dfpinvdpk2_per_model;

  RankTwoTensor ffeiginv = _temporary_deformation_gradient;

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
        dfedfpinv(i, j, k, j) = ffeiginv(i, k);

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _elastic_deformation_gradient(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _elastic_deformation_gradient(k, i) * 0.5;
      }

  dfpinvdpk2.zero();
  _jacobian =
      RankFourTensor::IdentityFour() - (_elasticity_tensor[_qp] * deedfe * dfedfpinv * dfpinvdpk2);
}

void
computePolyMultCPElasticStressCpl::calcTangentModuli(RankFourTensor & jacobian_mult)
{
  switch (_tan_mod_type)
  {
    case TangentModuliType::EXACT:
      elastoPlasticTangentModuli(jacobian_mult);
      break;
    default:
      elasticTangentModuli(jacobian_mult);
  }
}

void
computePolyMultCPElasticStressCpl::elastoPlasticTangentModuli(RankFourTensor & jacobian_mult)
{
  RankFourTensor tan_mod;
  RankTwoTensor pk2fet, fepk2, feiginvfpinv;
  RankFourTensor deedfe, dsigdpk2dfe, dfedf;

  // Fill in the matrix stiffness material property
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto k : make_range(Moose::dim))
      {
        deedfe(i, j, k, i) = deedfe(i, j, k, i) + _elastic_deformation_gradient(k, j) * 0.5;
        deedfe(i, j, k, j) = deedfe(i, j, k, j) + _elastic_deformation_gradient(k, i) * 0.5;
      }

  usingTensorIndices(i_, j_, k_, l_);
  dsigdpk2dfe = _elastic_deformation_gradient.times<i_, k_, j_, l_>(_elastic_deformation_gradient) *
                _elasticity_tensor[_qp] * deedfe;

  pk2fet = _pk2[_qp] * _elastic_deformation_gradient.transpose();
  fepk2 = _elastic_deformation_gradient * _pk2[_qp];

  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto l : make_range(Moose::dim))
      {
        tan_mod(i, j, i, l) += pk2fet(l, j);
        tan_mod(i, j, j, l) += fepk2(i, l);
      }

  tan_mod += dsigdpk2dfe;

  const auto je = _elastic_deformation_gradient.det();
  if (je > 0.0)
    tan_mod /= je;

  feiginvfpinv = _inverse_plastic_deformation_grad;
  for (const auto i : make_range(Moose::dim))
    for (const auto j : make_range(Moose::dim))
      for (const auto l : make_range(Moose::dim))
        dfedf(i, j, i, l) = feiginvfpinv(l, j);

  jacobian_mult = tan_mod * dfedf;
}

void
computePolyMultCPElasticStressCpl::elasticTangentModuli(RankFourTensor & jacobian_mult)
{
  // update jacobian_mult
  jacobian_mult = _elasticity_tensor[_qp];
}

bool
computePolyMultCPElasticStressCpl::lineSearchUpdate(const Real & rnorm_prev,
                                                         const RankTwoTensor & dpk2)
{
  if (_line_search_method == LineSearchMethod::CutHalf)
  {
    Real rnorm;
    Real step = 1.0;

    do
    {
      _pk2[_qp] = _pk2[_qp] - step * dpk2;
      step /= 2.0;
      _pk2[_qp] = _pk2[_qp] + step * dpk2;

      calculateResidual();
      rnorm = _residual_tensor.L2norm();
    } while (rnorm > rnorm_prev && step > _min_line_search_step_size);

    // has norm improved or is the step still above minumum search step size?
    return (rnorm <= rnorm_prev || step > _min_line_search_step_size);
  }
  else if (_line_search_method == LineSearchMethod::Bisection)
  {
    unsigned int count = 0;
    Real step_a = 0.0;
    Real step_b = 1.0;
    Real step = 1.0;
    Real s_m = 1000.0;
    Real rnorm = 1000.0;

    calculateResidual();
    auto s_b = _residual_tensor.doubleContraction(dpk2);
    const auto rnorm1 = _residual_tensor.L2norm();
    _pk2[_qp] = _pk2[_qp] - dpk2;
    calculateResidual();
    auto s_a = _residual_tensor.doubleContraction(dpk2);
    const auto rnorm0 = _residual_tensor.L2norm();
    _pk2[_qp] = _pk2[_qp] + dpk2;

    if ((rnorm1 / rnorm0) < _line_search_tolerance || s_a * s_b > 0)
    {
      calculateResidual();
      return true;
    }

    while ((rnorm / rnorm0) > _line_search_tolerance && count < _line_search_max_iterations)
    {
      _pk2[_qp] = _pk2[_qp] - step * dpk2;
      step = 0.5 * (step_b + step_a);
      _pk2[_qp] = _pk2[_qp] + step * dpk2;
      calculateResidual();
      s_m = _residual_tensor.doubleContraction(dpk2);
      rnorm = _residual_tensor.L2norm();

      if (s_m * s_a < 0.0)
      {
        step_b = step;
        s_b = s_m;
      }
      if (s_m * s_b < 0.0)
      {
        step_a = step;
        s_a = s_m;
      }
      count++;
    }

    // below tolerance and max iterations?
    return ((rnorm / rnorm0) < _line_search_tolerance && count < _line_search_max_iterations);
  }
  else
    mooseError("Line search method is not provided.");
}

void
computePolyMultCPElasticStressCpl::calculateEigenstrainDeformationGrad()
{
  _inverse_eigenstrain_deformation_grad.zero();
  _inverse_eigenstrain_deformation_grad.addIa(1.0);

  for (unsigned int i = 0; i < _num_eigenstrains; ++i)
  {
    _eigenstrains[i]->setSubstepDt(_substep_dt);
    _eigenstrains[i]->computeQpProperties();
    _inverse_eigenstrain_deformation_grad *= _eigenstrains[i]->getDeformationGradientInverse();
  }
  (*_eigenstrain_deformation_gradient)[_qp] = _inverse_eigenstrain_deformation_grad.inverse();
}
