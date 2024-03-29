//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVTurbulentTemperatureWallFunction.h"
#include "NavierStokesMethods.h"
#include "NS.h"

registerADMooseObject("NavierStokesApp", INSFVTurbulentTemperatureWallFunction);

InputParameters
INSFVTurbulentTemperatureWallFunction::validParams()
{
  InputParameters params = FVFluxBC::validParams();

  params.addClassDescription("Adds turbulent temperature wall function.");
  params.addRequiredParam<MooseFunctorName>("T_w", "The wall temperature.");
  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Density");
  params.addRequiredParam<MooseFunctorName>(NS::mu, "Dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>(NS::cp, "The spcific heat at constant pressure.");
  params.addParam<MooseFunctorName>(NS::kappa, "The thermal conductivity.");
  params.addParam<MooseFunctorName>("Pr_t", 0.58, "The turbulent Prandtl number.");
  params.addParam<bool>(
      "linearized_yplus",
      false,
      "Boolean to indicate if yplus must be estimate locally for the blending functions.");
  return params;
}

INSFVTurbulentTemperatureWallFunction::INSFVTurbulentTemperatureWallFunction(
    const InputParameters & parameters)
  : FVFluxBC(parameters),
    _dim(_subproblem.mesh().dimension()),
    _T_w(getFunctor<ADReal>("T_w")),
    _u_var(getFunctor<ADReal>("u")),
    _v_var(parameters.isParamValid("v") ? &(getFunctor<ADReal>("v")) : nullptr),
    _w_var(parameters.isParamValid("w") ? &(getFunctor<ADReal>("w")) : nullptr),
    _rho(getFunctor<ADReal>(NS::density)),
    _mu(getFunctor<ADReal>(NS::mu)),
    _cp(getFunctor<ADReal>(NS::cp)),
    _kappa(getFunctor<ADReal>(NS::kappa)),
    _Pr_t(getFunctor<ADReal>("Pr_t")),
    _linearized_yplus(getParam<bool>("linearized_yplus"))
{
}

ADReal
INSFVTurbulentTemperatureWallFunction::computeQpResidual()
{
  // Useful parameters
  const FaceInfo & fi = *_face_info;
  const Real wall_dist = std::abs((fi.elemCentroid() - fi.faceCentroid()) * fi.normal());
  const Elem & _current_elem = fi.elem();
  const auto current_argument = makeElemArg(&_current_elem);
  const auto state = determineState();
  const auto mu = _mu(current_argument, state);
  const auto rho = _rho(current_argument, state);
  const auto cp = _cp(current_argument, state);
  const auto kappa = _kappa(current_argument, state);

  // Get the velocity vector
  ADRealVectorValue velocity(_u_var(current_argument, state));
  if (_v_var)
    velocity(1) = (*_v_var)(current_argument, state);
  if (_w_var)
    velocity(2) = (*_w_var)(current_argument, state);

  // Compute the velocity and direction of the velocity component that is parallel to the wall
  const ADReal parallel_speed = (velocity - velocity * (fi.normal()) * (fi.normal())).norm();

  // Computing friction velocity
  ADReal u_tau;
  if (_linearized_yplus)
  {
    const ADReal a_c = 1 / NS::von_karman_constant;
    const ADReal b_c =
        1 / NS::von_karman_constant * (std::log(NS::E_turb_constant * wall_dist / mu) + 1.0);
    const ADReal c_c = parallel_speed;
    u_tau = (-b_c + std::sqrt(std::pow(b_c, 2) + 4.0 * a_c * c_c)) / (2.0 * a_c);
  }
  else
    u_tau = NS::findUStar(mu, rho, parallel_speed, wall_dist);

  // Computing non-dimensional wall distance
  const ADReal y_plus = wall_dist * u_tau * rho / mu;

  ADReal alpha;
  if (y_plus <= 5.0) // sub-laminar layer
  {
    alpha = kappa / (rho * cp);
  }
  else if (y_plus >= 30.0) // log-layer
  {
    const auto Pr = cp * mu / kappa;
    const auto Pr_ratio = Pr / _Pr_t(current_argument, state);
    const auto jayatilleke_P =
        9.24 * (std::pow(Pr_ratio, 0.75) - 1.0) * (1.0 + 0.28 * std::exp(-0.007 * Pr_ratio));
    const auto wall_scaling =
        1.0 / NS::von_karman_constant * std::log(NS::E_turb_constant * y_plus) + jayatilleke_P;
    alpha = u_tau * wall_dist / (_Pr_t(current_argument, state) * wall_scaling);
  }
  else // buffer layer
  {
    const auto alpha_lam = kappa / (rho * cp);
    const auto Pr = cp * mu / kappa;
    const auto Pr_ratio = Pr / _Pr_t(current_argument, state);
    const auto jayatilleke_P =
        9.24 * (std::pow(Pr_ratio, 0.75) - 1.0) * (1.0 + 0.28 * std::exp(-0.007 * Pr_ratio));
    const auto wall_scaling =
        1.0 / NS::von_karman_constant * std::log(NS::E_turb_constant * y_plus) + jayatilleke_P;
    const auto alpha_turb = u_tau * wall_dist / (_Pr_t(current_argument, state) * wall_scaling);
    const auto blending_function = (y_plus - 5.0) / 25.0;
    alpha = blending_function * alpha_turb + (1.0 - blending_function) * alpha_lam;
  }

  const auto face_arg = singleSidedFaceArg();
  return -rho * cp * alpha * (_T_w(face_arg, state) - _var(current_argument, state)) / wall_dist;
}
