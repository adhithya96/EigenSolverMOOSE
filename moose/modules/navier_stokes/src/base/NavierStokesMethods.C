//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NavierStokesMethods.h"
#include "MooseError.h"
#include "libmesh/vector_value.h"
#include "NS.h"

namespace NS
{
int
delta(unsigned int i, unsigned int j)
{
  if (i == j)
    return 1;
  else
    return 0;
}

int
computeSign(const Real & a)
{
  return a > 0 ? 1 : (a < 0 ? -1 : 0);
}

unsigned int
getIndex(const Real & p, const std::vector<Real> & bounds)
{
  if (p < bounds.front() || p > bounds.back())
    mooseError("Point exceeds bounds of domain!");

  for (unsigned int i = 1; i < bounds.size(); ++i)
    if (p <= bounds[i])
      return i - 1;

  return bounds.size() - 2;
}

Real
reynoldsPropertyDerivative(
    const Real & Re, const Real & rho, const Real & mu, const Real & drho, const Real & dmu)
{
  return Re * (drho / std::max(rho, 1e-6) - dmu / std::max(mu, 1e-8));
}

Real
prandtlPropertyDerivative(const Real & mu,
                          const Real & cp,
                          const Real & k,
                          const Real & dmu,
                          const Real & dcp,
                          const Real & dk)
{
  return (k * (mu * dcp + cp * dmu) - mu * cp * dk) / std::max(k * k, 1e-8);
}

ADReal
findUStar(const ADReal & mu, const ADReal & rho, const ADReal & u, const Real dist)
{
  // usually takes about 3-4 iterations
  constexpr int MAX_ITERS{50};
  constexpr Real REL_TOLERANCE{1e-6};

  const ADReal nu = mu / rho;

  ADReal u_star = std::sqrt(nu * u / dist);

  // Newton-Raphson method to solve for u_star (friction velocity).
  for (int i = 0; i < MAX_ITERS; ++i)
  {
    ADReal residual = u_star / NS::von_karman_constant * std::log(u_star * dist / (0.111 * nu)) - u;
    ADReal deriv = (1 + std::log(u_star * dist / (0.111 * nu))) / NS::von_karman_constant;
    ADReal new_u_star = std::max(1e-20, u_star - residual / deriv);

    Real rel_err = std::abs((new_u_star.value() - u_star.value()) / new_u_star.value());

    u_star = new_u_star;
    if (rel_err < REL_TOLERANCE)
      return u_star;
  }

  mooseException("Could not find the wall friction velocity (mu: ",
                 mu,
                 " rho: ",
                 rho,
                 " velocity: ",
                 u,
                 " wall distance: ",
                 dist,
                 ")");
}

ADReal
findyPlus(const ADReal & mu, const ADReal & rho, const ADReal & u, const Real dist)
{
  // Fixed point iteration method to find y_plus
  // It should take 3 or 4 iterations
  constexpr int MAX_ITERS{10};
  constexpr Real REL_TOLERANCE{1e-2};

  ADReal yPlusLast = 0.0;
  ADReal yPlus = dist * u * rho / mu; // Assign intitial value to laminar
  const ADReal rev_yPlusLam = 1.0 / yPlus;
  const ADReal kappa_time_Re = NS::von_karman_constant * u * dist / (mu / rho);
  unsigned int iters = 0;

  do
  {
    yPlusLast = yPlus;
    yPlus = (kappa_time_Re + yPlus) / (1.0 + std::log(NS::E_turb_constant * yPlus));
  } while (rev_yPlusLam * (yPlus - yPlusLast) > REL_TOLERANCE && ++iters < MAX_ITERS);

  return std::max(0.0, yPlus);
}

ADReal
computeSpeed(const ADRealVectorValue & velocity)
{
  // if the velocity is zero, then the norm function call fails because AD tries to calculate the
  // derivatives which causes a divide by zero - because d/dx(sqrt(f(x))) = 1/2/sqrt(f(x))*df/dx.
  // So add a bit of noise (based on hitchhiker's guide to the galaxy's meaning of life number) to
  // avoid this failure mode.
  return isZero(velocity) ? 1e-42 : velocity.norm();
}

/// Bounded element maps for wall treatement
void
getWallBoundedElements(const std::vector<BoundaryName> & wall_boundary_name,
                       const FEProblemBase & fe_problem,
                       const SubProblem & subproblem,
                       const std::set<SubdomainID> & block_ids,
                       std::map<const Elem *, bool> & wall_bounded_map)
{

  wall_bounded_map.clear();

  for (const auto & elem : fe_problem.mesh().getMesh().active_element_ptr_range())
  {
    if (block_ids.find(elem->subdomain_id()) != block_ids.end())
      for (const auto i_side : elem->side_index_range())
      {
        const auto & side_bnds = subproblem.mesh().getBoundaryIDs(elem, i_side);
        for (const auto & name : wall_boundary_name)
        {
          const auto wall_id = subproblem.mesh().getBoundaryID(name);
          for (const auto side_id : side_bnds)
            if (side_id == wall_id)
              wall_bounded_map[elem] = true;
        }
      }
  }
}

/// Bounded element face distances for wall treatement
void
getWallDistance(const std::vector<BoundaryName> & wall_boundary_name,
                const FEProblemBase & fe_problem,
                const SubProblem & subproblem,
                const std::set<SubdomainID> & block_ids,
                std::map<const Elem *, std::vector<Real>> & dist_map)
{

  dist_map.clear();

  for (const auto & elem : fe_problem.mesh().getMesh().active_element_ptr_range())
    if (block_ids.find(elem->subdomain_id()) != block_ids.end())
      for (const auto i_side : elem->side_index_range())
      {
        const auto & side_bnds = subproblem.mesh().getBoundaryIDs(elem, i_side);
        for (const auto & name : wall_boundary_name)
        {
          const auto wall_id = subproblem.mesh().getBoundaryID(name);
          for (const auto side_id : side_bnds)
            if (side_id == wall_id)
            {
              const FaceInfo * const fi = subproblem.mesh().faceInfo(elem, i_side);
              const Real dist = std::abs((fi->elemCentroid() - fi->faceCentroid()) * fi->normal());
              dist_map[elem].push_back(dist);
            }
        }
      }
}

/// Face arguments to wall-bounded faces for wall tretement
void
getElementFaceArgs(const std::vector<BoundaryName> & wall_boundary_name,
                   const FEProblemBase & fe_problem,
                   const SubProblem & subproblem,
                   const std::set<SubdomainID> & block_ids,
                   std::map<const Elem *, std::vector<const FaceInfo *>> & face_info_map)
{

  face_info_map.clear();

  for (const auto & elem : fe_problem.mesh().getMesh().active_element_ptr_range())
    if (block_ids.find(elem->subdomain_id()) != block_ids.end())
      for (const auto i_side : elem->side_index_range())
      {
        const auto & side_bnds = subproblem.mesh().getBoundaryIDs(elem, i_side);
        for (const auto & name : wall_boundary_name)
        {
          const auto wall_id = subproblem.mesh().getBoundaryID(name);
          for (const auto side_id : side_bnds)
            if (side_id == wall_id)
            {
              const FaceInfo * fi = subproblem.mesh().faceInfo(elem, i_side);
              face_info_map[elem].push_back(fi);
            }
        }
      }
}
}
