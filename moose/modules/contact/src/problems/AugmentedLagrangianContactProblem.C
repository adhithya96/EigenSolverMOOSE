//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AugmentedLagrangianContactProblem.h"

// MOOSE includes
#include "AuxiliarySystem.h"
#include "DisplacedProblem.h"
#include "MooseApp.h"
#include "MooseMesh.h"
#include "MooseVariable.h"
#include "NearestNodeLocator.h"
#include "NonlinearSystem.h"
#include "PenetrationLocator.h"

#include "SystemBase.h"
#include "Assembly.h"
#include "Executioner.h"
#include "AddVariableAction.h"
#include "ConstraintWarehouse.h"
#include "MortarUserObject.h"
#include "AugmentedLagrangeInterface.h"

registerMooseObject("ContactApp", AugmentedLagrangianContactProblem);
registerMooseObject("ContactApp", AugmentedLagrangianContactFEProblem);

InputParameters
AugmentedLagrangianContactProblemInterface::validParams()
{
  auto params = emptyInputParameters();
  params.addParam<int>("maximum_lagrangian_update_iterations",
                       100,
                       "Maximum number of update Lagrangian Multiplier iterations per step");
  return params;
}

AugmentedLagrangianContactProblemInterface::AugmentedLagrangianContactProblemInterface(
    const InputParameters & params)
  : _maximum_number_lagrangian_iterations(params.get<int>("maximum_lagrangian_update_iterations"))
{
}

template <class T>
InputParameters
AugmentedLagrangianContactProblemTempl<T>::validParams()
{
  InputParameters params = T::validParams();
  params += AugmentedLagrangianContactProblemInterface::validParams();
  params.addClassDescription("Manages nested solution for augmented Lagrange contact");
  return params;
}

template <class T>
AugmentedLagrangianContactProblemTempl<T>::AugmentedLagrangianContactProblemTempl(
    const InputParameters & params)
  : T(params), AugmentedLagrangianContactProblemInterface(params)
{
}

template <class T>
void
AugmentedLagrangianContactProblemTempl<T>::timestepSetup()
{
  _lagrangian_iteration_number = 0;
  T::timestepSetup();
}

template <class T>
MooseNonlinearConvergenceReason
AugmentedLagrangianContactProblemTempl<T>::checkNonlinearConvergence(std::string & msg,
                                                                     const PetscInt it,
                                                                     const Real xnorm,
                                                                     const Real snorm,
                                                                     const Real fnorm,
                                                                     const Real rtol,
                                                                     const Real divtol,
                                                                     const Real stol,
                                                                     const Real abstol,
                                                                     const PetscInt nfuncs,
                                                                     const PetscInt /*max_funcs*/,
                                                                     const Real ref_resid,
                                                                     const Real /*div_threshold*/)
{
  Real my_max_funcs = std::numeric_limits<int>::max();
  Real my_div_threshold = std::numeric_limits<Real>::max();

  MooseNonlinearConvergenceReason reason = T::checkNonlinearConvergence(msg,
                                                                        it,
                                                                        xnorm,
                                                                        snorm,
                                                                        fnorm,
                                                                        rtol,
                                                                        divtol,
                                                                        stol,
                                                                        abstol,
                                                                        nfuncs,
                                                                        my_max_funcs,
                                                                        ref_resid,
                                                                        my_div_threshold);

  _console << "Augmented Lagrangian contact iteration " << _lagrangian_iteration_number
           << std::endl;

  bool repeat_augmented_lagrange_step = false;

  if (reason == MooseNonlinearConvergenceReason::CONVERGED_FNORM_ABS ||
      reason == MooseNonlinearConvergenceReason::CONVERGED_FNORM_RELATIVE ||
      reason == MooseNonlinearConvergenceReason::CONVERGED_SNORM_RELATIVE)
  {
    if (_lagrangian_iteration_number < _maximum_number_lagrangian_iterations)
    {
      auto & nonlinear_sys = currentNonlinearSystem();
      nonlinear_sys.update();

      // Get the penetration locator from the displaced mesh if it exist, otherwise get
      // it from the undisplaced mesh.
      const auto displaced_problem = getDisplacedProblem();
      const auto & penetration_locators =
          (displaced_problem ? displaced_problem->geomSearchData() : geomSearchData())
              ._penetration_locators;

      // loop over contact pairs (penetration locators)
      const ConstraintWarehouse & constraints = nonlinear_sys.getConstraintWarehouse();
      std::list<std::shared_ptr<MechanicalContactConstraint>> mccs;
      for (const auto & pair : penetration_locators)
      {
        const auto & boundaries = pair.first;

        if (!constraints.hasActiveNodeFaceConstraints(boundaries.second, bool(displaced_problem)))
          continue;
        const auto & ncs =
            constraints.getActiveNodeFaceConstraints(boundaries.second, bool(displaced_problem));

        mccs.emplace_back(nullptr);
        for (const auto & nc : ncs)
          if (const auto mcc = std::dynamic_pointer_cast<MechanicalContactConstraint>(nc); !mcc)
            mooseError("AugmentedLagrangianContactProblem: dynamic cast of "
                       "MechanicalContactConstraint object failed.");
          else
          {
            // Return if this constraint does not correspond to the primary-secondary pair
            // prepared by the outer loops.
            // This continue statement is required when, e.g. one secondary surface constrains
            // more than one primary surface.
            if (mcc->secondaryBoundary() != boundaries.second ||
                mcc->primaryBoundary() != boundaries.first)
              continue;

            // save one constraint pointer for each contact pair
            if (!mccs.back())
              mccs.back() = mcc;

            // check if any of the constraints is not yet converged
            if (repeat_augmented_lagrange_step || !mcc->AugmentedLagrangianContactConverged())
            {
              repeat_augmented_lagrange_step = true;
              break;
            }
          }
      }

      // next loop over penalty mortar user objects
      const auto & pmuos = this->_app.template getInterfaceObjects<AugmentedLagrangeInterface>();
      for (auto * pmuo : pmuos)
      {
        // check if any of the constraints is not yet converged
        if (!repeat_augmented_lagrange_step && !pmuo->isAugmentedLagrangianConverged())
          repeat_augmented_lagrange_step = true;
      }

      // Communicate the repeat_augmented_lagrange_step in parallel.
      // If one proc needs to do another loop, all do.
      this->_communicator.max(repeat_augmented_lagrange_step);

      // repeat update step if necessary
      if (repeat_augmented_lagrange_step)
      {
        _lagrangian_iteration_number++;

        // Each contact pair will have constraints for all displacements, but those share the
        // Lagrange multipliers, which are stored on the penetration locator. We call update
        // only for the first constraint for each contact pair.
        for (const auto & mcc : mccs)
          mcc->updateAugmentedLagrangianMultiplier(/* beginning_of_step = */ false);

        // Update all penalty mortar user objects
        for (const auto & pmuo : pmuos)
          pmuo->updateAugmentedLagrangianMultipliers();

        // call AM setup again (e.g. to update active sets)
        for (const auto & pmuo : pmuos)
          pmuo->augmentedLagrangianSetup();

        // force it to keep iterating
        reason = MooseNonlinearConvergenceReason::ITERATING;
        _console << "Augmented Lagrangian Multiplier needs updating." << std::endl;
      }
      else
        _console << "Augmented Lagrangian contact constraint enforcement is satisfied."
                 << std::endl;
    }
    else
    {
      // maxed out
      _console << "Maximum Augmented Lagrangian contact iterations have been reached." << std::endl;
      reason = MooseNonlinearConvergenceReason::DIVERGED_FUNCTION_COUNT;
    }
  }

  return reason;
}

template class AugmentedLagrangianContactProblemTempl<ReferenceResidualProblem>;
template class AugmentedLagrangianContactProblemTempl<FEProblem>;
