//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADDirichletBCBaseTempl.h"

/**
 * Dirichlet boundary condition with functor inputs.
 */
class BlochDirichletBC : public ADDirichletBCBaseTempl<Real>
{
public:
  static InputParameters validParams();

  BlochDirichletBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpValue() override;

  const Real & _lattice_vec;
  const Real & _wave_num;
  /// The functor value to impose
  const Moose::Functor<ADReal> & _functor;
  /// Coefficient
  const Moose::Functor<ADReal> & _coef;
};