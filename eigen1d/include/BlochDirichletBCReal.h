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
class BlochDirichletBCReal : public ADDirichletBCBaseTempl<Real>
{
public:
  static InputParameters validParams();

  BlochDirichletBCReal(const InputParameters & parameters);

protected:
  virtual ADReal computeQpValue() override;

  double _lattice_vec;
  double _wave_num;
  /// The functor value to impose
  const Moose::Functor<ADReal> & _uim;
  const Moose::Functor<ADReal> & _ur;
  /// Coefficient
  const Moose::Functor<ADReal> & _coef;
};
