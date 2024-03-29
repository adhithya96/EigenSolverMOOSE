//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Reaction.h"

/**
 * Computes a finite element mass matrix meant for use in preconditioning schemes which require one
 */
class MassMatrix : public Reaction
{
public:
  static InputParameters validParams();

  MassMatrix(const InputParameters & parameters);

  virtual void computeResidual() override {}
};
