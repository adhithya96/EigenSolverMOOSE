//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlochDirichletBC.h"

registerMooseObject("MooseApp", BlochDirichletBC);

InputParameters
BlochDirichletBC::validParams()
{
  InputParameters params = ADDirichletBCBase::validParams();

  params.addClassDescription("Imposes the Dirichlet boundary condition "
                             "$u(t,\\vec{x})=h(t,\\vec{x})$, "
                             "where $h$ is a functor and can have complex dependencies.");

  params.addRequiredParam<MooseFunctorName>("functor", "The functor to impose");
  params.addParam<MooseFunctorName>(
      "coefficient", 1.0, "An optional functor coefficient to multiply the imposed functor");
  
  params.addParam<Real>("lattice length");
  params.declareControllable("lattice length");
  params.addParam<Real>("wave number");
  params.declareControllable("wave number");
  



  return params;
}

BlochDirichletBC::BlochDirichletBC(const InputParameters & parameters)
  : ADDirichletBCBaseTempl<Real>(parameters),
    _lattice_vec(getParam<Real>("lattice length")),
    _wave_num(getParam<Real>("wave number")),
    _functor(getFunctor<ADReal>("functor")),
    _coef(getFunctor<ADReal>("coefficient"))
{
}

ADReal
BlochDirichletBC::computeQpValue()
{
  //Get point  locator  
  const auto pl = _mesh.getPointLocator();
  //Locate node on the other side of the geometry
  //The translation vector will be a user parameter
  const auto new_node = pl.locate_node(_current_node + _translation_vector);

  if(!new_node)
    mooseError("Did not find the opposite side value");
  
  const Moose::NodeArg space_arg = {new_node, Moose::INVALID_BLOCK_ID};
  
  const Moose::StateArg time_arg = Moose::currentState();

  
  return _coef(space_arg, time_arg) * _functor(space_arg, time_arg) *  _u * exp(_wave_num * _lattice_vec);

}