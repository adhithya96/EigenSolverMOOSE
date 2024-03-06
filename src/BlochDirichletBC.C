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
  
  params.addParam<double>("lattice_length", 1.0, "length of the phononic lattice");
  params.declareControllable("lattice_length");
  params.addParam<double>("wave_number", 1.0, "wave number exp(ikx)");
  params.declareControllable("wave_number");

  return params;
}

BlochDirichletBC::BlochDirichletBC(const InputParameters & parameters)
  : ADDirichletBCBaseTempl<Real>(parameters),
    _lattice_vec(getParam<double>("lattice_length")),
    _wave_num(getParam<double>("wave_number")),
    _functor(getFunctor<ADReal>("functor")),
    _coef(getFunctor<ADReal>("coefficient"))
{
}

// u(x + h) = u(x) * exp(ikx)
//Translation vector  is used  to  define  the  left side of  the  boundary 
ADReal
BlochDirichletBC::computeQpValue()
{
  //Get point  locator  
  const auto pl = _mesh.getPointLocator();
  //Locate node on the other side of the geometry
  //The translation vector will be a user parameter
  const libMesh::Point current_point = pl->Point();
  const auto new_node = pl->locate_node(current_point);

  if(!new_node)
    mooseError("Did not find the opposite side value");
  
  const Moose::NodeArg space_arg = {new_node, Moose::INVALID_BLOCK_ID};
  
  const Moose::StateArg time_arg = Moose::currentState();

  
  return _coef(space_arg, time_arg) * _functor(space_arg, time_arg) * exp(_wave_num * _lattice_vec);

}