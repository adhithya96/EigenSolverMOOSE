//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BlochDirichletBCReal.h"

registerMooseObject("MooseApp", BlochDirichletBCReal);

InputParameters
BlochDirichletBCReal::validParams()
{
  InputParameters params = ADDirichletBCBase::validParams();

  params.addClassDescription("Imposes the Bloch boundary condition on displacements "
                             "$u(t,\\vec{x + h})=u(t,\\vec{x})exp(ikx)$, ");

  params.addRequiredParam<MooseFunctorName>("uim", "Functor used for uim");
  params.addRequiredParam<MooseFunctorName>("ur", "Functor used for ur");
  
  params.addParam<MooseFunctorName>(
      "coefficient", 1.0, "An optional functor coefficient to multiply the imposed functor");
  
  params.addParam<double>("lattice_length", 1.0, "length of the phononic lattice");
  params.declareControllable("lattice_length");
  params.addParam<double>("wave_number", 1.0, "wave number exp(ikx)");
  params.declareControllable("wave_number");

  return params;
}

BlochDirichletBCReal::BlochDirichletBCReal(const InputParameters & parameters)
  : ADDirichletBCBaseTempl<Real>(parameters),
    _lattice_vec(getParam<double>("lattice_length")),
    _wave_num(getParam<double>("wave_number")),
    _uim(getFunctor<ADReal>("uim")),
    _ur(getFunctor<ADReal>("ur")),
    _coef(getFunctor<ADReal>("coefficient"))
{
}

// u(x + h) = u(x) * exp(ikx)
//Translation vector  is used  to  define  the right side of  the boundary
ADReal
BlochDirichletBCReal::computeQpValue()
{
  //Get point  locator  
  const auto pl = _mesh.getPointLocator();
  //Locate node on the other side of the geometry
  //The translation vector will be a user parameter
  libMesh::Point _translation_vec= libMesh::Point(_lattice_vec, 0.0, 0.0);
  const auto new_node = pl->locate_node(*_current_node - _translation_vec);


  if(!new_node)
    mooseError("Did not find the opposite side value");
  
  const Moose::NodeArg space_arg = {new_node, Moose::INVALID_BLOCK_ID};
  
  const Moose::StateArg time_arg = Moose::currentState();

  return _coef(space_arg, time_arg) * (_ur(space_arg, time_arg) * cos (_wave_num * _lattice_vec)  -  _uim(space_arg, time_arg) * sin (_wave_num * _lattice_vec));

}