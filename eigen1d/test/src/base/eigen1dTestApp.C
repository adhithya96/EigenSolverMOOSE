//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "eigen1dTestApp.h"
#include "eigen1dApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
eigen1dTestApp::validParams()
{
  InputParameters params = eigen1dApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

eigen1dTestApp::eigen1dTestApp(InputParameters parameters) : MooseApp(parameters)
{
  eigen1dTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

eigen1dTestApp::~eigen1dTestApp() {}

void
eigen1dTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  eigen1dApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"eigen1dTestApp"});
    Registry::registerActionsTo(af, {"eigen1dTestApp"});
  }
}

void
eigen1dTestApp::registerApps()
{
  registerApp(eigen1dApp);
  registerApp(eigen1dTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
eigen1dTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  eigen1dTestApp::registerAll(f, af, s);
}
extern "C" void
eigen1dTestApp__registerApps()
{
  eigen1dTestApp::registerApps();
}
