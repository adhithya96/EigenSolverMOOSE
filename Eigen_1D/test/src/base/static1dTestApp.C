//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "static1dTestApp.h"
#include "static1dApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
static1dTestApp::validParams()
{
  InputParameters params = static1dApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

static1dTestApp::static1dTestApp(InputParameters parameters) : MooseApp(parameters)
{
  static1dTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

static1dTestApp::~static1dTestApp() {}

void
static1dTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  static1dApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"static1dTestApp"});
    Registry::registerActionsTo(af, {"static1dTestApp"});
  }
}

void
static1dTestApp::registerApps()
{
  registerApp(static1dApp);
  registerApp(static1dTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
static1dTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  static1dTestApp::registerAll(f, af, s);
}
extern "C" void
static1dTestApp__registerApps()
{
  static1dTestApp::registerApps();
}
