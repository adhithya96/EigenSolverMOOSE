#include "eigen_xfemApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
eigen_xfemApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

eigen_xfemApp::eigen_xfemApp(InputParameters parameters) : MooseApp(parameters)
{
  eigen_xfemApp::registerAll(_factory, _action_factory, _syntax);
}

eigen_xfemApp::~eigen_xfemApp() {}

void 
eigen_xfemApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<eigen_xfemApp>(f, af, s);
  Registry::registerObjectsTo(f, {"eigen_xfemApp"});
  Registry::registerActionsTo(af, {"eigen_xfemApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
eigen_xfemApp::registerApps()
{
  registerApp(eigen_xfemApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
eigen_xfemApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  eigen_xfemApp::registerAll(f, af, s);
}
extern "C" void
eigen_xfemApp__registerApps()
{
  eigen_xfemApp::registerApps();
}
