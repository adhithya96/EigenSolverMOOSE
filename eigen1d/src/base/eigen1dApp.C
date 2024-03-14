#include "eigen1dApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
eigen1dApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

eigen1dApp::eigen1dApp(InputParameters parameters) : MooseApp(parameters)
{
  eigen1dApp::registerAll(_factory, _action_factory, _syntax);
}

eigen1dApp::~eigen1dApp() {}

void 
eigen1dApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<eigen1dApp>(f, af, s);
  Registry::registerObjectsTo(f, {"eigen1dApp"});
  Registry::registerActionsTo(af, {"eigen1dApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
eigen1dApp::registerApps()
{
  registerApp(eigen1dApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
eigen1dApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  eigen1dApp::registerAll(f, af, s);
}
extern "C" void
eigen1dApp__registerApps()
{
  eigen1dApp::registerApps();
}
