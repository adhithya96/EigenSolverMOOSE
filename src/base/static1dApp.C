#include "static1dApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
static1dApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

static1dApp::static1dApp(InputParameters parameters) : MooseApp(parameters)
{
  static1dApp::registerAll(_factory, _action_factory, _syntax);
}

static1dApp::~static1dApp() {}

void 
static1dApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<static1dApp>(f, af, s);
  Registry::registerObjectsTo(f, {"static1dApp"});
  Registry::registerActionsTo(af, {"static1dApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
static1dApp::registerApps()
{
  registerApp(static1dApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
static1dApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  static1dApp::registerAll(f, af, s);
}
extern "C" void
static1dApp__registerApps()
{
  static1dApp::registerApps();
}
