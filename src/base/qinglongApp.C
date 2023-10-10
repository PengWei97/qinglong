#include "qinglongApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
qinglongApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

qinglongApp::qinglongApp(InputParameters parameters) : MooseApp(parameters)
{
  qinglongApp::registerAll(_factory, _action_factory, _syntax);
}

qinglongApp::~qinglongApp() {}

void 
qinglongApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<qinglongApp>(f, af, s);
  Registry::registerObjectsTo(f, {"qinglongApp"});
  Registry::registerActionsTo(af, {"qinglongApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
qinglongApp::registerApps()
{
  registerApp(qinglongApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
qinglongApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  qinglongApp::registerAll(f, af, s);
}
extern "C" void
qinglongApp__registerApps()
{
  qinglongApp::registerApps();
}
