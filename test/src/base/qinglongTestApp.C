//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "qinglongTestApp.h"
#include "qinglongApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
qinglongTestApp::validParams()
{
  InputParameters params = qinglongApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

qinglongTestApp::qinglongTestApp(InputParameters parameters) : MooseApp(parameters)
{
  qinglongTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

qinglongTestApp::~qinglongTestApp() {}

void
qinglongTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  qinglongApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"qinglongTestApp"});
    Registry::registerActionsTo(af, {"qinglongTestApp"});
  }
}

void
qinglongTestApp::registerApps()
{
  registerApp(qinglongApp);
  registerApp(qinglongTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
qinglongTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  qinglongTestApp::registerAll(f, af, s);
}
extern "C" void
qinglongTestApp__registerApps()
{
  qinglongTestApp::registerApps();
}
