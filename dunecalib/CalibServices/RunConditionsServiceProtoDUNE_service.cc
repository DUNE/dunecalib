////////////////////////////////////////////////////////////////////////
// \file RunConditionsProtoDUNE_service.cc
//
// \brief implementation of class for storing/accessing run conditions parameters for ProtoDUNE
//
// \author avizcaya@colostate.edu
// \date Jan 19, 2024
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>

#include "TTimeStamp.h"

// LArSoft includes
#include "dunecalib/CalibServices/RunConditionsServiceProtoDUNE.h"

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

//-----------------------------------------------
runc::RunConditionsServiceProtoDUNE::RunConditionsServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg)
{
  std::cout << "Service funtion 1" << std::endl;
  fProp.reset(new runc::RunConditionsProtoDUNE(pset));
  
  fProp->Configure(pset);
  fProp->GetRunNumber();
  fProp->Update(23300); //(run.id().run());
  fProp->LoadConditionsT();
  runc::RunCond_t rc = fProp->GetRunConditions(0);
  std::cout << "\tstart time = " << rc.start_time << std::endl;
  reg.sPreBeginRun.watch(this, &RunConditionsServiceProtoDUNE::preBeginRun);
}

//----------------------------------------------
void runc::RunConditionsServiceProtoDUNE::preBeginRun(const art::Run& run)
{
  art::Timestamp ts = run.beginTime();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  uint64_t  runtime = tts.AsDouble();
  
  std::cout << "db: runtime in runs conditions " << runtime << std::endl;
  std::cout << "El numero de run es: " << run.id() << std::endl;
  // one can also consider using event time through "sPreProcessEvent" by define a "preProcessEvent" function, for exmaple.
  
  //fProp->Update(runtime);
  fProp->Update(23300); //(run.id().run());
  fProp->LoadConditionsT();
  runc::RunCond_t rc = fProp->GetRunConditions(0);
  std::cout << "\tstart time = " << rc.start_time << std::endl;  
}

//------------------------------------------------
void runc::RunConditionsServiceProtoDUNE::reconfigure(fhicl::ParameterSet const& pset)
{
  fProp->Configure(pset);  
  return;
}

//------------------------------------------------
DEFINE_ART_SERVICE_INTERFACE_IMPL(runc::RunConditionsServiceProtoDUNE, runc::RunConditionsService)