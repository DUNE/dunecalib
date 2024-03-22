////////////////////////////////////////////////////////////////////////
//// \file RunConditionsProtoDUNE_service.cc
////
//// \brief implementation of class for storing/accessing run conditions parameters for ProtoDUNE
////
//// \author avizcaya@colostate.edu
//// \date Jan 19, 2024
//// 
//////////////////////////////////////////////////////////////////////////


// C++ language includes
#include <iostream>


#include "TTimeStamp.h"

// LArSoft includes
#include "dunecalib/ConIntServices/RunConditionsServicePDUNE.h"

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"


//-----------------------------------------------
condb::RunConditionsServicePDUNE::RunConditionsServicePDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry &reg) {
  std::cout << "Service funtion 1" << std::endl;
  fProp.reset(new condb::RunConditionsProtoDUNE(pset));
  
  fProp->Configure(pset);
  fProp->GetRunNumber();
  fProp->UpdateRN(23300); //(run.id().run());
  fProp->LoadConditionsT();
  condb::RunCond_t rc = fProp->GetRunConditions(0);
  std::cout << "\tstart time = " << rc.start_time << std::endl;
  reg.sPreBeginRun.watch(this, &RunConditionsServicePDUNE::preBeginRun);
}

//-----------------------------------------------
void condb::RunConditionsServicePDUNE::preBeginRun(const art::Run& run) {
  art::Timestamp ts = run.beginTime();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  uint64_t  runtime = tts.AsDouble();
  
  std::cout << "db: runtime in runs conditions Creo que si es este" << runtime << std::endl;
  std::cout << "El numero de run es: " << run.id() << std::endl;

  fProp->UpdateRN(23300); //(run.id().run());
  fProp->LoadConditionsT();
  condb::RunCond_t rc = fProp->GetRunConditions(0);
  std::cout << "\tstart time = " << rc.start_time << std::endl;  
}

//------------------------------------------------
void condb::RunConditionsServicePDUNE::reconfigure(fhicl::ParameterSet const& pset) {
  fProp->Configure(pset);  
  return;
}

DEFINE_ART_SERVICE_INTERFACE_IMPL(condb::RunConditionsServicePDUNE, condb::ConditionsService)
