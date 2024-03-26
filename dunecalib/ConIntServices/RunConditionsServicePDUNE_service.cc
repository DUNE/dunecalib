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

  //Conditions service functions
  fProp.reset(new condb::RunConditionsProtoDUNE(pset));
  fProp->Configure(pset);
  float rn = fProp->GetRunNumber();
  fProp->UpdateRN(rn); 
  fProp->LoadConditionsT();
  //Load table parameters in rc
  condb::RunCond_t rc = fProp->GetRunConditions(0);
  std::cout << "\tstart time = " << rc.start_time 
            << "\n\tstop time  = " << rc.stop_time
            << "\n\tdata type = " << rc.data_type << std::endl;
  reg.sPreBeginRun.watch(this, &RunConditionsServicePDUNE::preBeginRun);
}

//-----------------------------------------------
void condb::RunConditionsServicePDUNE::preBeginRun(const art::Run& run) {
  art::Timestamp ts = run.beginTime();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  uint64_t  runtime = tts.AsDouble();
  
  std::cout << "db: runtime " << runtime << std::endl;

}

//------------------------------------------------
void condb::RunConditionsServicePDUNE::reconfigure(fhicl::ParameterSet const& pset) {
  fProp->Configure(pset);  
  return;
}

DEFINE_ART_SERVICE_INTERFACE_IMPL(condb::RunConditionsServicePDUNE, condb::ConditionsService)
