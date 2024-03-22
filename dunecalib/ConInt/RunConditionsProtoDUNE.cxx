////////////////////////////////////////////////////////////////////////
//// \file RunConditionsProtoDUNE.cxx
////
//// \brief implementation of class for accessing Run Conditions 
////        constants for ProtoDUNE
////
//// \author avizcaya@colostate.edu
//// 
//////////////////////////////////////////////////////////////////////////


// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

// LArSoft includes
#include "dunecalib/ConInt/RunConditionsProtoDUNE.h"
#include "wda.h"

// nutools includes
#include "nuevdb/IFDatabase/Table.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------
//Initialize values
condb::RunConditionsProtoDUNE::RunConditionsProtoDUNE() {
  fRunConditionsLoaded = false;
  fCurrentTS = 0;
  fRunNumber = 0;
  fCurrentRN = 0;
  fVerbosity = 0;
  fTableName = "";
  fTableURL  = "";
  fDBTag   = "";
}

//-----------------------------------------------
//Initialize values
condb::RunConditionsProtoDUNE::RunConditionsProtoDUNE(fhicl::ParameterSet const& pset) {
  fRunConditionsLoaded = false;
  fCurrentTS = 0;
  fRunNumber = 0;
  fCurrentRN = 0;
  fVerbosity = 0;
  fTableName ="";
  fTableURL = "";
  fDBTag    = "";
}


//------------------------------------------------
// Configure the parameters with the fhicl values
bool condb::RunConditionsProtoDUNE::Configure(fhicl::ParameterSet const& pset) {
  fDBTag         = pset.get<std::string>("DBTag");
  fTableURL      = pset.get<std::string>("TableURL");
  fTableName     = pset.get<std::string>("TableName");
  fRunNumber     = pset.get<float>("RunNumber");
  fVerbosity     = pset.get<int>("Verbosity");
  return true;
}

//------------------------------------------------
// Get new time 
bool condb::RunConditionsProtoDUNE::UpdateTS(uint64_t ts) {
  if (fRunConditionsLoaded && ts != fCurrentTS) {
    fRunCond.clear();
    fRunConditionsLoaded = false;
  }
  fCurrentTS = ts;
  return true;
}

//------------------------------------------------
// Get new run number
bool condb::RunConditionsProtoDUNE::UpdateRN(float rn) {
  if (fRunConditionsLoaded && rn != fCurrentRN) {
    fRunCond.clear();
    fRunConditionsLoaded = false;
  }
  fCurrentRN = rn;
  return true;
}

//------------------------------------------------
// Get the run conditions for selected channel
condb::RunCond_t condb::RunConditionsProtoDUNE::GetRunConditions(int chanId) {
  if (!fRunConditionsLoaded) this->LoadConditionsT();
  if (fRunCond.find(chanId) == fRunCond.end()) {
    mf::LogError("RunConditionsProtoDUNE") << "Channel " << chanId << "not found!";
    std::abort();
  }
  return fRunCond[chanId];
}

//------------------------------------------------
//Load the table info to the variables
bool condb::RunConditionsProtoDUNE::LoadConditionsT() {
  Table ct;
  //URL to where the database is
  ct.SetFolderURL(fTableURL);
  // table name of the database with the format schema.table or if no schema just the table name
  ct.SetFolderName(fTableName);
  //How much feedback do you want. 0 is none, 2 is all
  ct.SetVerbosity(fVerbosity); 

  // So as not to interpolate -> Should update!
  ct.SetMinTSVld(fCurrentRN);
  ct.SetMaxTSVld(fCurrentRN); 

  //Add the column names and types of parameters that you want
  //example
  //int data_typeIdx  = t.AddCol("data_type","string");
  int data_typeIdx  = ct.AddCol("data_type","string");
  int upload_tIdx   = ct.AddCol("upload_time","float");
  int start_timeIdx = ct.AddCol("start_time","float");
  int stop_timeIdx  = ct.AddCol("stop_time","float");
  int software_versionIdx  = ct.AddCol("software_version","string");
  int run_typeIdx   = ct.AddCol("run_type","string");
  int bufferIdx        = ct.AddCol("buffer","float");
  int ac_coupleIdx     = ct.AddCol("ac_couple","bool");

  ct.LoadConditionsTable();
  if (ct.NRow() == 0) {
    mf::LogError("RunConditionsProtoDUNE") << "Number of rows in run conditions table is 0.  This should never be the case!";
    return false;
  }
  long int chan;
  nutools::dbi::Row* row;
  for (int i=0; i<ct.NRow(); ++i) {
    RunCond_t c;
    row = ct.GetRow(i);
    chan = row->Channel();
    c.run_number = row->VldTime();
    row->Col(data_typeIdx).Get(c.data_type);
    row->Col(upload_tIdx).Get(c.upload_t);
    row->Col(start_timeIdx).Get(c.start_time);
    row->Col(software_versionIdx).Get(c.software_version);
    row->Col(stop_timeIdx).Get(c.stop_time);    
    row->Col(run_typeIdx).Get(c.run_type);
    row->Col(bufferIdx).Get(c.buffer);
    row->Col(ac_coupleIdx).Get(c.ac_couple);

    fRunCond[chan] = c;
  } 
  fRunConditionsLoaded = true;
  return true; 
}














  
