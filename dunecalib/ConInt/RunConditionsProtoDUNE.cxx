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
#include <optional>

// LArSoft includes
#include "dunecalib/ConInt/RunConditionsProtoDUNE.h"
#include "wda.h"

// nutools includes
#include "nuevdb/IFDatabase/Table.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


//------------------------------------------------
// Configure the parameters with the fhicl values
bool condb::RunConditionsProtoDUNE::Configure(fhicl::ParameterSet const& pset) {
  fDBTag         = pset.get<std::string>("DBTag");
  fTableURL      = pset.get<std::string>("TableURL");
  fTableName     = pset.get<std::string>("TableName");
  fRunNumber     = pset.get<float>("RunNumber");
  fRunNumber1    = pset.get<float>("RunNumber1");
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
// Get Run Conditions values to Null or reset to some value that is clear 
// its not the set value but a default
condb::RunCond_t condb::ResetRunCond_t(condb::RunCond_t rct){
  std::optional<float> myFloat = std::nullopt;
  if (myFloat) {std::cout << "myFloat contains: " << *myFloat << std::endl;}
  rct ={
  -99999, //run_number
  "None", //data_type
  -99999, //upload_t
  -99999, //start_time
  -99999, //stop_time
  "None", //run_type
  "None", //detector_id;
  "None", //software_version
  -99999, //buffer
  "None", //ac_couple 
  -99999, //baseline
  false, //enabled
  -99999, //gain
  false, //gain_match
  -99999, //leak
  false, //leak_10x
  -99999, //leak_f
  -99999, //peak_time
  -99999, //pulse_dac
  -99999, //strobe_delay
  -99999, //strobe_length
  -99999, //strobe_skip
  false, //test_cap
  false, //adc_test_pattern
  false, //cold
  "None", //detector_type
  false, //pulser
  -99999, //beam_momentum
  -99999, //beam_polarity
  -99999, //detector_hv
  -99999, //detector_hvset
  -99999 //beam_setmomentum
  };
  return rct;
}


//------------------------------------------------
// Get the run conditions for selected channel
condb::RunCond_t condb::RunConditionsProtoDUNE::GetRunConditions(float chanId) { //(int chanId) {
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
  //Add table and define parameters
  Table ct;
  //URL to where the database is
  ct.SetFolderURL(fTableURL);
  // table name of the database with the format schema.table or if no schema just the table name
  ct.SetFolderName(fTableName);
  //How much feedback do you want. 0 is none, 2 is all
  ct.SetVerbosity(fVerbosity); 
  ct.GetRangeOfValues(fCurrentRN,fRunNumber1);
  ct.SetTTag(fDBTag);
  
 
  //Add the column names and types of parameters that you want
  //example
  //int data_typeIdx  = t.AddCol("data_type","string");
  int data_typeIdx  = ct.AddCol("data_type","string");
  int upload_tIdx   = ct.AddCol("upload_time","float");
  int start_timeIdx = ct.AddCol("start_time","float");
  int stop_timeIdx  = ct.AddCol("stop_time","float");
  int software_versionIdx  = ct.AddCol("software_version","string");
  int run_typeIdx   = ct.AddCol("run_type","string");
  int bufferIdx     = ct.AddCol("buffer","float");
  int ac_coupleIdx  = ct.AddCol("ac_couple","bool");
  int detector_idIdx   = ct.AddCol("detector_id","string");
  int baselineIdx   = ct.AddCol("baseline", "integer");
  int enabledIdx    = ct.AddCol("enabled", " bool");
  int gainIdx       = ct.AddCol("gain", "float");
  int gain_matchIdx = ct.AddCol("gain_match", "bool");
  int leakIdx       = ct.AddCol("leak",  "float");
  int leak_10xIdx   = ct.AddCol("leak_10x",  "bool");
  int leak_fIdx     = ct.AddCol("leak_f", "float");
  int peak_timeIdx  = ct.AddCol("peak_time", "float");
  int pulse_dacIdx  = ct.AddCol("pulse_dac", "integer");
  int strobe_delayIdx = ct.AddCol("strobe_delay", "integer");
  int strobe_lengthIdx= ct.AddCol("strobe_length", "integer");
  int strobe_skipIdx = ct.AddCol("strobe_skip", "integer");
  int test_capIdx    = ct.AddCol("test_cap", "bool");
  int adc_test_patternIdx = ct.AddCol("adc_test_pattern", "bool");
  int coldIdx       = ct.AddCol("cold", "bool");
  int detector_typeIdx = ct.AddCol("detector_type", "string");
  int pulserIdx     = ct.AddCol("pulser", "bool");
  int beam_momentumIdx = ct.AddCol("beam_momentum", "float");
  int beam_polarityIdx = ct.AddCol("beam_polarity", "float");
  int detector_hvIdx = ct.AddCol("detector_hv", "float");
  int detector_hvsetIdx = ct.AddCol("detector_hvset", "float");
  int beam_setmomentumIdx = ct.AddCol("beam_setmomentum", "float");

  ct.LoadConditionsTable();
  if (ct.NRow() == 0) {
    mf::LogError("RunConditionsProtoDUNE") << "Number of rows in run conditions table is 0.  This should never be the case!";
    return false;
  }
  //long int chan;
  nutools::dbi::Row* row;
  float run_row = fCurrentRN; 
  for (int i=0; i<ct.NRow()-1; ++i) {
    RunCond_t c; 
    c = ResetRunCond_t(c); 
    row = ct.GetRow(i);
    //chan = row->Channel();
    c.run_number = row->VldTime(); //for tables with run as key

    //For tables with run as key, if c.run_number is not the same as the run number then don't interpolate
    //and return that the run info is not on the db   
    if (c.run_number != run_row  ) {
      mf::LogError("RunConditionsProtoDUNE") << "Run number " << run_row << " is not on the database!";
      return false;
    }
    run_row = run_row + 1.0;
    
    // Fill the columns of each row with data
    row->Col(data_typeIdx).Get(c.data_type);
    row->Col(upload_tIdx).Get(c.upload_t);
    row->Col(start_timeIdx).Get(c.start_time);
    row->Col(software_versionIdx).Get(c.software_version);
    row->Col(stop_timeIdx).Get(c.stop_time);    
    row->Col(run_typeIdx).Get(c.run_type);
    row->Col(bufferIdx).Get(c.buffer);
    row->Col(ac_coupleIdx).Get(c.ac_couple);
    row->Col(detector_idIdx).Get(c.detector_id);
    row->Col(baselineIdx).Get(c.baseline);
    row->Col(enabledIdx).Get(c.enabled);
    row->Col(gainIdx).Get(c.gain);
    row->Col(gain_matchIdx).Get(c.gain_match);
    row->Col(leakIdx).Get(c.leak);
    row->Col(leak_10xIdx).Get(c.leak_10x);
    row->Col(leak_fIdx).Get(c.leak_f);
    row->Col(peak_timeIdx).Get(c.peak_time);
    row->Col(pulse_dacIdx).Get(c.pulse_dac);
    row->Col(strobe_delayIdx).Get(c.strobe_delay);
    row->Col(strobe_lengthIdx).Get(c.strobe_length);
    row->Col(strobe_skipIdx).Get(c.strobe_skip);
    row->Col(test_capIdx).Get(c.test_cap);
    row->Col(adc_test_patternIdx).Get(c.adc_test_pattern);
    row->Col(coldIdx).Get(c.cold);
    row->Col(detector_typeIdx).Get(c.detector_type);
    row->Col(pulserIdx).Get(c.pulser);
    row->Col(beam_momentumIdx).Get(c.beam_momentum);
    row->Col(beam_polarityIdx).Get(c.beam_polarity);
    row->Col(detector_hvIdx).Get(c.detector_hv);
    row->Col(detector_hvsetIdx).Get(c.detector_hvset);
    row->Col(beam_setmomentumIdx).Get(c.beam_setmomentum);

    //fRunCond[chan] = c;
    fRunCond[c.run_number] = c;
  } 
  fRunConditionsLoaded = true;
  return true; 
}


 
