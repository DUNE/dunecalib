////////////////////////////////////////////////////////////////////////
// \file RunConditionsProtoDUNE.cxx
//
// \brief implementation of class for accessing Run Conditions 
//        constants for ProtoDUNE
//
// \author avizcaya@colostate.edu
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

// LArSoft includes
#include "dunecalib/Calib/RunConditionsProtoDUNE.h"

// nutools includes
#include "nuevdb/IFDatabase/Table.h"

// Framework includes
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//-----------------------------------------------
runc::RunConditionsProtoDUNE::RunConditionsProtoDUNE()
{
  fUseCondb = true;
  fIsMC = false;
  fRunConditionsLoaded = false;
  fCurrentTS = 0; // I think this one has to be the run number
  fCSVFileName="";
  fDBTag="";
  fTableName="";
}


//-----------------------------------------------
runc::RunConditionsProtoDUNE::RunConditionsProtoDUNE(
  fhicl::ParameterSet const& pset
)
{
  fUseCondb = true;
  fIsMC = false;
  fRunConditionsLoaded = false;
  fCurrentTS = 0; //Run number?
  fCSVFileName="";
  fDBTag="";
  fTableName="";
}

//------------------------------------------------
bool runc::RunConditionsProtoDUNE::Configure(fhicl::ParameterSet const& pset)
{  
  fUseCondb      = pset.get<bool>("UseCondb");
  fCSVFileName   = pset.get<std::string>("CSVFileName");
  fDBTag         = pset.get<std::string>("DBTag");
  fTableName     = pset.get<std::string>("TableName");
  return true;
}

//------------------------------------------------
bool runc::RunConditionsProtoDUNE::Update(uint64_t ts) 
{

  if (fRunConditionsLoaded && ts != fCurrentTS) {
    fRunCond.clear();
    fRunConditionsLoaded = false;
  }

  fCurrentTS = ts;
  // all done! 

  return true;
}

//------------------------------------------------
runc::RunCond_t runc::RunConditionsProtoDUNE::GetRunConditions(int chanId) 
{
  if (!fRunConditionsLoaded) this->LoadRunConditions();

  if (fRunCond.find(chanId) == fRunCond.end()) {
    mf::LogError("RunConditionsProtoDUNE") << "Channel " << chanId 
				      << "not found!";
    std::abort();
  }

  return fRunCond[chanId];
}

//------------------------------------------------
bool runc::RunConditionsProtoDUNE::LoadRunConditions()
{
  if (!fUseCondb) return true;

  if (fRunConditionsLoaded) return true;

  nutools::dbi::Table t;
  
  t.SetDetector("pdunesp");
  t.SetTableName(fTableName);
  t.SetTableType(nutools::dbi::kConditionsTable);
  t.SetDataTypeMask(nutools::dbi::kDataOnly);
  if (fIsMC)
    t.SetDataTypeMask(nutools::dbi::kMCOnly);
  
  //int run_number = t.AddCol("run","float");
  //int data_type  = t.AddCol("data_type","string");
  //int offsetIdx = t.AddCol("offset","float");
  //int shapeIdx  = t.AddCol("shape","float");
  //int chi2Idx   = t.AddCol("chi2","float");
  //int adcLowIdx = t.AddCol("adc_low","int");
  //int adcHiIdx  = t.AddCol("adc_high","int");
  //int nlIdx[20];
  //char buff[64];
  //for (int i=0; i<20; ++i) {
  //  sprintf(buff,"nl%d",i);
  //  nlIdx[i] = t.AddCol(buff,"float");
  //}
  
  //So as not to interpolate
  t.SetMinTSVld(fCurrentTS);
  t.SetMaxTSVld(fCurrentTS);
  t.SetTag(fDBTag);

  t.SetVerbosity(100);

  bool readOk = false;
  if (!fCSVFileName.empty()) 
    readOk = t.LoadFromCSV(fCSVFileName);
  else
    readOk = t.Load();

  if (! readOk) {
    mf::LogError("RunConditionsProtoDUNE") << "Load from run conditions database table failed.";
    
    return false; //std::abort();

  }
  
  if (t.NRow() == 0) {
    mf::LogError("RunConditionsProtoDUNE") << "Number of rows in run conditions table is 0.  This should never be the case!";
    return false;
  }
  
  nutools::dbi::Row* row;
  uint64_t chan;
  for (int i=0; i<t.NRow(); ++i) {
    RunCond_t c;
    row = t.GetRow(i);
    std::cout << row << std::endl;      
    chan = row->Channel();
    //row->Col(statusIdx).Get(c.status);
    //row->Col(gainIdx).Get(c.gain);
    //row->Col(offsetIdx).Get(c.offset);
    //row->Col(shapeIdx).Get(c.shape);
    //row->Col(chi2Idx).Get(c.chi2);
    //row->Col(adcLowIdx).Get(c.adc_low);
    //row->Col(adcHiIdx).Get(c.adc_high);
    //for (int j=0; j<20; ++j) {
    //  row->Col(nlIdx[j]).Get(c.nl[j]);
    //}
    
    fRunCond[chan] = c;
  }    

  fRunConditionsLoaded = true;
  return true;
}


