////////////////////////////////////////////////////////////////////////
//// \file RunConditionsProtoDUNE.h
////
//// \brief header of class for accessing run conditions 
////        data for ProtoDUNE
////
//// \author avizcaya@colostate.edu
//// 
//////////////////////////////////////////////////////////////////////////
#ifndef CONDB_RUNCONDITIONSPROTODUNE_H
#define CONDB_RUNCONDITIONSPROTODUNE_H

// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT includes

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>

// dunecalib includes
#include "dunecalib/ConInt/Condb.h"

namespace condb {

  // All the parameters from the table in the conditions DB that you want to include
  // with the format: type name;
  // example
  // float run_number;
  typedef struct {
    float run_number;
    std::string data_type;
    float upload_t;
    float start_time;
    float stop_time;
    std::string run_type;
    std::string software_version;
    int buffer;
    bool ac_couple;
  } RunCond_t; //Name your structure


  class RunConditionsProtoDUNE : public Conditions {

  public:

    RunConditionsProtoDUNE();
    RunConditionsProtoDUNE(fhicl::ParameterSet const& pset);
    RunConditionsProtoDUNE(RunConditionsProtoDUNE const&) = delete;
    virtual ~RunConditionsProtoDUNE() = default;

    //I think most of this functions can be in the inheritance class
    bool Configure(fhicl::ParameterSet const& pset);
    bool UpdateTS(uint64_t ts=0);
    bool UpdateRN(float rn=0);

    bool LoadConditionsT();
    RunCond_t GetRunConditions(int chanId); 

    //void SetTag(std::string ta) { fDBTag = ta; }
    //void SetTableName(std::string tn) {fTableName = tn; }
    //void SetTableURL(std::string url) {fTableURL = url; }
    //void SetRunNumber(float rn) {fRunNumber = rn;}
    //void SetVerbosity(int vr) {fVerbosity = vr;}

    //void GetRunNumber() {std::cout << "Run Number: "<< fRunNumber << std::endl;}

  protected:
    bool LoadRunConditions();
    
    //bool fRunConditionsLoaded;
    //uint64_t fCurrentTS;
    //float fCurrentRN;
    //float fRunNumber;
    //std::string fDBTag;
    //std::string fTableName;
    //std::string fTableURL;
    //int fVerbosity;

    std::map<int,RunCond_t> fRunCond;    
    

  }; //class RunConditionsProtoDUNE
} //namespace condb
#endif // CONDB_RUNCONDITIONSPROTODUNE_H

