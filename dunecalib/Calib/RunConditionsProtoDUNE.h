////////////////////////////////////////////////////////////////////////
// \file RunConditionsProtoDUNE.h
//
// \brief header of class for accessing run conditions 
//        data for ProtoDUNE
//
// \author avizcaya@colostate.edu
// 
////////////////////////////////////////////////////////////////////////
#ifndef RUNC_RUNCONDITIONSPROTODUNE_H
#define RUNC_RUNCONDITIONSPROTODUNE_H


// FHiCL libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT includes
#include "TH1F.h"
#include "TH2F.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>

// dunetpc includes
#include "dunecalib/Calib/RunConditions.h"

namespace runc {

  typedef struct {
    // All the parameters from the run conditinos DB 
    float run_number;
    std::string data_type;
    float upload_t;
    float start_time;
    float stop_time;
    std::string run_type;
    std::string software_version;
    //int buffer;
    //bool ac_couple;
  } RunCond_t;

  class RunConditionsProtoDUNE : public RunConditions {
    
  public:

    RunConditionsProtoDUNE();
    RunConditionsProtoDUNE(fhicl::ParameterSet const& pset);
    RunConditionsProtoDUNE(RunConditionsProtoDUNE const&) = delete;
    virtual ~RunConditionsProtoDUNE() = default;
      
    bool Configure(fhicl::ParameterSet const& pset);
    bool Update(uint64_t ts=0);

    bool LoadConditionsT();

    
    RunCond_t GetRunConditions(int chanId); 
    
    void SetIsMC(bool v) { fIsMC = v; }
    void SetUseCondb(bool v) { fUseCondb = v; }
    void SetCSVFileName(std::string f) { fCSVFileName = f; }
    void SetTag(std::string ta) { fDBTag = ta; }
    void SetTableName(std::string tn) {fTableName = tn; }
    void SetTableURL(std::string url) {fTableURL = url; }
    //void SetRunNumber()

  protected:
      bool LoadRunConditions();
      bool fUseCondb;
      bool fRunConditionsLoaded;
      bool fIsMC;
      uint64_t fCurrentTS;
      uint64_t fRunNumber;
      std::string fCSVFileName;
      std::string fDBTag;
      std::string fTableName;
      std::string fTableURL;

      std::map<int,RunCond_t> fRunCond;

  }; // class RunConditionsProtoDUNE
} //namespace runc
#endif // RUNC_RUNCONDITIONSPROTODUNE_H
