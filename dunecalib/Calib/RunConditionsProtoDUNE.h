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
    float run_number;
    std::string data_type;
    float upload_t;
    //std::string run_type;
    //float chi2;
    //int adc_low;
    //int adc_high;
    //float nl[20];
  } RunCond_t;

  class RunConditionsProtoDUNE : public RunConditions {
    
  public:

    RunConditionsProtoDUNE();
    RunConditionsProtoDUNE(fhicl::ParameterSet const& pset);
    RunConditionsProtoDUNE(RunConditionsProtoDUNE const&) = delete;
    virtual ~RunConditionsProtoDUNE() = default;
      
    bool Configure(fhicl::ParameterSet const& pset);
    bool Update(uint64_t ts=0);
    
    RunCond_t GetRunConditions(int chanId); 
    
    void SetIsMC(bool v) { fIsMC = v; }
    void SetUseCondb(bool v) { fUseCondb = v; }
    void SetCSVFileName(std::string f) { fCSVFileName = f; }
    void SetTag(std::string ta) { fDBTag = ta; }
    void SetTableName(std::string tn) {fTableName = tn; }

  protected:
      bool LoadRunConditions();

      bool fUseCondb;
      bool fRunConditionsLoaded;
      bool fIsMC;
      uint64_t fCurrentTS;
      std::string fCSVFileName;
      std::string fDBTag;
      std::string fTableName;

      std::map<int,RunCond_t> fRunCond;

  }; // class RunConditionsProtoDUNE
} //namespace runc
#endif // RUNC_RUNCONDITIONSPROTODUNE_H
