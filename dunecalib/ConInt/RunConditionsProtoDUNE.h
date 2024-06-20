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
    std::string detector_id;
    std::string software_version;
    int buffer;
    std::string ac_couple;
    int baseline;
    bool enabled;
    float gain;
    bool gain_match;
    float leak;
    bool leak_10x;
    float leak_f;
    float peak_time;
    int pulse_dac;
    int strobe_delay;
    int strobe_length;
    int strobe_skip;
    bool test_cap;
    bool adc_test_pattern;
    bool cold;
    std::string detector_type;
    bool pulser;
  } RunCond_t; //Name your structure

  RunCond_t ResetRunCond_t(RunCond_t rct);

  class RunConditionsProtoDUNE : public Conditions {

  public:

    RunConditionsProtoDUNE() {InitialfVal();};
    RunConditionsProtoDUNE(fhicl::ParameterSet const& pset) {InitialfVal();};
    RunConditionsProtoDUNE(RunConditionsProtoDUNE const&) = delete;
    virtual ~RunConditionsProtoDUNE() = default;

    //I think most of this functions can be in the inheritance class
    bool Configure(fhicl::ParameterSet const& pset);
    bool UpdateTS(uint64_t ts=0);
    bool UpdateRN(float rn=0);

    bool LoadConditionsT();
    RunCond_t GetRunConditions(float chanId);


  protected:
    bool LoadRunConditions();
    
    std::map<int,RunCond_t> fRunCond;    
    std::map<float,RunCond_t> fRunCondR;
  }; //class RunConditionsProtoDUNE
} //namespace condb
#endif // CONDB_RUNCONDITIONSPROTODUNE_H

