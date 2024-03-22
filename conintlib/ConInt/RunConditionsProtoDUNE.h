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
#include "TH1F.h"
#include "TH2F.h"

// C/C++ standard libraries
#include <string>
#include <vector>
#include <map>

// dunecalib includes
#include "conintlib/ConInt/Condb.h"

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


  } //class RunConditionsProtoDUNE
} //namespace condb
#endif // CONDB_RUNCONDITIONSPROTODUNE_H

