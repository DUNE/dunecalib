////////////////////////////////////////////////////////////////////////
//// ConditionsService.h
////
//// Pure virtual service interface for run conditions functions
////
////  avizcaya@colostate.edu
////  date: Jan 19, 2024
////
//////////////////////////////////////////////////////////////////////////
#ifndef CONDITIONSSERVICE_H
#define CONDITIONSSERVICE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "dunecalib/ConInt/Condb.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace condb{
  class ConditionsService {
  public:
    virtual ~ConditionsService() = default;
    typedef condb::Conditions provider_type;
    virtual condb::Conditions* provider() const = 0; 
   
    virtual void   reconfigure(fhicl::ParameterSet const& pset) = 0;
 
  };
}

DECLARE_ART_SERVICE_INTERFACE(condb::ConditionsService, LEGACY)
#endif // CONDITIONSSERVICE_H
