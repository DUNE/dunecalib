////////////////////////////////////////////////////////////////////////
// \file RunConditionsServiceProtoDUNE.h
//
// \brief header of service for storing/accessing run conditions parameters for ProtoDUNE
//
// \author avizcaya@colostate.edu
// \date Jan 19, 2024
// 
////////////////////////////////////////////////////////////////////////
#ifndef RUNCONDITIONSSERVICEPROTODUNE_H
#define RUNCONDITIONSSERVICEPROTODUNE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "dunecalib/Calib/RunConditionsProtoDUNE.h"
#include "dunecalib/CalibServices/RunConditionsService.h"

namespace runc{
  class RunConditionsServiceProtoDUNE : public RunConditionsService {
    public:
      
      RunConditionsServiceProtoDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);

      virtual void   reconfigure(fhicl::ParameterSet const& pset)  override;
      void   preBeginRun(const art::Run& run); 

      virtual provider_type* provider() const override { return fProp.get();}

    private:

      std::unique_ptr<runc::RunConditionsProtoDUNE> fProp;

    }; // class RunConditionsServiceProtoDUNE
} //namespace runc
DECLARE_ART_SERVICE_INTERFACE_IMPL(runc::RunConditionsServiceProtoDUNE, runc::RunConditionsService, LEGACY)
#endif // RUNCONDITIONSSERVICEPROTODUNE_H
