////////////////////////////////////////////////////////////////////////
//// \file RunConditionsServiceProtoDUNE.h
////
//// \brief header of service for storing/accessing run conditions parameters for ProtoDUNE
////
//// \author avizcaya@colostate.edu
//// \date Jan 19, 2024
//// 
//////////////////////////////////////////////////////////////////////////
#ifndef RUNCONDITIONSSERVICEPDUNE_H
#define RUNCONDITIONSSERVICEPDUNE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "dunecalib/ConInt/RunConditionsProtoDUNE.h"
#include "dunecalib/ConIntServices/ConditionsService.h"

namespace condb{

  class RunConditionsServicePDUNE : public ConditionsService {
  public:
    RunConditionsServicePDUNE(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    virtual void reconfigure(fhicl::ParameterSet const& pset)  override;
    void   preBeginRun(const art::Run& run); 
    virtual provider_type* provider() const override { return fProp.get();}

  private:
    std::unique_ptr<condb::RunConditionsProtoDUNE> fProp;
  };

} // namespace condb
DECLARE_ART_SERVICE_INTERFACE_IMPL(condb::RunConditionsServicePDUNE, condb::ConditionsService, LEGACY)
#endif // RUNCONDITIONSSERVICEPDUNE_H
