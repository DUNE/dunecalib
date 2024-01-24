////////////////////////////////////////////////////////////////////////
// RunConditionsService.h
//
// Pure virtual service interface for run conditions functions
//
//  avizcaya@colostate.edu
//  date: Jan 19, 2024
//
////////////////////////////////////////////////////////////////////////
#ifndef RUNCONDITIONSSERVICE_H
#define RUNCONDITIONSSERVICE_H

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "dunecalib/Calib/LifetimeCalib.h"
#include "larcore/CoreUtils/ServiceUtil.h"

namespace runc{
    class RunConditionsService {
      public:
      //typedef runc::RunConditions provider_type;

      public:
      virtual ~RunConditionsService() = default;

      virtual void   reconfigure(fhicl::ParameterSet const& pset) = 0;
      //virtual runc::RunConditions* provider() const = 0;

      }; // class RunConditionsService
    } //namespace detinfo
DECLARE_ART_SERVICE_INTERFACE(runc::RunConditionsService, LEGACY)
#endif // RUNCONDITIONSSERVICE_H
