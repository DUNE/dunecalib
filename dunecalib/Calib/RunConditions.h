////////////////////////////////////////////////////////////////////////
// \file RunConditions.h
//
// \brief pure virtual base interface for run conditions
//
// \author avizcaya@colostate.edu
// 
////////////////////////////////////////////////////////////////////////
#ifndef RUNC_RUNCONDITIONS_H
#define RUNC_RUNCONDITIONS_H

namespace runc {
  
  class RunConditions {
    
  public:
    
    RunConditions(const RunConditions &) = delete;
    RunConditions(RunConditions &&) = delete;
    RunConditions& operator = (const RunConditions &) = delete;
    RunConditions& operator = (RunConditions &&) = delete;
    virtual ~RunConditions() = default;

  protected:
    RunConditions() = default;

  }; // class RunConditions
} //namespace runc
#endif // RUNC_RUNCONDITIONS_H
