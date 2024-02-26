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

#include "wda.h"

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

  // crate a class to access the database and get info
  class Table {
  
  public: 

    Table();
    Table(std::string tableName, std::string tableURL);
    virtual ~Table() = default;


    void SetFolderName(std::string folder);
    void SetFolderURL(std::string url) {fConDBURL = url;};
    bool LoadConditionsTable();

  private:

    //bool LoadConditionsTable();
    bool GetDataFromWebService(Dataset&, std::string);

    std::string fFolder;
    std::string fConDBURL;
    short       fVerbosity;
    std::string fTag;

    int     fConnectionTimeout;
    double  fMaxTSVld;
    double  fMinTSVld;

  }; // class Table
} //namespace runc
#endif // RUNC_RUNCONDITIONS_H
