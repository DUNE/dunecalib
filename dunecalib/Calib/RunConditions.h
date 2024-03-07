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
#include "nuevdb/IFDatabase/Table.h"

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



    nutools::dbi::Row* const GetRow(int i);
    int  AddCol(std::string cname, std::string ctype);
    void AddEmptyRows(unsigned int nrow);
    nutools::dbi::Row* const NewRow() { nutools::dbi::Row* r = new nutools::dbi::Row(fCol); return r;}
    int NRow() {return fRow.size();}

    void SetFolderName(std::string folder) {fFolder = folder;};
    void SetFolderURL(std::string url) {fConDBURL = url;};
    void SetMinTSVld(double t) { fMinTSVld = t;}
    void SetMaxTSVld(double t) { fMaxTSVld = t;}
    void SetVerbosity(int i) { fVerbosity = i;}

    bool LoadConditionsTable();

  private:

    //bool LoadConditionsTable();
    bool GetDataFromWebService(Dataset&, std::string);
    bool GetDataFromWebService_dev(Dataset&, std::string);

    std::string fFolder;
    std::string fConDBURL;
    short       fVerbosity;
    std::string fTag;

    int     fConnectionTimeout;
    double  fMaxTSVld;
    double  fMinTSVld;

    std::vector<nutools::dbi::ColumnDef> fCol;
    std::vector<nutools::dbi::Row>    fRow;

  }; // class Table
} //namespace runc
#endif // RUNC_RUNCONDITIONS_H
