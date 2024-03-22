////////////////////////////////////////////////////////////////////////
//// \file Condb.h
////
//// \Interface to the conditions database condb2
////
//// \author avizcaya@colostate.edu
//// 
//////////////////////////////////////////////////////////////////////////

#ifndef CONDB_CONDITIONS_H
#define CONDB_CONDITIONS_H

#include "wda.h"
#include "nuevdb/IFDatabase/Table.h"

namespace condb {

  //General conditions class
  class Conditions {
    
  public:
    
    Conditions(const Conditions &) = delete;
    Conditions(Conditions &&) = delete;
    Conditions& operator = (const Conditions &) = delete;
    Conditions& operator = (Conditions &&) = delete;
    virtual ~Conditions() = default;

  protected:
    Conditions() = default;

  }; // class Conditions


  //Table class to access data from specific condb2 table
  class Table {
  
  public: 

    Table();
    Table(std::string tableName, std::string tableURL);
    virtual ~Table() = default;
    
    // using nutools to create/get row with info from all selected columns
    nutools::dbi::Row* const GetRow(int i);
    nutools::dbi::Row* const NewRow() { nutools::dbi::Row* r = new nutools::dbi::Row(fCol); return r;}
    void AddEmptyRows(unsigned int nrow);
    int NRow() {return fRow.size();}
    int AddCol(std::string cname, std::string ctype);
    
    // functions for c++ interface
    void SetFolderName(std::string folder) {fFolder = folder;};
    void SetFolderURL(std::string url) {fConDBURL = url;};
    void SetVerbosity(int i) { fVerbosity = i;}
    void SetMinTSVld(double t) { fMinTSVld = t;}
    void SetMaxTSVld(double t) { fMaxTSVld = t;}
    
    // Load data to the table
    bool LoadConditionsTable();
  
  private:
    // get data from condb2
    bool GetDataFromWebService(Dataset&, std::string);
    bool CheckForNulls();

    std::string fFolder;
    std::string fConDBURL;
    short       fVerbosity;
    int     fConnectionTimeout;
    double  fMaxTSVld;
    double  fMinTSVld;
    std::string fTag;

    std::vector<nutools::dbi::ColumnDef> fCol;
    std::vector<nutools::dbi::Row>    fRow;
    std::vector<std::pair<int,int> > fNullList;

   }; // class Table

} //namespace condb
#endif // CONDB_CONDITIONS_H
