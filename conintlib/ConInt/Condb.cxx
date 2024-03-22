#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <cctype>
#include <algorithm>
#include <ctime>
#include <cinttypes>
#include <cstdio>

#include "wda.h"
#include "Condb.h"

namespace condb {

  Table::Table() {

    fFolder="";      //table name schema.table
    fConDBURL="";    //url to the db 
    fVerbosity = 1;  //how much output 
    fMaxTSVld = 0;   //Max range to get the data
    fMinTSVld = 0;   //Min range to get the data
    fTag = "";       //Tag of the table, to get older or new version

    // default total time to attempt to connect to the db/web server will be 
    // ~4 minutes (some randomness is introduced by libwda)
    fConnectionTimeout = 4*60; 
    }

  //************************************************************
  //Add just the columns that you want to retrieve data from
  int Table::AddCol(std::string cname, std::string ctype) {
    if (fVerbosity >= 1) {
      std::cout << "Adding column: " << cname << std::endl;
    }
    for (size_t i=0; i<fCol.size(); ++i) {
      if (fCol[i].Name() == cname) {
	std::cerr << "Table::AddCol: column \'" << cname << "\' already exists!  Fatal, aborting..." << std::endl;
	abort();
	}
      }
    
    nutools::dbi::ColumnDef cdef(cname,ctype);
    fCol.push_back(cdef);

    return fCol.size()-1;

  }  

  //************************************************************
  //Adds one row per run number or time
  void Table::AddEmptyRows(unsigned int nrow) {
    nutools::dbi::Row* row = this->NewRow();
    fRow.resize(fRow.size()+nrow,*row);
  }

  //************************************************************
  //Get data from the specified row
  nutools::dbi::Row* const Table::GetRow(int i) { 
    if (i >= 0 && i < (int)fRow.size())
      return &fRow[i];
    else
      return 0;
  }

  //************************************************************
  bool Table::CheckForNulls() {
    bool isOk = fNullList.empty();
    if (!isOk) // print out list of null columns
      for (unsigned int i=0; i<fNullList.size(); ++i)
        if (fVerbosity>0)
          std::cerr << fCol[fNullList[i].second].Name() << " is NULL in row "
                    << fNullList[i].first << std::endl;
    return isOk;
  }

  //************************************************************
  // Gets data from the condb database and returns tuples
  bool Table::GetDataFromWebService(Dataset& ds, std::string myss) {
    Tuple tu;
    char ss[1024]; 
    char ss2[1024];
    int wda_err, err;
    const char* uagent = NULL;

    if(fVerbosity > 0)
      std::cout << "DBWeb query: " << myss << std::endl;
    
    ds = getDataWithTimeout(myss.c_str(), uagent, 
			      fConnectionTimeout, &wda_err);
    int httpStatus = getHTTPstatus(ds);
    if (httpStatus != 200) {
        std::cerr << "Table::Load: Web Service returned HTTP status " 
		  << httpStatus << ": " << getHTTPmessage(ds) << std::endl;
	return false;
    }
    int ntup = getNtuples(ds);

    // Getting no rows back can be legitimate
    if(ntup == 0){
	if(fVerbosity > 0)
	  std::cout << "Got zero rows from database. Is that expected?" << std::endl;
        fRow.clear();
	return true;
    }
    if(fVerbosity > 0)
	std::cout << "Got " << ntup-1 << " rows from database" << std::endl;
    int ioff=fRow.size();
    AddEmptyRows(ntup);

    tu = getFirstTuple(ds);
    if (tu == NULL) {
	std::cerr << "Table::Load(" << fFolder << ") has NULL first tuple!"
		  << std::endl;
	return false;
    }
    
    int ncol2 = getNfields(tu);
    std::vector<int> colMap(ncol2);
    std::vector<bool> isString(ncol2);
    std::vector<bool> isKnownField(ncol2);
    std::string chanStr = "channel";
    std::string tvStr = "tv";
    std::string tvEndStr = "tvend";
    int chanIdx=-1;
    int tvIdx=-1;
    int tvEndIdx=-1;

    for (int i=0; i<ncol2; ++i) {
      getStringValue(tu,i,ss,sizeof(ss),&err);
      if (chanStr == ss)  { chanIdx=i;  continue;}
      if (tvStr == ss)    { tvIdx=i;    continue;}
      if (tvEndStr == ss) { tvEndIdx=i; continue;}
      bool foundMatch=false; 
      for (unsigned int icol=0; icol<fCol.size(); ++icol) {
        if (fCol[icol].Name() == ss) {
           colMap[i] = icol;
           isString[i] = false;
           if (fCol[icol].Type() == "string" || fCol[icol].Type() == "text")
	      isString[i] = true;
           foundMatch=true;
	   break;
        }
      } 
      if (!foundMatch) // this means this field was unexpected, so ignore it downstream
        isKnownField[i] = false;
      else
	isKnownField[i] = true; 	  
    }

    releaseTuple(tu);
    tu = getNextTuple(ds);
    int irow=0;
    std::string nn = "None";
    std::string nu = "NULL";

    while (tu != NULL) {
      for (int i=0; i<ncol2; ++i) {	
        getStringValue(tu,i,ss,sizeof(ss),&err);
        
        if (i == chanIdx) {
	    uint64_t chan = strtoull(ss,NULL,10);
	    fRow[ioff+irow].SetChannel(chan);
	    continue;
        }
        else if (i == tvIdx) {
	    double t1 = strtod(ss,NULL);
	    fRow[ioff+irow].SetVldTime(t1);
	}
	else if (i == tvEndIdx) {
	    double t1 = strtod(ss,NULL);
	    fRow[ioff+irow].SetVldTimeEnd(t1);	    
	}
        else {
          if (isKnownField[i]) {
            if (isString[i] && (ss[0]=='\'' || ss[0]=='\"')) { // remove quotes
              int k = strlen(ss);
	      strncpy(ss2,&ss[1],k-2);
	      ss2[k-2] = '\0';
	      fRow[ioff+irow].Col(colMap[i]).FastSet(ss2);
	    }
            else if (ss==nn) {
               fRow[ioff+irow].Col(colMap[i]).FastSet(nn);  
            }
            else {
	      fRow[ioff+irow].Col(colMap[i]).FastSet(ss); 
            }
          }
        }
      }
      releaseTuple(tu);
      tu = getNextTuple(ds);
      ++irow;
    }
    releaseDataset(ds);
    return true;
  } //Table::GetDataFromWebService

  //************************************************************
  //Get the complete url used to get the data
  bool Table::LoadConditionsTable() {
    if (fMinTSVld == 0 || fMaxTSVld == 0) {
        std::cerr << "Table::LoadConditionsTable-run: Validity time or Run number is not set!" << std::endl;
        return false;
    }

    std::stringstream myss;
    myss << fConDBURL << "get?folder=" << fFolder << "&";
    if (fTag != "") myss << "tag=" << fTag << "&";
    if (fMinTSVld == fMaxTSVld) 
    {
      myss << "t=" << std::setprecision(12) << fMinTSVld;
    }
    else 
    {
      myss << "t0=" << std::setprecision(12) << fMinTSVld << "&t1=" << std::setprecision(12) << fMaxTSVld;
    }
      
    Dataset ds; 
 
    GetDataFromWebService(ds,myss.str());
    return true;    
  } //Table::LoadConditionsTable

} //namespace condb
