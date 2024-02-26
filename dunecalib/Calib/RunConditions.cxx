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

#include "RunConditions.h"


namespace runc {


  Table::Table() {

    fFolder="";
    fConDBURL="";
    fVerbosity = 1;
    fMaxTSVld = 0;
    fMinTSVld = 0;
    fTag = "";

    // default total time to attempt to connect to the db/web server will be 
    // ~4 minutes (some randomness is introduced by libwda)
    fConnectionTimeout = 4*60; 
  }
 
  //************************************************************
  Table::Table(std::string tableName, std::string tableURL)
  {
    std::cout << "This is the Table Name: " << tableName << std::endl;

    // default total time to attempt to connect to the db/web server will be 
    // ~4 minutes (some randomness is introduced by libwda)
    fConnectionTimeout = 4*60; 
  }


  //************************************************************
  bool Table::GetDataFromWebService(Dataset& ds, std::string myss) 
  {
    Tuple tu;
    //char ss[1024]; 
    //char ss2[1024];
    int wda_err; //, err;
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

	return true;
    }
    if(fVerbosity > 0)
	std::cout << "Got " << ntup-1 << " rows from database" << std::endl;

    tu = getFirstTuple(ds);
    if (tu == NULL) {
	std::cerr << "Table::Load(" << fFolder << ") has NULL first tuple!"
		  << std::endl;
	return false;
    }

    //int ncol2 = getNfields(tu);
    releaseTuple(tu);
    releaseDataset(ds);

    return true;
  }


  //************************************************************
  bool Table::LoadConditionsTable()
  {
    std::stringstream myss;
    //This is just to try, take it out afterwards
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
     
    return GetDataFromWebService(ds,myss.str());
  

  }


}
