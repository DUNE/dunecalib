#include "dunecalib/ConInt/RunConditionsProtoDUNE.h"
#include "nuevdb/IFDatabase/Table.h"
#include <getopt.h>
#include <iostream>

// Set the variables here or some in terminal
float gRun = -1; //first run number
float gRun1 = 0; //not necesary, last run number
std::string gDBTag = "v1.4"; // version of the database table
std::string gTableName = "pdunesp.run_conditionstest"; //name of the database table as: schema.tablename pdunesp.test
int gVerbosity = 1; //how much output to get where 0 is none, 1 more...
std::string gTableURL = "https://dbdata0vm.fnal.gov:9443/dune_runcon_prod/"; // table url

//------------------------------------------------------------
void PrintUsage() {
std::cout << "Usage: getRunConditionsPDUNE -r|--run [run number] -f|--runf [run number f] -t|--table [db table name]" << std::endl;
}

//------------------------------------------------------------
bool ParseCLArgs(int argc, char* argv[]) {
  struct option long_options[] = {
    {"help",  0, 0, 'h'},
    {"run",   0, 0, 'r'},
    {"runf",  0, 0, 'f'},
    {"table", 0, 0, 't'},
    {0,0,0,0}
  };

  while (1) {
    int optindx;

    int c = getopt_long(argc,argv,"hr:d:f:t:",long_options,&optindx);
        
    if (c==-1) break;
    
    switch(c) {
    case 'r':
      {
	float run = atoi(optarg);
	if (run < 0) {
	  std::cout << "Invalid run number." << std::endl;
	  exit(0);
	}
	gRun = run;
	break;
      }
    case 'f':
      {
        float runf = atoi(optarg);
        if (runf < 0) {
          std::cout << "Invalid runf number." << std::endl;
          exit(0);
        }
        gRun1 = runf;
        break;
      }
    case 't':
      {
        gTableName = optarg;
        break;
      }
    case 'h':
    default:
      {
	return false;
      }
    }
  }
  if (gRun<0)
    return false;
  return true;
}


//------------------------------------------------------------
int main(int argc, char **argv)
{
  if (!ParseCLArgs(argc,argv)) {
    PrintUsage();
    return 1;
  }
  condb::RunConditionsProtoDUNE* runCond = new condb::RunConditionsProtoDUNE();
  runCond->SetTableURL(gTableURL);
  runCond->SetTableName(gTableName);
  runCond->SetVerbosity(gVerbosity);
  runCond->SetRunNumber1(gRun1);
  runCond->UpdateRN(gRun); //When using tables with run as key (UpdateRN), when using tables with time as key (UpdateVT)

  //runCond->SetTag(gDBTag);
  runCond->LoadConditionsT();
  float run1;
  if (gRun1 > 0) {run1 = gRun1;}
  else {run1 = gRun;}
  for (int run = gRun; run <= run1; run = run + 1) {
    std::cout << "Run Conditions for run number: " <<  std::endl;
    condb::RunCond_t rc = runCond->GetRunConditions(run);
    std::cout << "\tStart time = " << rc.start_time
            << "\n\tdata type = " << rc.data_type
            << "\n\trun Number/sofw = " << rc.run_number
  	    << "\n\tupload time = " << rc.upload_t
            << "\n\tsoftware version = " << rc.software_version
            << "\n\tstop_time = " << rc.stop_time 
            << "\n\tbuffer = " << rc.buffer
            << "\n\tac_couple = " << rc.ac_couple
            << "\n\thv = " << rc.detector_hv
            << "\n\thv set value= " << rc.detector_hvset
            << "\n\trun type = " << rc.run_type << std::endl;
  }
  delete runCond;
  
  return 0;
}





