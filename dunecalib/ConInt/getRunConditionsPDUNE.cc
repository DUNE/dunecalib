#include "dunecalib/ConInt/RunConditionsProtoDUNE.h"
#include "nuevdb/IFDatabase/Table.h"
#include <getopt.h>
#include <iostream>


int gRun = -1;
std::string gTableName = "pdunesp.test";
int gVerbosity = 1;
//------------------------------------------------------------
void PrintUsage() {
std::cout << "Usage: getRunConditionsPDUNE -r|--run [run number] -t|--table [db table name]" << std::endl;
}

//------------------------------------------------------------
bool ParseCLArgs(int argc, char* argv[]) {
  struct option long_options[] = {
    {"help",   0, 0, 'h'},
    {"run",   0, 0, 'r'},
    {"table",   0, 0, 't'},
    {0,0,0,0}
  };

  while (1) {
    int optindx;

    int c = getopt_long(argc,argv,"hr:d:f:t:",long_options,&optindx);
        
    if (c==-1) break;
    
    switch(c) {
    case 'r':
      {
	int run = atoi(optarg);
	if (run < 0) {
	  std::cout << "Invalid run number." << std::endl;
	  exit(0);
	}
	gRun = run;
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
  runCond->SetTableURL("https://dbdata0vm.fnal.gov:9443/dune_runcon_prod/");
  runCond->SetTableName(gTableName);
  runCond->SetVerbosity(gVerbosity);
  runCond->UpdateRN(gRun); //When using run based tables VT when using time based
  runCond->LoadConditionsT();
  runCond->GetRunNumber();

  std::cout << "Run Conditions for channel 0:" << std::endl;
  condb::RunCond_t rc = runCond->GetRunConditions(gRun);
  std::cout << "\tStart time = " << rc.start_time
            << "\n\tdata type = " << rc.data_type
            << "\n\trun Number/sofw = " << rc.run_number
  	    << "\n\tupload time = " << rc.upload_t
            << "\n\tsoftware version = " << rc.software_version
            << "\n\tstop_time = " << rc.stop_time 
            << "\n\tbuffer = " << rc.buffer
            << "\n\tac_couple = " << rc.ac_couple
            << "\n\trun type = " << rc.run_type << std::endl;

  delete runCond;
  
  return 0;
}





