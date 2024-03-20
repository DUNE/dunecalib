#include "dunecalib/Calib/RunConditionsProtoDUNE.h"

#include "nuevdb/IFDatabase/Table.h"
#include <getopt.h>
#include <iostream>

int gRun = -1;
std::string gDataType = "data";
std::string gCSVFile = "";
std::string gTableName = "";

//------------------------------------------------------------

void PrintUsage()
{
  std::cout << "Usage: getRunConditionsProtoDUNE -r|--run [run number] -d|--datatype [data|mc] -f|--file [csv file] -t|--table [db table name]" << std::endl;
}

//------------------------------------------------------------

bool ParseCLArgs(int argc, char* argv[])
{

  struct option long_options[] = {
    {"help",   0, 0, 'h'},
    {"run",   0, 0, 'r'},
    {"datatype",   0, 0, 'd'},
    {"file",   0, 0, 'f'},
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
    case 'd':
      {
	gDataType = optarg;
	break;
      }
    case 'f':
      {
	gCSVFile = optarg;
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

  if ( gDataType != "mc" && gDataType != "data")
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

  runc::RunConditionsProtoDUNE* runCond = new runc::RunConditionsProtoDUNE();
  runCond->SetTableURL("https://dbdata0vm.fnal.gov:9443/dune_runcon_prod/");
  runCond->SetTableName("pdunesp.test");
  runCond->Update(gRun);
  runCond->LoadConditionsT();


  //Old code
  //runCond->SetUseCondb(true);
  //if (! gCSVFile.empty())
  //  runCond->SetCSVFileName(gCSVFile);
  //runCond->SetTableName(gTableName);
  //runCond->Update(gRun);
  std::cout << "Run Conditions for channel 0:"
            << std::endl;
  runc::RunCond_t rc = runCond->GetRunConditions(0);
  std::cout << "\tstart time = " << rc.start_time
            << "\n\tdata type = " << rc.data_type
            << "\n\tRun Number/sofw = " << rc.run_number
  	    << "\n\tupload time = " << rc.upload_t
            << "\n\tsoftware version = " << rc.software_version
            << "\n\tstop_time = " << rc.stop_time 
            << "\n\tbuffer = " << rc.buffer
            << "\n\tac_couple = " << rc.ac_couple
            << "\n\trun type = " << rc.run_type << std::endl;

  delete runCond;
  
  return 0;
}

