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

  //runCond->SetIsMC(false); //(gDataType == "data"));
  runCond->SetUseCondb(true);
  if (! gCSVFile.empty())
    runCond->SetCSVFileName(gCSVFile);
  runCond->SetTableName(gTableName);
  runCond->Update(gRun);
  runc::RunCond_t rc = runCond->GetRunConditions(0);
  std::cout << "Run Conditions for channel 0:" 
	    << std::endl; 
  std::cout << "\tstart time = " << rc.start_time
            << "\n\t data type = " << rc.data_type
            << "\n\tRun Number/sofw = " << rc.run_number
  	    << "\n\tupload time = " << rc.upload_t << std::endl;

  delete runCond;
  
  return 0;
}

