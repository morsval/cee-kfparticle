#include "AnaMaker/AnaMaker.h"
#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TDatime.h"
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char *argv[])
{
  // check parameter
  std::cout << "Number of command line arguments: " << argc << std::endl;
  std::cout << "Command line arguments:" << std::endl;
  for (int i = 0; i < argc; i++) 
    std::cout << "Argument " << i << ": " << argv[i] << std::endl;

  if(argc < 2){
    std::cout<< "The program can be excute by:"<<std::endl;
    std::cout<<"./RunAna.bin /RootFile/ /file number/ /debug_level/"<<std::endl;
		std::cout << "Example: ./RunAna.bin ./input 1 0" << std::endl;
    return 0;
  }

  //prepare the parameter
  TString inFile = argv[1];
	int nFiles = 1;
	if(argc >2 )
		nFiles = std::stoi(argv[2]);
  int debg =0;
	if(argc == 4)
		debg = std::stoi(argv[3]);

  //// set out file name
  TDatime ctime;
 // TString outfile = Form("AnaFullsim_y%dm%02dd%02d.root", ctime.GetYear(), ctime.GetMonth(), ctime.GetDay(), ctime.GetHour(), ctime.GetMinute());
  TString outfile = Form("AnaFullsim_y%dm%02dd%02d_h%02dm%02d.root", ctime.GetYear(), ctime.GetMonth(), ctime.GetDay(), ctime.GetHour(), ctime.GetMinute());
  AnaMaker* dataAnalyzer = new AnaMaker(inFile, nFiles, debg); // the 3rd parm is debug level
	dataAnalyzer->DebugLevel(debg); // default:0, highest: 2 (export per event and track)
  dataAnalyzer->SetOutputFile(outfile);
	//dataAnalyzer->SetEventNumber(10);

  TStopwatch timer;
  timer.Start();

  dataAnalyzer->ProcessData();

  timer.Stop();
  Double_t elapsedTime = timer.RealTime();
  int hours = elapsedTime / 3600;
  int minutes = (elapsedTime / 60) - (hours * 60);
  int seconds = elapsedTime - (hours * 3600) - (minutes * 60);

  std::cout << "Run time: " << std::setfill('0') << std::setw(2) << hours << ":"
        << std::setfill('0') << std::setw(2) << minutes << ":"
        << std::setfill('0') << std::setw(2) << seconds << std::endl;

  return 0;
}
