#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>


void FileListProcess(){
  string line;
  ifstream fstream ("FileList_q2_10to100.txt");
  int FileCount = 0;
  //TChain *AnalysisChain = new TChain("events"); // We could define a chain to process our files too and add them as we scan over our list
  while(getline(fstream, line)){
    if (FileCount > 5) continue; // Stop loop after 5 files, comment out to read full file
    // Check file exists
    TString tmpFile{line};
    auto RootFile = TFile::Open(tmpFile);
    if(!RootFile){ // Check file exists
      cout << "File not found:"<<tmpFile << endl;
      continue;
    }
    cout << "Found file - " << line << endl;
    //AnalysisChain->Add(tmpFile) // Add to our chain if we want
    FileCount++;
  }
}

void read_in_files(){
  FileListProcess();
}
