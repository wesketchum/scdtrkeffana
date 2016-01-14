#ifndef TEST_USERANALYSIS_CXX
#define TEST_USERANALYSIS_CXX

#include "TrkEffAlg.hh"
#include <iostream>

trkeff::TrkEffAlg::TrkEffAlg()
  : fAlgName("TrkEffAlg")
{}

void trkeff::TrkEffAlg::SetupOutputTree(TTree* tfs_tree){
  fTree = tfs_tree;

  std::string title = fAlgName + " Tree";
  fTree->SetObject(fTree->GetName(),title.c_str());
}

void trkeff::TrkEffAlg::RunAnalysis(){
  PrintInfo();
}

void trkeff::TrkEffAlg::PrintInfo(){
  std::cout << "\n================================== TrkEffAlg ==========================" << std::endl;
  std::cout << "This is a ub_TrkEffAlg class called " << fAlgName << std::endl;
  std::cout << "\tThere is an output tree called "
	    << fTree->GetName() << " (" << fTree->GetTitle() << ")" << std::endl;
  std::cout << "==========================================================================\n" << std::endl;
}

#endif
