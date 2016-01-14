#ifndef TRKEFF_TAGCREATORALG_CXX
#define TRKEFF_TAGCREATORALG_CXX

#include "TagCreatorAlg.hh"
#include <iostream>

trkeff::TagCreatorAlg::TagCreatorAlg()
{}

void trkeff::TagCreatorAlg::SetupOutputTree(TTree* tfs_tree){
  fTree = tfs_tree;
  fTree->SetObject(fTree->GetName(),"TagCreatorAlg Tree");
}

#endif
