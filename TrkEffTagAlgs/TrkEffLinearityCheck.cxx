#ifndef TRKEFF_TRKEFFLINEARITYCHECK_CXX
#define TRKEFF_TRKEFFLINEARITYCHECK_CXX

#include <iostream>

#include "cetlib/exception.h"

#include "RecoBase/Hit.h"

#include "TrkEffLinearityCheck.hh"

// ROOT includes 
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TPrincipal.h>
#include <TGraph.h>
#include <TF1.h>

trkeff::TrkEffLinearityCheck::TrkEffLinearityCheck()
{}


void trkeff::TrkEffLinearityCheck::LineFitCheck( std::vector<recob::Hit> const& hit_collection, 
							std::vector<size_t> const& hit_index){

  std::vector<Float_t> wireVector;
  std::vector<Float_t> tickVector;
  // Look at a fit to a line via a TGraph
  for(size_t i_h = 0; i_h < hit_index.size(); i_h++){
    wireVector.push_back(hit_collection[i_h].WireID().Wire);
    tickVector.push_back(0.5*(hit_collection[i_h].PeakTimeMinusRMS()+hit_collection[i_h].PeakTimePlusRMS()));
  }
  TGraph *gClusterHitsVsWires = new TGraph(wireVector.size(),&wireVector[0],&tickVector[0]);
  // Fit the graph
  TF1 *fLinear = new TF1("fLinear","pol1");
  gClusterHitsVsWires->Fit("fLinear");
  fLinear->GetChisquare();
  fLinear->GetNDF();
  delete gClusterHitsVsWires;

}




#endif
