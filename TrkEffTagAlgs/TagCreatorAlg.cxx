#ifndef TRKEFF_TAGCREATORALG_CXX
#define TRKEFF_TAGCREATORALG_CXX

#include <iostream>

#include "cetlib/exception.h"

#include "RecoBase/Hit.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/GeometryCore.h"

#include "TagCreatorAlg.hh"

trkeff::TagCreatorAlg::TagCreatorAlg()
{}

void trkeff::TagCreatorAlg::SetupOutputTree(TTree* tfs_tree){
  fTree = tfs_tree;
  fTree->SetObject(fTree->GetName(),"TagCreatorAlg Tree");
}

void trkeff::TagCreatorAlg::FillConfigParameters(fhicl::ParameterSet const& p){
  fSearchRegions = p.get< std::vector< std::vector<float> > >("SearchRegions");
  fTagWiresPerPlane = p.get< std::vector<unsigned int> >("TagWiresPerPlane");
  fLineMaxChiSquare = p.get< double >("LineMaxChiSquare");
  fPlaneCombinations = p.get< std::vector<std::string> >("PlaneCombinations");
  fTimeMatch = p.get< double >("TimeMatch");
  fMinHitAmplitudes = p.get< std::vector<double> >("MinHitAmplitudes");
  fMaxHitAmplitudes = p.get< std::vector<double> >("MaxHitAmplitudes");
}

void trkeff::TagCreatorAlg::ProcessConfigParameters(geo::GeometryCore const& geo){

  if(fLineMaxChiSquare<=0)
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "LineMaxChiSquare must be greater than zero\n";
  if(fTimeMatch<=0)
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "TimeMatch must be greater than zero\n";

  if(fTagWiresPerPlane.size()==1)
    fTagWiresPerPlane = std::vector<unsigned int>(geo.Nplanes(),fTagWiresPerPlane[0]);
  if(fTagWiresPerPlane.size()!=geo.Nplanes())
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "TagWiresPerPlane must have size Nplanes=" << geo.Nplanes() << ".\n";

  if(fMinHitAmplitudes.size()==1)
    fMinHitAmplitudes = std::vector<double>(geo.Nplanes(),fMinHitAmplitudes[0]);
  if(fMinHitAmplitudes.size()!=geo.Nplanes())
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "fMinHitAmplitudes must have size Nplanes=" << geo.Nplanes() << ".\n";

  if(fMaxHitAmplitudes.size()==1)
    fMaxHitAmplitudes = std::vector<double>(geo.Nplanes(),fMaxHitAmplitudes[0]);
  if(fMaxHitAmplitudes.size()!=geo.Nplanes())
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "fMinHitAmplitudes must have size Nplanes=" << geo.Nplanes() << ".\n";

  for(size_t i_p=0; i_p<geo.Nplanes(); ++i_p){
    if(fMaxHitAmplitudes[i_p]<fMinHitAmplitudes[i_p])
      throw cet::exception("trkeff::TagCreatorAlg::Configure")
	<< "MaxHitAmplitude must be greater than MinHitAmplitude\n";
    if(fMinHitAmplitudes[i_p]<0)
      throw cet::exception("trkeff::TagCreatorAlg::Configure")
	<< "MinHitAmplitude must be greater than zero\n";
  }
  
  TranslateSearchRegions(geo);
  //ValidatePlaneCombinations();
  
    
  
}

void trkeff::TagCreatorAlg::TranslateSearchRegions(geo::GeometryCore const&){
}

void trkeff::TagCreatorAlg::Configure( fhicl::ParameterSet const& p, geo::GeometryCore const& geo){
  FillConfigParameters(p);
  ProcessConfigParameters(p);
}


#endif
