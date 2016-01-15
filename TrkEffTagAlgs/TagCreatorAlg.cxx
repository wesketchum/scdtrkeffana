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
  fSearchRegions = p.get< std::vector< std::vector<double> > >("SearchRegions");
  fTagWiresPerPlane = p.get< std::vector<unsigned int> >("TagWiresPerPlane");
  fLineMaxChiSquare = p.get< double >("LineMaxChiSquare");
  fViewCombinations_str = p.get< std::vector<std::string> >("ViewCombinations");
  fTimeMatch = p.get< double >("TimeMatch");
  fMinHitAmplitudes = p.get< std::vector<double> >("MinHitAmplitudes");
  fMaxHitAmplitudes = p.get< std::vector<double> >("MaxHitAmplitudes");
}

void trkeff::TagCreatorAlg::ProcessConfigParameters(geo::GeometryCore const& geo){
  fNplanes = geo.Nplanes();
  
  if(fLineMaxChiSquare<=0)
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "LineMaxChiSquare must be greater than zero\n";
  if(fTimeMatch<=0)
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "TimeMatch must be greater than zero\n";

  if(fTagWiresPerPlane.size()==1)
    fTagWiresPerPlane = std::vector<unsigned int>(fNplanes,fTagWiresPerPlane[0]);
  if(fTagWiresPerPlane.size()!=fNplanes)
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "TagWiresPerPlane must have size Nplanes=" << fNplanes << ".\n";

  if(fMinHitAmplitudes.size()==1)
    fMinHitAmplitudes = std::vector<double>(fNplanes,fMinHitAmplitudes[0]);
  if(fMinHitAmplitudes.size()!=fNplanes)
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "fMinHitAmplitudes must have size Nplanes=" << fNplanes << ".\n";

  if(fMaxHitAmplitudes.size()==1)
    fMaxHitAmplitudes = std::vector<double>(fNplanes,fMaxHitAmplitudes[0]);
  if(fMaxHitAmplitudes.size()!=fNplanes)
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "fMinHitAmplitudes must have size Nplanes=" << fNplanes << ".\n";

  for(size_t i_p=0; i_p<geo.Nviews(); ++i_p){
    if(fMaxHitAmplitudes[i_p]<fMinHitAmplitudes[i_p])
      throw cet::exception("trkeff::TagCreatorAlg::Configure")
	<< "MaxHitAmplitude must be greater than MinHitAmplitude\n";
    if(fMinHitAmplitudes[i_p]<0)
      throw cet::exception("trkeff::TagCreatorAlg::Configure")
	<< "MinHitAmplitude must be greater than zero\n";
  }

  for(auto const& views : fViewCombinations_str)
    SetViewCombination(views,geo);
  
  fSearchRegionsWires.resize(fSearchRegions.size());
  for (size_t i_s=0; i_s<fSearchRegions.size(); ++i_s)
    TranslateSearchRegion(i_s,geo);
  
  fSortedHitsIndex.resize(fSearchRegionsWires.size(),std::vector< std::vector<size_t> >(geo.Nplanes()));
    
  
}

void trkeff::TagCreatorAlg::SetViewCombination(std::string const& str, geo::GeometryCore const& geo){

  fViewCombinations.emplace_back(geo.Nviews());
  
  for(auto const& i_char : str){
    geo::View_t my_view;
    if(i_char=='U' || i_char=='u')
      my_view = geo::View_t::kU;
    else if(i_char=='V' || i_char=='v')
      my_view = geo::View_t::kV;
    else if(i_char=='Y' || i_char=='y'
	    || i_char=='W' || i_char=='w'
	    || i_char=='Z' || i_char=='z')
      my_view = geo::View_t::kZ;
    else
      throw cet::exception("trkeff::TagCreatorAlg::SetViewCombinations")
	<< "Don't understand view type " << i_char << ".\n";

    if(geo.Views().count(my_view)==0)
      throw cet::exception("trkeff::TagCreatorAlg::SetViewCombinations")
	<< "Geometry doesn't have view " << my_view << ".\n";
      
    fViewCombinations.back().push_back(my_view);
  }

  if(fViewCombinations.size()<2)
      throw cet::exception("trkeff::TagCreatorAlg::SetViewCombinations")
	<< "Less than 2 views currently not supported for tag creation.\n";
    
}

void trkeff::TagCreatorAlg::TranslateSearchRegion(size_t i_s, geo::GeometryCore const& geo){

  if(fSearchRegions[i_s].size()!=4)
    throw cet::exception("trkeff::TagCreatorAlg::TranslateSearchRegion")
      << "Search region must be of format [y1, z1, y2, z2].\n";

  if(fSearchRegions[i_s][0]>geo.DetHalfHeight())
    fSearchRegions[i_s][0] = geo.DetHalfHeight();
  else if(fSearchRegions[i_s][0]<-1*geo.DetHalfHeight())
    fSearchRegions[i_s][0] = -1*geo.DetHalfHeight();
  if(fSearchRegions[i_s][2]>geo.DetHalfHeight())
    fSearchRegions[i_s][2] = geo.DetHalfHeight();
  else if(fSearchRegions[i_s][2]<-1*geo.DetHalfHeight())
    fSearchRegions[i_s][2] = -1*geo.DetHalfHeight();
  if(fSearchRegions[i_s][1]>geo.DetLength())
    fSearchRegions[i_s][1] = geo.DetLength();
  else if(fSearchRegions[i_s][1]<0)
    fSearchRegions[i_s][1] = 0;
  if(fSearchRegions[i_s][3]>geo.DetLength())
    fSearchRegions[i_s][3] = geo.DetLength();
  else if(fSearchRegions[i_s][3]<0)
    fSearchRegions[i_s][3] = 0;

  if(fSearchRegions[i_s][0]>=fSearchRegions[i_s][2])
    throw cet::exception("trkeff::TagCreatorAlg::TranslateSearchRegion")
      << "Search region must be of format [y1, z1, y2, z2]. y2<y1 not valid.\n";
  if(fSearchRegions[i_s][1]>=fSearchRegions[i_s][3])
    throw cet::exception("trkeff::TagCreatorAlg::TranslateSearchRegion")
      << "Search region must be of format [y1, z1, y2, z2]. z2<z1 not valid.\n";
    
  std::vector<double> start_point{geo.DetHalfWidth(),fSearchRegions[i_s][0],fSearchRegions[i_s][1]};
  std::vector<double> end_point{geo.DetHalfWidth(),fSearchRegions[i_s][2],fSearchRegions[i_s][3]};

  fSearchRegionsWires[i_s].resize(geo.Nplanes());
  for(size_t i_p=0; i_p<geo.Nplanes(); ++i_p){
    geo::WireID start_wire = geo.NearestWireID(start_point,i_p);
    geo::WireID end_wire = geo.NearestWireID(end_point,i_p);
    fSearchRegionsWires[i_s][i_p].push_back(std::min(start_wire,end_wire));
    fSearchRegionsWires[i_s][i_p].push_back(std::max(start_wire,end_wire));
  }
  
}

void trkeff::TagCreatorAlg::Configure( fhicl::ParameterSet const& p, geo::GeometryCore const& geo){
  FillConfigParameters(p);
  ProcessConfigParameters(geo);
}

void trkeff::TagCreatorAlg::Cleanup(){
  for(auto & hits_per_search_reg : fSortedHitsIndex)
    for(auto & hits_per_plane : hits_per_search_reg)
      hits_per_plane.clear();

}

void trkeff::TagCreatorAlg::CreateTags( std::vector<recob::Hit> const& hit_collection){
  SortHitsBySearchRegion(hit_collection);
  
}

void trkeff::TagCreatorAlg::SortHitsBySearchRegion(std::vector<recob::Hit> const& hit_collection){

  for(size_t i_h=0; i_h<hit_collection.size(); i_h++){
      unsigned int planeID = hit_collection[i_h].WireID().planeID().Plane;
      for(size_t i_s=0; i_s<fSearchRegionsWires.size(); ++i_s){
	if(hit_collection[i_h].WireID().Wire>=fSearchRegionsWires[i_s][planeID][0].Wire &&
	   hit_collection[i_h].WireID().Wire<=fSearchRegionsWires[i_s][planeID][1].Wire &&
	   hit_collection[i_h].PeakAmplitude() > fMinHitAmplitudes[planeID] &&
	   hit_collection[i_h].PeakAmplitude() < fMaxHitAmplitudes[planeID])
	  fSortedHitsIndex[i_s][planeID].push_back(i_h);
      }
  }

}



#endif
