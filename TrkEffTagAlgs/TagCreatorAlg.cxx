#ifndef TRKEFF_TAGCREATORALG_CXX
#define TRKEFF_TAGCREATORALG_CXX

#include <iostream>

#include "cetlib/exception.h"

#include "RecoBase/Hit.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/GeometryCore.h"
#include "RecoAlg/DBScanAlg.h"
#include "Utilities/DetectorProperties.h"

#include "TagCreatorAlg.hh"


trkeff::TagCreatorAlg::TagCreatorAlg()
{}

void trkeff::TagCreatorAlg::SetupOutputTree(TTree* tfs_tree){
  fTree = tfs_tree;
  fTree->SetObject(fTree->GetName(),"TagCreatorAlg Tree");
}

void trkeff::TagCreatorAlg::FillConfigParameters(fhicl::ParameterSet const& p)
{
  fSearchRegions = p.get< std::vector< std::vector<double> > >("SearchRegions");
  fTagWiresPerPlane = p.get< std::vector<unsigned int> >("TagWiresPerPlane");
  fLineMaxChiSquare = p.get< double >("LineMaxChiSquare");
  fTimeMatch = p.get< double >("TimeMatch");
  fMinHitAmplitudes = p.get< std::vector<double> >("MinHitAmplitudes");
  fMaxHitAmplitudes = p.get< std::vector<double> >("MaxHitAmplitudes");
  fDBScan.reconfigure(p.get< fhicl::ParameterSet >("DBScanAlg"));
  fDebug = p.get<bool>("Debug",false);
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
      << "TagWiresPerPlane must have size as nplanes =" << geo.Nplanes() << ".\n";

  if(fMinHitAmplitudes.size()==1)
    fMinHitAmplitudes = std::vector<double>(geo.Nplanes(),fMinHitAmplitudes[0]);
  if(fMinHitAmplitudes.size()!=geo.Nplanes())
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "fMinHitAmplitudes must have same size as nplanes =" << geo.Nplanes() << ".\n";

  if(fMaxHitAmplitudes.size()==1)
    fMaxHitAmplitudes = std::vector<double>(geo.Nplanes(),fMaxHitAmplitudes[0]);
  if(fMaxHitAmplitudes.size()!=geo.Nplanes())
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "fMinHitAmplitudes must have same size as nplanes =" << geo.Nplanes() << ".\n";
  
  for(size_t i_p=0; i_p<fMaxHitAmplitudes.size(); ++i_p){
    if(fMaxHitAmplitudes[i_p]<fMinHitAmplitudes[i_p])
      throw cet::exception("trkeff::TagCreatorAlg::Configure")
	<< "MaxHitAmplitude must be greater than MinHitAmplitude\n";
    if(fMinHitAmplitudes[i_p]<0)
      throw cet::exception("trkeff::TagCreatorAlg::Configure")
	<< "MinHitAmplitude must be greater than zero\n";
  }

  fSearchRegionsWires.resize(fSearchRegions.size(),WireIDRegionByPlane_t(geo.Nplanes()));
  for (size_t i_s=0; i_s<fSearchRegions.size(); ++i_s)
    TranslateSearchRegion(i_s,geo);
  
  fSortedHitsIndex.resize(fSearchRegionsWires.size(),HitMapByPlane_t(geo.Nplanes()));
    
  
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

  for(size_t i_p=0; i_p<fSearchRegionsWires[i_s].size(); ++i_p){
    geo::WireID start_wire = geo.NearestWireID(start_point,i_p);
    geo::WireID end_wire = geo.NearestWireID(end_point,i_p);
    fSearchRegionsWires[i_s][i_p][0] = std::min(start_wire,end_wire);
    fSearchRegionsWires[i_s][i_p][1] = std::max(start_wire,end_wire);
  }
  
}

void trkeff::TagCreatorAlg::PrintSearchRegionsWires(){

  std::cout << "YZ Search regions: " << fSearchRegions.size()
	    << "\tWire search regions: " << fSearchRegionsWires.size()
	    << std::endl;

  if(fSearchRegions.size()!=fSearchRegionsWires.size())
    throw cet::exception("trkeff::TagCreatorAlg::PrintSearchRegionsWires")
      << "Search regions not of similar size\n";

  
  for(size_t i_s=0; i_s<fSearchRegions.size(); ++i_s){

    std::cout << "Search region " << i_s << ": "
	      << "\ty1,z1,y2,z2, = ";
    for(auto const& var : fSearchRegions[i_s])
      std::cout << var << ",";
    for(size_t i_p=0; i_p<fSearchRegionsWires[i_s].size(); ++i_p)
      std::cout << "\n\tPlane " << i_p << ": ["
		<< fSearchRegionsWires[i_s][i_p][0] << ","
		<< fSearchRegionsWires[i_s][i_p][1];
    std::cout << std::endl;
       
  }

}

void trkeff::TagCreatorAlg::PrintHitsBySearchRegion(){

  std::cout << "Wire search regions: " << fSearchRegionsWires.size()
	    << "Sorted hits collection: " << fSortedHitsIndex.size()
	    << std::endl;

  if(fSearchRegionsWires.size()!=fSortedHitsIndex.size())
    throw cet::exception("trkeff::TagCreatorAlg::PrintHitsBySearchRegion")
      << "Sorted hits not same size as search regions\n";

  
  for(size_t i_s=0; i_s<fSearchRegionsWires.size(); ++i_s){

    std::cout << "Search region " << i_s << ": " << std::endl;
    for(size_t i_p=0; i_p<fSearchRegionsWires[i_s].size(); ++i_p){
      std::cout << "\tPlane " << i_p << " ["
		<< fSearchRegionsWires[i_s][i_p][0] << ","
		<< fSearchRegionsWires[i_s][i_p][1] << "] : " << std::endl;
      for(size_t i_h=0; i_h<fSortedHitsIndex[i_s][i_p].size(); ++i_h)
	std::cout << "\t\tHit index " << fSortedHitsIndex[i_s][i_p][i_h] << std::endl;
    }
       
  }

}


void trkeff::TagCreatorAlg::Configure( fhicl::ParameterSet const& p, geo::GeometryCore const& geo){
  FillConfigParameters(p);
  ProcessConfigParameters(geo);
  if(fDebug) PrintSearchRegionsWires();
}

void trkeff::TagCreatorAlg::Cleanup(){
  for(auto & hits_per_search_reg : fSortedHitsIndex)
    for(auto & hits_per_plane : hits_per_search_reg)
      hits_per_plane.clear();

}

void trkeff::TagCreatorAlg::CreateTags( std::vector<recob::Hit>  const& hit_collection,
					geo::GeometryCore        const& geom,
					util::DetectorProperties & detprop,
					util::LArProperties      const& larprop){
  
  SortHitsBySearchRegion(hit_collection);
  if(fDebug) {
    std::cout << "Hits per search region before pruning by time." << std::endl;
    PrintHitsBySearchRegion();
  }

  for(auto & sr : fSortedHitsIndex){
    RemoveHitsWithoutTimeMatch(hit_collection,sr,detprop);
    if(fDebug) {
      std::cout << "Hits per search region after pruning by time." << std::endl;
      PrintHitsBySearchRegion();
    }
  }

  
  for(auto const& sr : fSortedHitsIndex){
    for(auto const& hm : sr){
      
      //create temp vector ...
      std::vector<size_t> hit_indices; hit_indices.reserve(hm.size());
      for(auto const& ih : hm)
	hit_indices.emplace_back(ih.second);
      ClusterHits(hit_collection,hit_indices,geom,detprop,larprop);
    }
  }

}

void trkeff::TagCreatorAlg::SortHitsBySearchRegion(std::vector<recob::Hit> const& hit_collection){

  for(size_t i_h=0; i_h<hit_collection.size(); i_h++){
    
    size_t i_p = hit_collection[i_h].WireID().planeID().Plane;

    for(size_t i_s=0; i_s<fSearchRegionsWires.size(); ++i_s){
      if(hit_collection[i_h].WireID().Wire>=fSearchRegionsWires[i_s][i_p][0].Wire &&
	 hit_collection[i_h].WireID().Wire<=fSearchRegionsWires[i_s][i_p][1].Wire &&
	 hit_collection[i_h].PeakAmplitude() > fMinHitAmplitudes[i_p] &&
	 hit_collection[i_h].PeakAmplitude() < fMaxHitAmplitudes[i_p])
	fSortedHitsIndex[i_s][i_p][hit_collection[i_h].PeakTime()] = i_h;
    }
  }
  
}

void trkeff::TagCreatorAlg::RemoveHitsWithoutTimeMatch(std::vector<recob::Hit> const& hit_collection,
						       HitMapByPlane_t & hitmaps,
						       util::DetectorProperties & detprop){

  if(hitmaps.size()<2) return;

  HitMapByPlane_t new_hitmaps(hitmaps.size()-1);
  
  for(HitMap_t::iterator iter_hit=hitmaps.back().begin();
      iter_hit!=hitmaps.back().end();
      ++iter_hit){

    for(size_t i_p=0; i_p<new_hitmaps.size(); ++i_p){

      for(HitMap_t::reverse_iterator iter_hit_p =
	    HitMap_t::reverse_iterator(std::next(hitmaps[i_p].lower_bound(iter_hit->first)));
	  iter_hit_p!=hitmaps[i_p].rend();
	  ++iter_hit_p)
	{
	  if(std::abs( (iter_hit->first-detprop.GetXTicksOffset((int)(hitmaps.size()-1),0,0)) - 
		       (iter_hit_p->first-detprop.GetXTicksOffset((int)i_p,0,0))
		       )>fTimeMatch)
	    break;
	  new_hitmaps[i_p].insert(*iter_hit_p);
	}
      for(HitMap_t::iterator iter_hit_p =
	    hitmaps[i_p].upper_bound(iter_hit->first);
	  iter_hit_p!=hitmaps[i_p].end();
	  ++iter_hit_p)
	{
	  if(std::abs( (iter_hit->first-detprop.GetXTicksOffset((int)(hitmaps.size()-1),0,0)) - 
		       (iter_hit_p->first-detprop.GetXTicksOffset((int)i_p,0,0))
		       )>fTimeMatch)
	    break;
	  new_hitmaps[i_p].insert(*iter_hit_p);
	}


    }
    
  }

  for(size_t i_p=0; i_p<new_hitmaps.size(); ++i_p)
    hitmaps[i_p].swap(new_hitmaps[i_p]);
  

}

std::vector<size_t> trkeff::TagCreatorAlg::ClusterHits( std::vector<recob::Hit> const& hit_collection, std::vector<size_t> const& hit_index,
							geo::GeometryCore        const& geom,
							util::DetectorProperties const& detprop,
							util::LArProperties      const& larprop){

  // Select the hits to cluster, the hit indices are marked by the hit_index vector
  // I'm doing this here temporarily, I'm not sure how much of DBScanAlg should be modified
  std::vector<recob::Hit> hit_collection_cluster;
  for(size_t i_h=0; i_h<hit_index.size(); i_h++){
    hit_collection_cluster.push_back(hit_collection[hit_index[i_h]]);
  }

  // Initialize the fDBScan object, this also clears out any data from previous scans
  fDBScan.InitScan(hit_collection, hit_index, std::set<uint32_t>(),
		   geom,larprop,detprop);

  // Run the algorithm
  fDBScan.run_cluster();

  // Now we fill the hit_cluster vector with cluster numbers for each hit
  // The fpointId_to_clusterId function provides the cluster id
  std::vector<size_t> hit_cluster;
  hit_cluster.resize(hit_collection_cluster.size());
  for(size_t i_h = 0; i_h < fDBScan.fpointId_to_clusterId.size(); ++i_h){	  
    hit_cluster[i_h] = fDBScan.fpointId_to_clusterId[i_h];
  }

  return hit_cluster;

}







#endif
