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

#include "TF1.h"
#include "TGraphErrors.h"

trkeff::TagCreatorAlg::TagCreatorAlg()
{
  if(fDebugCanvas){
    fCanvas = new TCanvas("canvas","Tag Creator Debug Canvas",600,600);
  }
}

void trkeff::TagCreatorAlg::SetupOutputTree(TTree* tfs_tree)
{
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
  fMaxHitWidths     = p.get< std::vector<double> >("MaxHitWidths");
  fDBScan.reconfigure(p.get< fhicl::ParameterSet >("DBScanAlg"));
  fDebug = p.get<bool>("Debug",false);
  fDebugCanvas = p.get<bool>("DebugCanvas",false);
  
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
  
  if(fMaxHitWidths.size()==1)
    fMaxHitWidths = std::vector<double>(geo.Nplanes(),fMaxHitWidths[0]);
  if(fMaxHitWidths.size()!=geo.Nplanes())
    throw cet::exception("trkeff::TagCreatorAlg::Configure")
      << "fMaxHitWidths must have same size as nplanes =" << geo.Nplanes() << ".\n";

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

void trkeff::TagCreatorAlg::PrintHitsBySearchRegion(std::vector<recob::Hit> const& hit_collection){

  std::cout << "Wire search regions: " << fSearchRegionsWires.size()
	    << std::endl
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
      for(auto const& i_h : fSortedHitsIndex[i_s][i_p]){
	std::cout << "\t\tHit: (time,index) = (" << i_h.first << "," << i_h.second << ")" << std::endl;
	recob::Hit const& hit(hit_collection[i_h.second]);
	std::cout << "\t\t\t(PeakTime,RMS,Wire) = (" << hit.PeakTime() << "," << hit.RMS() << "," << hit.WireID().Wire << ")" << std::endl;
      }
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

void trkeff::TagCreatorAlg::DebugCanvas(std::vector<recob::Hit>  const& hit_collection,
					std::vector<size_t> const& hit_indices,
					std::string title){

  fCanvas->cd();
  TGraphErrors graph(hit_indices.size());
  for(size_t ih=0; ih<hit_indices.size(); ++ih){
    graph.SetPoint(ih,
		   hit_collection[hit_indices[ih]].WireID().Wire,
		   hit_collection[hit_indices[ih]].PeakTime());
    graph.SetPointError(ih,
			0.5,
			hit_collection[hit_indices[ih]].RMS());
  }
  title += ";Wire number;Time (ticks)";
  graph.SetTitle(title.c_str());
  graph.SetMarkerStyle(7);
  fCanvas->cd();
  graph.Draw("AP");
  fCanvas->Update();

  std::string response;
  std::cout << "Proceed? " << std::endl;
  std::cin >> response;
}


void trkeff::TagCreatorAlg::DebugCanvas(std::vector<recob::Hit>  const& hit_collection,
					std::vector<size_t> const& hit_indices,
					trkeff::LinearLeastSquaresFit::LeastSquaresResult_t const& result,
					trkeff::LinearLeastSquaresFit::LeastSquaresResult_t const& result_invert,
					std::string title){

  fCanvas->cd();
  TGraphErrors graph(hit_indices.size());
  for(size_t ih=0; ih<hit_indices.size(); ++ih){
    graph.SetPoint(ih,
		   hit_collection[hit_indices[ih]].WireID().Wire,
		   hit_collection[hit_indices[ih]].PeakTime());
    graph.SetPointError(ih,
			0.5,
			hit_collection[hit_indices[ih]].RMS());
  }
  title += ";Wire number;Time (ticks)";
  graph.SetTitle(title.c_str());
  fCanvas->cd();
  graph.Draw("AP");
  graph.SetMarkerStyle(7);
  TF1 line("line","pol1",0,4000);
  line.SetParameters(result.intercept,result.slope);
  TF1 line_invert("line_invert","pol1",0,1e4);
  line_invert.SetParameters(-1*result_invert.intercept/result_invert.slope,1./result_invert.slope);
  line.Draw("sames");
  line_invert.Draw("sames");
  fCanvas->Update();

  std::string response;
  std::cout << "Proceed? " << std::endl;
  std::cin >> response;
}

void trkeff::TagCreatorAlg::CreateTags( std::vector<recob::Hit>  const& hit_collection,
					geo::GeometryCore        const& geom,
					util::DetectorProperties & detprop,
					util::LArProperties      const& larprop){

  if(fDebug) std::cout << "In CreateTags..." << std::endl;

  Cleanup();

  if(fDebug) std::cout << "\tCleanup finished..." << std::endl;

  
  SortHitsBySearchRegion(hit_collection,detprop);
  if(fDebug) {
    std::cout << "Hits per search region before pruning by time." << std::endl;
    PrintHitsBySearchRegion(hit_collection);
  }


  for(auto & sr : fSortedHitsIndex){

    if(fDebugCanvas) {
      if(fDebug) std::cout << "\tPlotting hits in search region..." << std::endl;
      for(auto const& hm : sr){
	if(fDebug) std::cout << "\tHits size is " << hm.size() << std::endl;
	
	std::vector<size_t> hit_indices; hit_indices.reserve(hm.size());
	for(auto const& ih : hm)
	  hit_indices.emplace_back(ih.second);	

	if(fDebug) std::cout << "\tCalling plotter..." << std::endl;
	DebugCanvas(hit_collection,hit_indices,"Hits before pruning");
      }
    }
    
    
    if(fDebug) std::cout << "\tPruning hits..." << std::endl;

  RemoveHitsWithoutTimeMatch(hit_collection,sr);
    if(fDebug) {
      std::cout << "Hits per search region after pruning by time." << std::endl;
	PrintHitsBySearchRegion(hit_collection);
    }
    if(fDebugCanvas) {
      for(auto const& hm : sr){
	std::vector<size_t> hit_indices; hit_indices.reserve(hm.size());
	for(auto const& ih : hm)
	  hit_indices.emplace_back(ih.second);	
	DebugCanvas(hit_collection,hit_indices,"Hits after pruning");
      }
    }
  }
  
  for(auto & sr : fSortedHitsIndex){
    for(auto & hm : sr){
      
      if(fDebug)
	std::cout << "NEW PLANE REGION!" << std::endl;

      while(1){

	//create temp vector ...
	std::vector<size_t> hit_indices; hit_indices.reserve(hm.size());
	for(auto const& ih : hm)
	  hit_indices.emplace_back(ih.second);
	
	//Here is the clustering call...
	//ClusterHits(hit_collection,hit_indices,geom,detprop,larprop);
	
	//or could do least squares fits
	auto result = RawLeastSquaresFit(hit_collection,hit_indices);
	auto result_invert = RawLeastSquaresFit(hit_collection,hit_indices,true);
	
	if(fDebugCanvas)
	  DebugCanvas(hit_collection,hit_indices,result,result_invert,"Least squares fit");
	
	if(result.chi2/result.npts < fLineMaxChiSquare)
	  break;
	
	RemoveHitsBadMatch(hit_collection,hm,result,result_invert);
	
      }
    }
  }
}

void trkeff::TagCreatorAlg::SortHitsBySearchRegion(std::vector<recob::Hit> const& hit_collection,
						   util::DetectorProperties & detprop){
  
  for(size_t i_h=0; i_h<hit_collection.size(); i_h++){
    
    size_t i_p = hit_collection[i_h].WireID().planeID().Plane;

    for(size_t i_s=0; i_s<fSearchRegionsWires.size(); ++i_s){
      if(hit_collection[i_h].WireID().Wire>=fSearchRegionsWires[i_s][i_p][0].Wire &&
	 hit_collection[i_h].WireID().Wire<=fSearchRegionsWires[i_s][i_p][1].Wire &&
	 hit_collection[i_h].PeakAmplitude() > fMinHitAmplitudes[i_p] &&
	 hit_collection[i_h].PeakAmplitude() < fMaxHitAmplitudes[i_p])
	fSortedHitsIndex[i_s][i_p][hit_collection[i_h].PeakTime()-detprop.GetXTicksOffset(i_p,0,0)] = i_h;
    }
  }
  
}

void trkeff::TagCreatorAlg::RemoveHitsWithoutTimeMatch(std::vector<recob::Hit> const& hit_collection,
						       HitMapByPlane_t & hitmaps){

  if(hitmaps.size()<2) return;

  HitMapByPlane_t new_hitmaps(hitmaps.size());

  for(size_t i_p1=0; i_p1<new_hitmaps.size()-1; ++i_p1){
    for(size_t i_p2=1; i_p2<new_hitmaps.size(); ++i_p2){
      
      for(auto const& i_h : hitmaps[i_p1]){
	for(auto iter_hit=hitmaps[i_p2].lower_bound(i_h.first-fTimeMatch);
	    iter_hit!=hitmaps[i_p2].end();
	    iter_hit++){

	  if(std::abs(iter_hit->first-i_h.first)>fTimeMatch)
	    break;
	  new_hitmaps[i_p1].insert(i_h);
	  new_hitmaps[i_p2].insert(*iter_hit);
	}
      }
      
      hitmaps[i_p1].swap(new_hitmaps[i_p1]);
      hitmaps[i_p2].swap(new_hitmaps[i_p2]);
      new_hitmaps[i_p1].clear();
      new_hitmaps[i_p2].clear();
      
    }
  }

}

std::vector<unsigned int> trkeff::TagCreatorAlg::ClusterHits( std::vector<recob::Hit> const& hit_collection,
							      std::vector<size_t> const& hit_index,
							      geo::GeometryCore        const& geom,
							      util::DetectorProperties const& detprop,
							      util::LArProperties      const& larprop){

  // Initialize the fDBScan object, this also clears out any data from previous scans
  fDBScan.InitScan(hit_collection, hit_index, std::set<uint32_t>(),
		   geom,larprop,detprop);

  // Run the algorithm
  fDBScan.run_cluster();

  return fDBScan.fpointId_to_clusterId;
  
}

trkeff::LinearLeastSquaresFit::LeastSquaresResult_t
trkeff::TagCreatorAlg::RawLeastSquaresFit(std::vector<recob::Hit> const& hit_collection,
					  std::vector<size_t> const& hit_index,
					  bool invert)
{
  fLSqFit.Clear();
  auto result = fLSqFit.LinearFit(hit_collection,hit_index,invert);
  if(fDebug){
    if(!invert){
      std::cout << "\t\tLeastSquares result: time = " << result.slope
		<< " * wire + " << result.intercept
		<< "\tchi2=" << result.chi2 << " / " << result.npts << std::endl;
    }
    else{
      std::cout << "\t\tLeastSquares result: wire = " << result.slope
		<< " * time + " << result.intercept
		<< "\tchi2=" << result.chi2 << " / " << result.npts << std::endl;
    }
  }
  return result;
}

void trkeff::TagCreatorAlg::RemoveHitsBadMatch(std::vector<recob::Hit> const& hit_collection,
					       HitMap_t & hitmap,
					       trkeff::LinearLeastSquaresFit::LeastSquaresResult_t const& lsq_result,
					       trkeff::LinearLeastSquaresFit::LeastSquaresResult_t const& lsq_result_invert){

  double max_diff=0;
  double max_diff_invert=0;
  double line_value=0;
  double line_value_invert=0;
  auto max_iter = hitmap.begin();
  auto max_iter_invert = hitmap.begin();
  for(auto ih=hitmap.begin(); ih!=hitmap.end(); ++ih){

    line_value = hit_collection[ih->second].WireID().Wire*lsq_result.slope + lsq_result.intercept;
    if(max_diff < std::abs(hit_collection[ih->second].PeakTime()-line_value)){
      max_diff = std::abs(hit_collection[ih->second].PeakTime()-line_value);
      max_iter = ih;
    }

    line_value_invert = hit_collection[ih->second].PeakTime()*lsq_result_invert.slope + lsq_result_invert.intercept;
    if(max_diff_invert < std::abs((double)hit_collection[ih->second].WireID().Wire-line_value_invert)){
      max_diff_invert = std::abs((double)hit_collection[ih->second].WireID().Wire-line_value_invert);
      max_iter_invert = ih;
    }

  }

  if(max_iter!=max_iter_invert){
    hitmap.erase(max_iter);
    hitmap.erase(max_iter_invert);
  }
  else
    hitmap.erase(max_iter);
    
}

#endif
