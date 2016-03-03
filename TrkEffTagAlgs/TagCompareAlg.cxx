#ifndef TRKEFF_TAGCOMPARE_CXX
#define TRKEFF_TAGCOMPARE_CXX

#include "TagCompareAlg.hh"

void trkeff::TagCompareAlg::SetupOutputTree(TTree* tree){

      fTree = tree;
      fTree->Branch("tag",&fTag,fTag.Leaflist().c_str());
      fTree->Branch("mctrk",&fMCTrack,"x_start/D:y_start/D:z_start/D:x_end/D:y_end/D:z_end/D");
      fTree->Branch("diff",&fClosestPts,"x_start/D:y_start/D:z_start/D:dist_start/D:x_end/D:y_end/D:z_end/D:dist_end/D");
  
}

void trkeff::TagCompareAlg::Compare(TrkEffTag const& tag, sim::MCTrack const& mctrk){
  fTag = tag.GetRootTreeType();
  FillMCTrackStruct(mctrk);
  FillClosestPtsStruct(tag,mctrk);
  fTree->Fill();
}

void trkeff::TagCompareAlg::FillMCTrackStruct(sim::MCTrack const& mctrk){
  fMCTrack.x_start = mctrk.Start().X();
  fMCTrack.y_start = mctrk.Start().Y();
  fMCTrack.z_start = mctrk.Start().Z();
  fMCTrack.x_end = mctrk.End().X();
  fMCTrack.y_end = mctrk.End().Y();
  fMCTrack.z_end = mctrk.End().Z();
}

void trkeff::TagCompareAlg::FillClosestPtsStruct(TrkEffTag const& tag, sim::MCTrack const& mctrk){

  double mindsq_start=9e9;
  double mindsq_end=9e9;

  double dist;
  
  for(auto const& step : mctrk){

    dist = ((step.X()-tag.StartPoint()[0])*(step.X()-tag.StartPoint()[0]) +
	    (step.Y()-tag.StartPoint()[1])*(step.Y()-tag.StartPoint()[1]) +
	    (step.Z()-tag.StartPoint()[2])*(step.Z()-tag.StartPoint()[2]));
    
    if( dist < mindsq_start){
      mindsq_start = dist;
      fClosestPts.x_start = step.X();
      fClosestPts.y_start = step.Y();
      fClosestPts.z_start = step.Z();
    }
	
    dist = ((step.X()-tag.EndPoint()[0])*(step.X()-tag.EndPoint()[0]) +
	    (step.Y()-tag.EndPoint()[1])*(step.Y()-tag.EndPoint()[1]) +
	    (step.Z()-tag.EndPoint()[2])*(step.Z()-tag.EndPoint()[2]));
    
    if( dist < mindsq_end){
      mindsq_end = dist;
      fClosestPts.x_end = step.X();
      fClosestPts.y_end = step.Y();
      fClosestPts.z_end = step.Z();
    }
    
  }

  fClosestPts.dist_start = std::sqrt(mindsq_start);
  fClosestPts.dist_end = std::sqrt(mindsq_end);

  std::cout << "closest pts distances = " << fClosestPts.dist_start << ", " << fClosestPts.dist_end << std::endl;
  
}
  
#endif

