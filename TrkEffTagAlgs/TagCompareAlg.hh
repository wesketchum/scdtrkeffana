#ifndef TRKEFF_TAGCOMPARE_H
#define TRKEFF_TAGCOMPARE_H

#include <string>
#include <vector>

#include "TTree.h"
#include "TrkEffTagObjects/TrkEffTag.h"
#include "MCBase/MCTrack.h"

#include "uboone/MuCS/MuCSRecoData.h"
//#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"

namespace trkeff{
  class TagCompareAlg;
}

class trkeff::TagCompareAlg{
  
public:
  
  /// Default constructor
  TagCompareAlg(){};
  void SetupOutputTree(TTree* tree);
  
  void Compare(TrkEffTag const&,sim::MCTrack const&);
  
  /// Default destructor
  virtual ~TagCompareAlg(){};

 private:



  typedef TrkEffTag::TrkEffTag_Tree_t TagStruct_t;

  typedef struct MCTrackStruct{
    double x_start;
    double y_start;
    double z_start;
    double x_end;
    double y_end;
    double z_end;
  } MCTrackStruct_t;

  void FillMCTrackStruct(sim::MCTrack const&);

  typedef struct ClosestPtStruct{
    double x_start;
    double y_start;
    double z_start;
    double dist_start;
    double x_end;
    double y_end;
    double z_end;
    double dist_end;
  } ClosestPtStruct_t;
  
  void FillClosestPtsStruct(TrkEffTag const&,sim::MCTrack const&);

  TTree*            fTree;
  TagStruct_t       fTag;
  MCTrackStruct_t   fMCTrack;
  ClosestPtStruct_t fClosestPts;
  
};

#endif
