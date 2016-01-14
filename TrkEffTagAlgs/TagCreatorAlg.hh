#ifndef TRKEFF_TAGCREATORALG_H
#define TRKEFF_TAGCREATORALG_H

#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"

#include "TTree.h"

#include "RecoBase/Hit.h"

namespace trkeff{
  class TagCreatorAlg;
}

class trkeff::TagCreatorAlg{
  
public:
  
  /// Default constructor
  TagCreatorAlg();

  void Configure( fhicl::ParameterSet const& );

  void CreateTags( std::vector<recob::Hit> const&);
  
  /// Default destructor
  virtual ~TagCreatorAlg(){};

  void SetupOutputTree(TTree*);
  
 private:

  std::vector< std::vector<float> > fSearchRegions; //vector of [y1, z1, y2, z2]
  std::vector< int >                fTagWiresPerPlane; //search size, in wires, for each plane
  float                             fLinearity; //min linearity for tag
  std::vector< std::string >        fPlaneCombinations; //allowed plane combinations
  double                            fTimeMatch; //time matching condition for hits across planes

  std::vector< std::vector< std::vector<size_t> > > fSortedHitsIndex;
  void SortHitsBySearchRegion(std::vector<recob::Hit> const&);

  std::vector<size_t> ClusterHits( std::vector<size_t> const& );
  
  TTree*      fTree;
  
};

#endif
