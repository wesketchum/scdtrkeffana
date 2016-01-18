#ifndef TRKEFF_TAGCREATORALG_H
#define TRKEFF_TAGCREATORALG_H

#include <string>
#include <vector>

#include "fhiclcpp/ParameterSet.h"

#include "TTree.h"

//#include "RecoBase/Hit.h"
#include "SimpleTypesAndConstants/geo_types.h"

namespace recob{ class Hit; }
namespace geo{ class GeometryCore; }

namespace trkeff{
  class TagCreatorAlg;
}

class trkeff::TagCreatorAlg{
  
public:
  
  /// Default constructor
  TagCreatorAlg();

  void Configure( fhicl::ParameterSet const&, geo::GeometryCore const& );

  void CreateTags( std::vector<recob::Hit> const&);
  
  /// Default destructor
  virtual ~TagCreatorAlg(){};

  void SetupOutputTree(TTree*);
  
 private:

  //configs from fhicl file
  std::vector< std::vector<double> > fSearchRegions;        // vector of [y1, z1, y2, z2]
  std::vector< unsigned int >        fTagWiresPerPlane;     // search size, in wires, for each plane
  double                             fLineMaxChiSquare;     // max chi squared allowed, per dof, for tag hits being in a line
  std::vector< std::string >         fViewCombinations_str; // allowed plane combinations
  double                             fTimeMatch;            // time matching condition for hits across planes
  std::vector< double >              fMinHitAmplitudes;     // min hit amplitudes per plane
  std::vector< double >              fMaxHitAmplitudes;     // max hit ampltidues per plane
  bool                               fDebug; //run functions for debugging

  
  //internal data
  unsigned int fNplanes;
  std::vector< std::vector<geo::View_t> > fViewCombinations;
  std::vector< std::vector< std::vector<geo::WireID> > > fSearchRegionsWires; //per search region, per plane, [start,end]
  std::vector< std::vector< std::vector<size_t> > > fSortedHitsIndex; //per search region, per plane

  //internal functions
  void FillConfigParameters(fhicl::ParameterSet const&);
  void ProcessConfigParameters(geo::GeometryCore const&);
  void TranslateSearchRegion(size_t, geo::GeometryCore const&);
  void SetViewCombination(std::string const&,geo::GeometryCore const&);
  void SortHitsBySearchRegion(std::vector<recob::Hit> const&);
  std::vector<size_t> ClusterHits( std::vector<recob::Hit> const&, std::vector<size_t> const& );
  void Cleanup();

  //checking functions
  void PrintSearchRegionsWires();
  void PrintHitsBySearchRegion();

  
  TTree*      fTree;
  
};

#endif
