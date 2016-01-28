#ifndef TRKEFF_TAGCREATORALG_H
#define TRKEFF_TAGCREATORALG_H

#include <string>
#include <vector>
#include <array>
#include <set>

#include "fhiclcpp/ParameterSet.h"

#include "TTree.h"

//#include "RecoBase/Hit.h"
#include "RecoAlg/DBScanAlg.h"
#include "SimpleTypesAndConstants/geo_types.h"

namespace recob{ class Hit; }
namespace geo{ class GeometryCore; }
//namespace cluster{ class DBScanAlg; }
namespace util{ class DetectorProperties; }

namespace trkeff{
  class TagCreatorAlg;
}

class trkeff::TagCreatorAlg{
  
public:
  
  /// Default constructor
  TagCreatorAlg();

  void Configure( fhicl::ParameterSet const&, geo::GeometryCore const& );

  void CreateTags( std::vector<recob::Hit>  const&,
		   geo::GeometryCore        const&,
		   util::DetectorProperties &,
		   util::LArProperties      const& );
  
  /// Default destructor
  virtual ~TagCreatorAlg(){};

  void SetupOutputTree(TTree*);
  
 private:

  //configs from fhicl file
  std::vector< std::vector<double> > fSearchRegions;         // vector of [y1, z1, y2, z2]
  std::vector< unsigned int >        fTagWiresPerPlane;      // search size, in wires, for each plane in combination
  double                             fLineMaxChiSquare;      // max chi squared allowed, per dof, for tag hits being in a line
  double                             fTimeMatch;             // time matching condition for hits across planes
  std::vector< double >              fMinHitAmplitudes;      // min hit amplitudes per plane in combination
  std::vector< double >              fMaxHitAmplitudes;      // max hit ampltidues per plane in combination
  bool                               fDebug; //run functions for debugging


  typedef std::array<geo::WireID,2> WireIDRegion_t;
  typedef std::vector<WireIDRegion_t> WireIDRegionByPlane_t;

  typedef std::map<float,size_t> HitMap_t;
  typedef std::vector<HitMap_t> HitMapByPlane_t;

  
  //internal data
  std::vector<WireIDRegionByPlane_t> fSearchRegionsWires; //per search region, per plane, [start,end]
  std::vector<HitMapByPlane_t>   fSortedHitsIndex; //per search region, per plane
  std::vector<double> fTimeCorrections; //time corrections per plane

  // dbscan object
  cluster::DBScanAlg fDBScan;

  
  //internal functions
  void FillConfigParameters(fhicl::ParameterSet const&);
  void ProcessConfigParameters(geo::GeometryCore const&);
  void TranslateSearchRegion(size_t, geo::GeometryCore const&);
  void SetPlaneCombination(geo::GeometryCore const&);
  void SortHitsBySearchRegion(std::vector<recob::Hit> const&, util::DetectorProperties &);
  void RemoveHitsWithoutTimeMatch(std::vector<recob::Hit> const&,
				  HitMapByPlane_t & hitmaps);
  std::vector<size_t> ClusterHits( std::vector<recob::Hit> const&, std::vector<size_t> const&,
				   geo::GeometryCore const&,
				   util::DetectorProperties const&,
				   util::LArProperties const&);
  void Cleanup();

  //checking functions
  void PrintSearchRegionsWires();
  void PrintHitsBySearchRegion();

  
  TTree*      fTree;
  
};

#endif
