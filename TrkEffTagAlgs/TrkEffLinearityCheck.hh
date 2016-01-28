#ifndef TRKEFF_TRKEFFLINEARITYCHECK_HH
#define TRKEFF_TRKEFFLINEARITYCHECK_HH


#include <string>
#include <vector>
#include <array>
#include <set>

#include "fhiclcpp/ParameterSet.h"

#include "TTree.h"

//#include "RecoBase/Hit.h"
#include "SimpleTypesAndConstants/geo_types.h"

namespace recob{ class Hit; }

namespace trkeff{
  class TrkEffLinearityCheck;
}

class trkeff::TrkEffLinearityCheck{
  
public:
  
  /// Default constructor
  TrkEffLinearityCheck();


  
  /// Default destructor
  virtual ~TrkEffLinearityCheck(){};
  
 private:


  
  //internal functions
  void LineFitCheck( std::vector<recob::Hit> const&, 
      std::vector<size_t> const&);

  
  
};

#endif
