#ifndef TRKEFF_TAGCREATORALG_H
#define TRKEFF_TAGCREATORALG_H

#include <string>

#include "TTree.h"

namespace trkeff{
  class TagCreatorAlg;
}

class trkeff::TagCreatorAlg{
  
public:
  
  /// Default constructor
  TagCreatorAlg();

  /// Default destructor
  virtual ~TagCreatorAlg(){};

  void SetupOutputTree(TTree*);
  
 private:

  TTree*      fTree;
  
};

#endif
