/**
 * \file TrkEffAlg.h
 *
 * 
 * \brief Little sample program for establishing a user analysis space.
 *
 * @author wketchum
*/

#ifndef TEST_USERANALYSIS_H
#define TEST_USERANALYSIS_H

#include <string>

#include "TTree.h"

namespace trkeff{
  class TrkEffAlg;
}

class trkeff::TrkEffAlg{
  
public:
  
  /// Default constructor
  TrkEffAlg();

  /// Default destructor
  virtual ~TrkEffAlg(){};

  void RunAnalysis();
  void SetupOutputTree(TTree*);
  
 private:

  std::string fAlgName;
  TTree*      fTree;
  
  void PrintInfo();

  
};

#endif
