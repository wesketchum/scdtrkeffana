////////////////////////////////////////////////////////////////////////
// Class:       TrkEffTagAnalyzer
// Module Type: analyzer
// File:        TrkEffTagAnalyzer_module.cc
//
// Generated at Tue Nov 10 13:06:09 2015 by Wesley Ketchum using artmod
// from cetpkgsupport v1_08_07.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "TTree.h"

#include "TrkEffTagAlgs/TagCompareAlg.hh"

#include <string>
#include <vector>

namespace trkeff {
  class TrkEffTagAnalyzer;
}

class trkeff::TrkEffTagAnalyzer : public art::EDAnalyzer {
public:
  explicit TrkEffTagAnalyzer(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrkEffTagAnalyzer(TrkEffTagAnalyzer const &) = delete;
  TrkEffTagAnalyzer(TrkEffTagAnalyzer &&) = delete;
  TrkEffTagAnalyzer & operator = (TrkEffTagAnalyzer const &) = delete;
  TrkEffTagAnalyzer & operator = (TrkEffTagAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  // Declare member data here.
  trkeff::TagCompareAlg fAlg;
  std::string fTagCollectionLabel;
  std::string fMCTrkCollectionLabel;  

  typedef struct EventSum{
    unsigned int run;
    unsigned int event;
    unsigned int n_tags;
    unsigned int n_muons;
    double muon_length;
  } EventSum_t;
  EventSum_t fEventSum;
  TTree *fTreeEventSum;
};


trkeff::TrkEffTagAnalyzer::TrkEffTagAnalyzer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  this->reconfigure(p);
  art::ServiceHandle<art::TFileService> tfs;
  fAlg.SetupOutputTree(tfs->make<TTree>("tagcmpana","TagComparison Tree"));

  fTreeEventSum = tfs->make<TTree>("evsum","Event Summary");
  fTreeEventSum->Branch("sum",&fEventSum,"run/i:event/i:n_tags/i:n_muons/i:muon_length/D");

}

void trkeff::TrkEffTagAnalyzer::analyze(art::Event const & e)
{

  art::Handle< std::vector<trkeff::TrkEffTag> > tagHandle;
  e.getByLabel(fTagCollectionLabel,tagHandle);
  auto const& tagVector(*tagHandle);

  art::Handle< std::vector<sim::MCTrack> > mctrkHandle;
  e.getByLabel(fMCTrkCollectionLabel,mctrkHandle);
  auto const& mctrkVector(*mctrkHandle);

  fEventSum.run = e.run();
  fEventSum.event = e.event();
  fEventSum.n_tags = tagVector.size();

  fEventSum.n_muons=0;  size_t i_p=99;
  for(size_t i=0; i<mctrkVector.size(); ++i)
    if(mctrkVector[i].PdgCode()==13) {
      i_p = i; ++fEventSum.n_muons;
      fEventSum.muon_length = (mctrkVector[i_p].End().Position() - mctrkVector[i_p].Start().Position()).Vect().Mag();
    }
  
  fTreeEventSum->Fill();
  
  
  if(tagVector.size()==0 || mctrkVector.size()==0)
    return;
  
  if(tagVector.size()!=1)
    throw cet::exception("TrkEffTagAnalyzer")
      << "Tag vector has size greater than 1, which was unexpected...";

  if(fEventSum.n_muons!=1){
      throw cet::exception("TrkEffTagAnalyzer")
	<< "MCTrack vector size greater than 1, which was unexpected...";
  }
  
  fAlg.Compare(tagVector[0],mctrkVector[i_p]);
}

void trkeff::TrkEffTagAnalyzer::reconfigure(fhicl::ParameterSet const & p)
{
  fTagCollectionLabel = p.get<std::string>("TagCollectionLabel");
  fMCTrkCollectionLabel = p.get<std::string>("MCTrackCollectionLabel");
}

DEFINE_ART_MODULE(trkeff::TrkEffTagAnalyzer)