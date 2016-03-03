////////////////////////////////////////////////////////////////////////
// Class:       TrkEffTagCreator
// Module Type: producer
// File:        TrkEffTagCreator_module.cc
//
// Generated at Thu Jan 28 07:50:04 2016 by Wesley Ketchum using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include <memory>

#include "RecoBase/Hit.h"
#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"

#include "TrkEffTagAlgs/TagCreatorAlg.hh"
#include "TrkEffTagObjects/TrkEffTag.h"

namespace trkeff {
  class TrkEffTagCreator;
}

class trkeff::TrkEffTagCreator : public art::EDProducer {
public:
  explicit TrkEffTagCreator(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TrkEffTagCreator(TrkEffTagCreator const &) = delete;
  TrkEffTagCreator(TrkEffTagCreator &&) = delete;
  TrkEffTagCreator & operator = (TrkEffTagCreator const &) = delete;
  TrkEffTagCreator & operator = (TrkEffTagCreator &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  std::string   fHitCollectionLabel;
  TagCreatorAlg fTagCreator;

};


trkeff::TrkEffTagCreator::TrkEffTagCreator(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  produces< std::vector<trkeff::TrkEffTag> >();
  this->reconfigure(p);

  art::ServiceHandle<art::TFileService> tfs;
  fTagCreator.SetupOutputTree(tfs->make<TTree>("tag_tree","Tag tree"),
			      tfs->make<TTree>("event_tree","Event tree"));
			    
}

void trkeff::TrkEffTagCreator::produce(art::Event & e)
{

  art::Handle< std::vector<recob::Hit> > hitHandle;
  e.getByLabel(fHitCollectionLabel,hitHandle);

  std::unique_ptr< std::vector<trkeff::TrkEffTag> > tagHandle(new std::vector<trkeff::TrkEffTag>);
  
  art::ServiceHandle<geo::Geometry> geoHandle;
  art::ServiceHandle<util::DetectorProperties> detpHandle;
  art::ServiceHandle<util::LArProperties> larpHandle;

  fTagCreator.CreateTags(*hitHandle,
			 *tagHandle,
			 *geoHandle,
			 *detpHandle,
			 *larpHandle,
			 e.run(),
			 e.event());

  e.put(std::move(tagHandle));
  
}

void trkeff::TrkEffTagCreator::reconfigure(fhicl::ParameterSet const & p)
{
  art::ServiceHandle<geo::Geometry> geoHandle;

  fHitCollectionLabel = p.get<std::string>("HitCollectionLabel");
  fTagCreator.Configure(p.get<fhicl::ParameterSet>("TagCreatorAlgParams"),
			*geoHandle);
}

DEFINE_ART_MODULE(trkeff::TrkEffTagCreator)
