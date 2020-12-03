//*****************************************************************************
// File:      ModifiedTowerIsolation.cc
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer
// Institute: IIHE-VUB
//=============================================================================
//*****************************************************************************

//CMSSW includes
#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedTowerIsolation.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include <cassert>
#include <memory>

#ifdef ETISTATDEBUG
// #include<iostream>
namespace etiStat {
  Count::~Count() {
    //    std::cout << "\nModifiedTowerIsolationNew " << create << "/" << comp << "/" << float(span)/float(comp)
    //	      << std::endl<< std::endl;
    }

  Count Count::count;
}
#endif

namespace {
  struct TLS {
    std::unique_ptr<ModifiedTowerIsolationNew<1>> newAlgo=nullptr;;
    const CaloTowerCollection* oldTowers=nullptr;;
    uint32_t id15=0;
  };
  thread_local TLS tls;
}

ModifiedTowerIsolation::ModifiedTowerIsolation (float extRadiusI,
					    float intRadiusI,
					    float etLow,
					    signed int depth,
					    const CaloTowerCollection* towers ) :
  depth_(depth),
  extRadius(extRadiusI),
  intRadius(intRadiusI)
{
  assert(0==etLow);

  // extremely poor in quality  (test of performance)
  if (tls.newAlgo.get()==nullptr ||  towers!=tls.oldTowers || towers->size()!=tls.newAlgo->nt || (towers->size()>15 && (*towers)[15].id()!=tls.id15)) {
    tls.newAlgo = std::make_unique<ModifiedTowerIsolationNew<1>>(&extRadius,&intRadius,*towers);
    tls.oldTowers=towers;
    tls.id15 = towers->size()>15 ? (*towers)[15].id() : 0;
  }
}


double  ModifiedTowerIsolation::getSum (bool et, reco::SuperCluster const & sc, const reco::TrackBase& addTrk, const std::vector<CaloTowerDetId> * detIdToExclude) const{

  if (nullptr!=detIdToExclude) assert(0==intRadius);

  // hack
  tls.newAlgo->setRadius(&extRadius,&intRadius);

  ModifiedTowerIsolationNew<1>::Sum sum;
  tls.newAlgo->compute(et, sum, sc, addTrk,
		  (detIdToExclude==nullptr) ? nullptr : &((*detIdToExclude).front()),
		  (detIdToExclude==nullptr) ? nullptr : (&(*detIdToExclude).back())+1
		  );

  switch(depth_){
  case AllDepths: return detIdToExclude==nullptr ? sum.he[0] : sum.heBC[0];
  case Depth1: return detIdToExclude==nullptr ? sum.he[0]-sum.h2[0] : sum.heBC[0]-sum.h2BC[0];
  case Depth2:return detIdToExclude==nullptr ? sum.h2[0] : sum.h2BC[0];
  default: return 0;
  }
  return 0;
}
