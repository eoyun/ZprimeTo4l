#ifndef PairingHelper_h
#define PairingHelper_h 1

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

namespace PairingHelper {
  bool pairBoostedElectrons(const pat::ElectronRef& a,
                            const pat::ElectronRef& b,
                            const edm::Handle<edm::ValueMap<reco::GsfTrackRef>>& addGsfTrkMap);
  bool scanBoostedElectrons(const std::vector<pat::ElectronRef>& eles,
                            const edm::Handle<edm::ValueMap<reco::GsfTrackRef>>& addGsfTrkMap,
                            std::pair<pat::ElectronRef,pat::ElectronRef>& pair);
  bool pair4E(const std::vector<pat::ElectronRef>& eles,
              const edm::Handle<edm::ValueMap<reco::GsfTrackRef>>& addGsfTrkMap,
              std::pair<pat::ElectronRef,pat::ElectronRef>& pair1,
              std::pair<pat::ElectronRef,pat::ElectronRef>& pair2);

  bool pairByInvM(const std::vector<pat::ElectronRef>& eles,
                  std::pair<pat::ElectronRef,pat::ElectronRef>& pair1,
                  std::pair<pat::ElectronRef,pat::ElectronRef>& pair2);

  bool pairByDR(const std::vector<pat::ElectronRef>& eles,
                std::pair<pat::ElectronRef,pat::ElectronRef>& pair1,
                std::pair<pat::ElectronRef,pat::ElectronRef>& pair2);
};

#endif
