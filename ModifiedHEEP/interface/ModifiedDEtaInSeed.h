#ifndef ModifiedDEtaInSeed_h
#define ModifiedDEtaInSeed_h 1

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

class ModifiedDEtaInSeed {
public:
  ModifiedDEtaInSeed(PositionCalc calc);
  ~ModifiedDEtaInSeed()=default;

  struct variables {
    variables()
    : modifiedDEtaInSeed(std::numeric_limits<float>::max()),
      alphaTrack(std::numeric_limits<float>::max()),
      alphaCalo(std::numeric_limits<float>::max()),
      union5x5ratioDR(std::numeric_limits<float>::max()) {}

    variables(const double alCalo)
    : modifiedDEtaInSeed(std::numeric_limits<float>::max()),
      alphaTrack(std::numeric_limits<float>::max()),
      alphaCalo(alCalo),
      union5x5ratioDR(std::numeric_limits<float>::max()) {}

    variables(const double dEta, const double alTrk, const double alCalo, const double ratioDR)
    : modifiedDEtaInSeed(dEta),
      alphaTrack(alTrk),
      alphaCalo(alCalo),
      union5x5ratioDR(ratioDR) {}

    double modifiedDEtaInSeed;
    double alphaTrack;
    double alphaCalo;
    double union5x5ratioDR;
  };

  variables value(const reco::GsfElectron& aEle,
                  const EcalRecHitCollection* ecalRecHits);

  variables value(const reco::GsfElectron& aEle,
                  const EcalRecHitCollection* ecalRecHits,
                  const reco::TrackBase& addTrk,
                  const reco::BeamSpot& beamSpot,
                  const edm::EventSetup& iSetup);

private:
  PositionCalc posCalcLog_;
};

#endif
