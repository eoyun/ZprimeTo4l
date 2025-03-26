#ifndef ModifiedDEtaInSeed_h
#define ModifiedDEtaInSeed_h 1

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronUtilities.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

class ModifiedDEtaInSeed {
public:
  ModifiedDEtaInSeed(PositionCalc calc, edm::ConsumesCollector iC);
  ModifiedDEtaInSeed(PositionCalc calc);
  ModifiedDEtaInSeed();
  ~ModifiedDEtaInSeed()=default;

  struct variables {
    variables()
    : dPerpIn(std::numeric_limits<float>::max()),
      dEtaInSeed2nd(std::numeric_limits<float>::max()),
      dPhiInSC2nd(std::numeric_limits<float>::max()),
      alphaTrack(std::numeric_limits<float>::max()),
      normalizedDParaIn(std::numeric_limits<float>::max()) {}

    variables(const double dPerp, const double dEta2nd, const double dPhi2nd,
              const double alTrk, const double dPara)
    : dPerpIn(dPerp),
      dEtaInSeed2nd(dEta2nd),
      dPhiInSC2nd(dPhi2nd),
      alphaTrack(alTrk),
      normalizedDParaIn(dPara) {}

    double dPerpIn;
    double dEtaInSeed2nd;
    double dPhiInSC2nd;
    double alphaTrack;
    double normalizedDParaIn;
  };

  variables value(const reco::GsfElectron& aEle,
                  const EcalRecHitCollection* ecalRecHits,
                  const reco::TrackBase& addTrk,
                  const reco::BeamSpot& beamSpot,
                  const edm::EventSetup& iSetup);

  bool extrapolate(const reco::GsfElectron& aEle, const reco::TrackBase& addTrk,
                   const math::XYZPoint& beamSpotPos, const edm::EventSetup& iSetup,
                   EleRelPointPair& scAtVtx, EleRelPointPair& seedAtCalo);

private:
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken_;
  edm::ESGetToken<CaloTopology, CaloTopologyRecord> topologyToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticToken_;
  edm::ESGetToken<GeometricSearchTracker, TrackerRecoGeometryRecord> geotrkToken_;
  
  PositionCalc posCalcLog_;
};

#endif
