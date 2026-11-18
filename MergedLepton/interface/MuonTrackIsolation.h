#ifndef ZprimeTo4l_MergedLepton_MuonTrackIsolation_h
#define ZprimeTo4l_MergedLepton_MuonTrackIsolation_h

#include <cmath>
#include <vector>
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"

namespace mergedMuonIso {
// Shared with Analysis/MuonCorrectionHelper::checkIso. Preserve its exact
// cone, vertex, and SIGNED dxy conditions (including boundary acceptance).
inline bool checkIso(const reco::Muon& mu, const reco::Track& trk,
                     const reco::BeamSpot& bs) {
  const double dr2 = reco::deltaR2(mu.eta(), mu.phi(), trk.eta(), trk.phi());
  if (dr2 > 0.09 || dr2 < 0.0001) return false;
  if (std::abs(mu.vz() - trk.vz()) > 0.2) return false;
  if (trk.dxy(bs.position()) > 0.1) return false;
  return true;
}

// Absolute modified R03 isolation (GeV), computed once per input muon.
// As in the Analysis boosted-isolation loops, subtract only ONE candidate:
// highest TuneP pT among tracks passing checkIso, using its inner-track pT.
// Unlike the Analysis event selection, no pT/HighPt-ID cuts are imposed on
// companions here: these ntuples also serve low-pT J/psi studies.
// Invalid isolation or missing beam spot -> -999; no companion -> raw iso.
// Negative corrected values are retained, matching the Analysis convention.
template <typename MuonCollection>
std::vector<float> modifiedTrackIsoValues(const MuonCollection& muons,
                                          const reco::BeamSpot* bs) {
  std::vector<float> values(muons.size(), -999.f);
  if (!bs) return values;
  for (unsigned i = 0; i < muons.size(); ++i) {
    const auto& mu = muons[i];
    if (!mu.isIsolationValid()) continue;
    values[i] = mu.isolationR03().sumPt;
    double highestPt = -1.;
    double subtractPt = 0.;
    for (unsigned j = 0; j < muons.size(); ++j) {
      if (i == j) continue;
      const auto& other = muons[j];
      const auto inner = other.innerTrack();
      const auto tuneP = other.tunePMuonBestTrack();
      if (inner.isNull() || !inner.isAvailable() ||
          tuneP.isNull() || !tuneP.isAvailable()) continue;
      if (!checkIso(mu, *inner, *bs)) continue;
      if (tuneP->pt() > highestPt) {
        highestPt = tuneP->pt();
        subtractPt = inner->pt();
      }
    }
    values[i] -= subtractPt;
  }
  return values;
}
}  // namespace mergedMuonIso
#endif
