#ifndef ZprimeTo4l_MergedAnalysisRun3_MuonIsolation_h
#define ZprimeTo4l_MergedAnalysisRun3_MuonIsolation_h

#include <cmath>
#include <cstddef>
#include <vector>

namespace mergedRun3 {
// Inputs retain the ORIGINAL collection ordering, including unstored companions.
struct IsolationMuon {
  double eta = 0., phi = 0., vz = 0., rawIso = 0.;
  double innerEta = 0., innerPhi = 0., innerVz = 0., innerDxyBS = 0.;
  double innerPt = 0., tunePt = 0.;
  bool isoValid = false, tracksValid = false, inResolvedPool = false;
};
struct IsolationResult {
  double sumPt = -999.;
  int neighborIndex = -1;
  bool valid = false;
};

// Same cone, inclusive boundaries and SIGNED dxy as Run-2 checkIso.
// resolvedPool=true: acceptance + HighPt/TrackerHighPt pool from ResolvedMuCRanalyzer.
// false: all track-valid companions, matching MergedMuon/MuonJpsi study inputs.
inline IsolationResult modifiedIsolation(const std::vector<IsolationMuon>& muons,
                                         std::size_t self, bool beamSpotValid,
                                         bool resolvedPool) {
  const auto& mu = muons.at(self);
  if (!beamSpotValid || !mu.isoValid) return {};
  IsolationResult result{mu.rawIso, -1, true};
  double highestPt = -1.;
  for (std::size_t j = 0; j < muons.size(); ++j) {
    const auto& other = muons[j];
    if (j == self || !other.tracksValid || (resolvedPool && !other.inResolvedPool)) continue;
    const double dphi = std::remainder(mu.phi - other.innerPhi, 2. * std::acos(-1.));
    const double deta = mu.eta - other.innerEta;
    const double dr2 = deta * deta + dphi * dphi;
    if (dr2 > 0.09 || dr2 < 0.0001 || std::abs(mu.vz - other.innerVz) > 0.2 ||
        other.innerDxyBS > 0.1) continue;
    if (other.tunePt > highestPt) {
      highestPt = other.tunePt;
      result.neighborIndex = static_cast<int>(j);
    }
  }
  if (result.neighborIndex >= 0) result.sumPt -= muons[result.neighborIndex].innerPt;
  return result;  // Negative isolation is intentionally retained.
}
}  // namespace mergedRun3
#endif
