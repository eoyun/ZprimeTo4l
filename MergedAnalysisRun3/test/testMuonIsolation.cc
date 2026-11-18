#include "ZprimeTo4l/MergedAnalysisRun3/interface/MuonIsolation.h"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cmath>
#include <iostream>

int main() {
  using namespace mergedRun3;
  std::vector<IsolationMuon> mu(3);
  mu[0].isoValid = true; mu[0].rawIso = 12.;
  // Highest TuneP neighbor fails the resolved pool, but is kept for study isolation.
  mu[1].tracksValid = true; mu[1].innerEta = 0.1;
  mu[1].tunePt = 80.; mu[1].innerPt = 18.;
  mu[2].tracksValid = true; mu[2].innerEta = 0.2;
  mu[2].tunePt = 40.; mu[2].innerPt = 7.; mu[2].inResolvedPool = true;
  const auto study = modifiedIsolation(mu, 0, true, false);
  const auto resolved = modifiedIsolation(mu, 0, true, true);
  assert(study.valid && study.neighborIndex == 1 && study.sumPt == -6.);
  assert(resolved.valid && resolved.neighborIndex == 2 && resolved.sumPt == 5.);
  assert(!modifiedIsolation(mu, 0, false, false).valid);
  mu[0].isoValid = false;
  assert(!modifiedIsolation(mu, 0, true, true).valid);
  mu[0].isoValid = true;
  // Signed dxy and inclusive boundaries reproduce legacy checkIso.
  mu.resize(2);
  mu[1].innerDxyBS = -2.; mu[1].innerEta = 0.3; mu[1].innerVz = 0.2;
  assert(modifiedIsolation(mu, 0, true, false).neighborIndex == 1);
  mu[1].innerDxyBS = 0.10001;
  assert(modifiedIsolation(mu, 0, true, false).sumPt == 12.);
  mu[1].innerDxyBS = 0.1; mu[1].innerEta = 0.01;
  assert(modifiedIsolation(mu, 0, true, false).neighborIndex == 1);
  mu[1].innerEta = 0.009;
  assert(modifiedIsolation(mu, 0, true, false).neighborIndex == -1);
  // Angular wrapping across +/-pi and explicit self exclusion.
  mu[0].phi = std::acos(-1.) - 0.02; mu[1].innerPhi = -std::acos(-1.) + 0.02;
  mu[1].innerEta = 0.;
  assert(modifiedIsolation(mu, 0, true, false).neighborIndex == 1);
  mu.resize(1); mu[0].tracksValid = true; mu[0].innerEta = 0.1; mu[0].innerPhi = mu[0].phi;
  mu[0].tunePt = 100.; mu[0].innerPt = 90.;
  assert(modifiedIsolation(mu, 0, true, false).sumPt == 12.);
  std::cout << "Muon isolation regression checks passed\n";
}
