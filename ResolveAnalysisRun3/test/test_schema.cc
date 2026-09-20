#include "ZprimeTo4l/ResolveAnalysisRun3/interface/NtupleSchema.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"
#include <string>
#include <vector>
#include <set>

using namespace raRun3::schema;

// flat 이름 규칙(prefix) + 버전 확인.
static void test_names_prefix_and_version() {
  CHECK(std::string(mu::pt).rfind("muon_", 0) == 0);
  CHECK(std::string(ele::pt).rfind("electron_", 0) == 0);
  CHECK(std::string(ele::modTrkIso).rfind("electron_", 0) == 0);
  CHECK(std::string(trig::pt).rfind("trigObj_", 0) == 0);
  CHECK(kVersion >= 2);
}

// 모든 branch 이름을 모아 중복이 없는지(set 크기 == 나열 개수) 확인. 매직넘버 없이 무결성만.
static void test_names_unique() {
  std::vector<std::string> names = {
    ev::run, ev::lumi, ev::event, ev::genWeight, ev::puTrue, ev::nPV, ev::rho,
    ev::hltFired, ev::passMETfilters,
    mu::index, mu::charge, mu::pt, mu::ptRaw, mu::eta, mu::phi, mu::recoEta, mu::recoPhi,
    mu::recoVz, mu::innerPt, mu::innerEta, mu::innerPhi, mu::innerVz, mu::innerDxyBS,
    mu::isHighPt, mu::isTrackerHighPt, mu::isTrackerMuon, mu::trackIso, mu::relPtErr,
    mu::trkLayers, mu::pixelHits, mu::matchedStations, mu::dxy, mu::dz,
    ele::index, ele::charge, ele::et, ele::ptRaw, ele::etaRaw, ele::phiRaw, ele::energyRaw,
    ele::pt, ele::eta, ele::phi, ele::mass, ele::energyCorr, ele::etaSC,
    ele::passModHeep, ele::modHeepBitmap, ele::addGsfIdx,
    ele::hOverE, ele::sigmaIeta, ele::dEtaInSeed, ele::dPhiIn, ele::e2x5, ele::e5x5,
    ele::dr03TkSumPt, ele::dr03EcalIso, ele::dr03HcalIso, ele::missingHits, ele::dxy, ele::dz,
    ele::ecalDriven,
    ele::modTrkIso, ele::modEcalHcalIso, ele::union5x5covIeIe, ele::union5x5covIeIp,
    ele::union5x5covIpIp, ele::union5x5dEtaIn, ele::union5x5dPhiIn, ele::union5x5Energy,
    ele::dEtaInSeed2nd, ele::dPhiInSC2nd, ele::dPerpIn, ele::alphaTrack, ele::alphaCalo,
    ele::normDParaIn,
    trig::pt, trig::eta, trig::phi, trig::filterBits,
  };
  std::set<std::string> uniq(names.begin(), names.end());
  CHECK(uniq.size() == names.size());  // 중복(오타) 없으면 동일
}

int main() {
  RUN(test_names_prefix_and_version);
  RUN(test_names_unique);
  REPORT();
}
