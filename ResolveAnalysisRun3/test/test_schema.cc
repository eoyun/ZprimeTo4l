#include "ZprimeTo4l/ResolveAnalysisRun3/interface/NtupleSchema.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"
#include <string>
#include <set>

using namespace raRun3::schema;

// 이름 규칙(prefix) + 버전 확인.
static void test_names_prefix_and_version() {
  CHECK(std::string(mu::corrTunePpt).rfind("mu_", 0) == 0);
  CHECK(std::string(ele::selEt).rfind("ele_", 0) == 0);
  CHECK(std::string(trig::pt).rfind("trigObj_", 0) == 0);
  CHECK(kVersion >= 1);
}

// 대표 branch 이름들이 서로 유일한지(오타로 겹치지 않는지) 확인.
static void test_names_unique() {
  std::set<std::string> names = {
    ev::run, ev::lumi, ev::event, ev::genWeight, ev::puTrue,
    mu::index, mu::charge, mu::corrTunePpt, mu::rawTunePpt, mu::tunePeta, mu::tunePphi,
    mu::isHighPt, mu::isTrackerHighPt, mu::isTrackerMuon, mu::trackIso, mu::innerPt,
    mu::recoEta, mu::recoPhi, mu::recoVz, mu::innerEta, mu::innerPhi, mu::innerVz,
    mu::innerDxyBS, mu::relPtErr, mu::trkLayers, mu::pixelHits, mu::matchedStations,
    mu::dxy, mu::dz,
    ele::index, ele::charge, ele::selEt, ele::etaSC, ele::rawPt, ele::rawEta, ele::rawPhi,
    ele::rawEnergy, ele::corrPt, ele::corrEta, ele::corrPhi, ele::corrM,
    ele::passModHeep, ele::modHeepBitmap, ele::addGsfIdx,
    trig::pt, trig::eta, trig::phi, trig::filterBits,
  };
  // 위 목록에 넣은 항목 수 (중복 있으면 set 크기가 작아짐)
  CHECK(names.size() == 48);
}

int main() {
  RUN(test_names_prefix_and_version);
  RUN(test_names_unique);
  REPORT();
}
