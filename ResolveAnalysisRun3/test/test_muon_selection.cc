#include "ZprimeTo4l/ResolveAnalysisRun3/interface/ObjectSelector.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Config.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"
#include <vector>

using namespace raRun3;

// 테스트용 config: neighbor 창 0.01<dR<0.3, dz<0.2, dxyBS<0.1, iso/pt<0.1.
static Config makeCfg() {
  Config c;
  c.muon.tunePptMin = 20.0;
  c.muon.etaMax = 2.4;
  c.muon.modIsoRelMax = 0.1;
  c.muon.neighbor = {0.01, 0.3, 0.2, 0.1};
  c.massCuts = {1.0, 200.0};
  return c;
}

// 방향이 (eta,phi)=(0,0) 근처인 기본 muon 하나 (ID 통과, 격리됨).
static Muon baseMuon() {
  Muon m;
  m.charge = 1;
  m.corrTunePpt = 60.0; m.rawTunePpt = 60.0; m.tunePeta = 0.0; m.tunePphi = 0.0;
  m.isHighPt = true; m.isTrackerHighPt = false;
  m.trackIso = 1.0; m.innerPt = 59.0;
  m.recoEta = 0.0; m.recoPhi = 0.0; m.recoVz = 0.0;
  m.innerEta = 0.0; m.innerPhi = 0.0; m.innerVz = 0.0; m.innerDxyBS = 0.0;
  return m;
}

// 1) neighbor 판정: dR 창 안/밖.
static void test_isNeighbor_dr_window() {
  Config cfg = makeCfg();
  Muon self = baseMuon();
  Muon near = baseMuon();  near.innerPhi = 0.1;   // dR=0.1 → 창 안
  Muon far  = baseMuon();  far.innerPhi = 0.5;    // dR=0.5 → 창 밖
  Muon coincident = baseMuon(); coincident.innerPhi = 0.0;  // dR=0 → drMin 미만, 배제
  CHECK(ObjectSelector::isMuonNeighbor(self, near, cfg) == true);
  CHECK(ObjectSelector::isMuonNeighbor(self, far, cfg) == false);
  CHECK(ObjectSelector::isMuonNeighbor(self, coincident, cfg) == false);
}

// 2) neighbor 판정: dz, dxyBS 컷.
static void test_isNeighbor_dz_dxy() {
  Config cfg = makeCfg();
  Muon self = baseMuon();
  Muon farZ = baseMuon(); farZ.innerPhi = 0.1; farZ.innerVz = 0.5;      // |dz|=0.5>0.2
  Muon bigDxy = baseMuon(); bigDxy.innerPhi = 0.1; bigDxy.innerDxyBS = 0.2; // >0.1
  CHECK(ObjectSelector::isMuonNeighbor(self, farZ, cfg) == false);
  CHECK(ObjectSelector::isMuonNeighbor(self, bigDxy, cfg) == false);
}

// 3) modified iso: neighbor 있으면 가장 높은 TuneP 후보의 innerPt 를 뺀다.
static void test_modifiedIso_subtracts_highest() {
  Config cfg = makeCfg();
  Muon self = baseMuon(); self.trackIso = 30.0;  // raw iso 큼
  Muon nb1 = baseMuon(); nb1.innerPhi = 0.1; nb1.corrTunePpt = 40.0; nb1.innerPt = 25.0;
  Muon nb2 = baseMuon(); nb2.innerPhi = 0.1; nb2.corrTunePpt = 55.0; nb2.innerPt = 28.0; // TuneP 최고
  std::vector<Muon> pool = {self, nb1, nb2};
  // 가장 높은 TuneP(nb2)의 innerPt=28 를 뺌 → 30-28=2
  CHECK_CLOSE(ObjectSelector::modifiedMuonIso(self, pool, cfg), 2.0, 1e-9);
}

// 4) neighbor 없으면 raw iso 그대로.
static void test_modifiedIso_no_neighbor() {
  Config cfg = makeCfg();
  Muon self = baseMuon(); self.trackIso = 5.0;
  Muon farOne = baseMuon(); farOne.innerPhi = 1.0;  // 창 밖
  std::vector<Muon> pool = {self, farOne};
  CHECK_CLOSE(ObjectSelector::modifiedMuonIso(self, pool, cfg), 5.0, 1e-9);
}

// 5) P 선택: ID 통과 + acceptance + modified iso/pt < 컷.
static void test_selectMuonsP() {
  Config cfg = makeCfg();
  Muon good = baseMuon();                 // iso 1/60 < 0.1 → P
  Muon lowPt = baseMuon(); lowPt.corrTunePpt = 10.0; lowPt.innerPt = 9.0; // pt<20 → 탈락
  Muon noId = baseMuon(); noId.isHighPt = false; noId.isTrackerHighPt = false; // ID 실패
  Muon notIso = baseMuon(); notIso.trackIso = 30.0; // 30/60=0.5>0.1, neighbor 없음 → 탈락
  std::vector<Muon> all = {good, lowPt, noId, notIso};
  std::vector<Muon> passP = ObjectSelector::selectMuonsP(all, cfg);
  CHECK(passP.size() == 1);
}

int main() {
  RUN(test_isNeighbor_dr_window);
  RUN(test_isNeighbor_dz_dxy);
  RUN(test_modifiedIso_subtracts_highest);
  RUN(test_modifiedIso_no_neighbor);
  RUN(test_selectMuonsP);
  REPORT();
}
