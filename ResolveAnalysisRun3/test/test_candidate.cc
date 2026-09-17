#include "ZprimeTo4l/ResolveAnalysisRun3/interface/CandidateBuilder.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"
#include <vector>

using namespace raRun3;
using CandidateBuilder::PairIdx;
using CandidateBuilder::Candidate;

// pair 가 위치 {x,y} 를 (순서 무관) 담는가.
static bool pairIs(const PairIdx& p, int x, int y) {
  return (p.first == x && p.second == y) || (p.first == y && p.second == x);
}

static Muon muAt(double eta, double phi) {
  Muon m; m.recoEta = eta; m.recoPhi = phi; return m;
}

static Electron eleAt(int index, double eta, double phi, int addGsfIdx) {
  Electron e; e.index = index; e.rawEta = eta; e.rawPhi = phi; e.addGsfIdx = addGsfIdx;
  return e;
}

// 4mu: (0,1) 이 가깝고 (2,3) 이 가까운 배치 → 그 조합 선택.
static void test_fourMuon_dr() {
  std::vector<Muon> mus = {
    muAt(0.0, 0.0), muAt(0.0, 0.1),   // 0,1 가까움
    muAt(2.0, 0.0), muAt(2.0, 0.1),   // 2,3 가까움
  };
  Candidate c = CandidateBuilder::buildFourMuon(mus);
  CHECK(c.valid);
  CHECK(pairIs(c.pair1, 0, 1));
  CHECK(pairIs(c.pair2, 2, 3));
}

// 4mu: 입력 개수 오류 → valid=false.
static void test_fourMuon_wrong_count() {
  std::vector<Muon> mus = {muAt(0, 0), muAt(1, 1), muAt(2, 2)};
  CHECK(CandidateBuilder::buildFourMuon(mus).valid == false);
}

// 4e: reciprocal-GSF 가 dR 보다 우선. dR 로는 (0,2)가 가깝지만 reciprocal 은 (0,1).
static void test_fourElectron_reciprocal_priority() {
  std::vector<Electron> eles = {
    eleAt(10, 0.0, 0.0, 11),   // pos0, index10, addGsf→11
    eleAt(11, 2.0, 0.0, 10),   // pos1, index11, addGsf→10  (0,1 reciprocal)
    eleAt(12, 0.0, 0.05, -1),  // pos2, index12, dR 로는 pos0 과 매우 가까움
    eleAt(13, 2.0, 0.05, -1),  // pos3, index13
  };
  Candidate c = CandidateBuilder::buildFourElectron(eles);
  CHECK(c.valid);
  CHECK(pairIs(c.pair1, 0, 1));  // reciprocal 강제
  CHECK(pairIs(c.pair2, 2, 3));
}

// 4e: reciprocal 없음 → dR fallback.
static void test_fourElectron_dr_fallback() {
  std::vector<Electron> eles = {
    eleAt(10, 0.0, 0.0, -1), eleAt(11, 0.0, 0.1, -1),  // 0,1 가까움
    eleAt(12, 3.0, 0.0, -1), eleAt(13, 3.0, 0.1, -1),  // 2,3 가까움
  };
  Candidate c = CandidateBuilder::buildFourElectron(eles);
  CHECK(c.valid);
  CHECK(pairIs(c.pair1, 0, 1));
  CHECK(pairIs(c.pair2, 2, 3));
}

// 2e2mu: pair1=ee(0,1), pair2=μμ(0,1).
static void test_twoETwoMu() {
  std::vector<Electron> eles = {eleAt(10, 0, 0, -1), eleAt(11, 1, 0, -1)};
  std::vector<Muon> mus = {muAt(0, 0), muAt(1, 0)};
  Candidate c = CandidateBuilder::buildTwoETwoMu(eles, mus);
  CHECK(c.valid);
  CHECK(pairIs(c.pair1, 0, 1));
  CHECK(pairIs(c.pair2, 0, 1));
}

int main() {
  RUN(test_fourMuon_dr);
  RUN(test_fourMuon_wrong_count);
  RUN(test_fourElectron_reciprocal_priority);
  RUN(test_fourElectron_dr_fallback);
  RUN(test_twoETwoMu);
  REPORT();
}
