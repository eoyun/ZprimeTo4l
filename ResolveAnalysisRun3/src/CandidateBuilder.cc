// 무엇: 채널별 pairing 구현. 기존 PairingHelper::pair4E / pairByDR 재현.
// 어떻게: dR pairing 은 3조합 (01)(23)/(02)(13)/(03)(12) 중 min(dR1²,dR2²) 최소 선택.
//   4e reciprocal-GSF: addGsfIdx 가 서로의 ntuple index 를 가리키는 첫 쌍.
// 의존: CandidateBuilder.h, Kinematics.h (deltaR2), <cmath>.
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/CandidateBuilder.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Kinematics.h"

#include <algorithm>

namespace {

using raRun3::CandidateBuilder::PairIdx;

// 4개 점(eta,phi)에 대해 3조합 중 min(dR1²,dR2²) 이 가장 작은 조합을 고른다.
// 기존 PairingHelper::pairByDR 재현.
void chooseByDR(const double eta[4], const double phi[4], PairIdx& p1, PairIdx& p2) {
  // 3가지 조합: {(0,1),(2,3)} / {(0,2),(1,3)} / {(0,3),(1,2)}
  const int combos[3][4] = {
    {0, 1, 2, 3},
    {0, 2, 1, 3},
    {0, 3, 1, 2},
  };

  double bestKey = 1e30;
  int    best    = 0;
  for (int c = 0; c < 3; ++c) {
    const int a1 = combos[c][0], b1 = combos[c][1];
    const int a2 = combos[c][2], b2 = combos[c][3];
    const double dr1 = raRun3::deltaR2(eta[a1], phi[a1], eta[b1], phi[b1]);
    const double dr2 = raRun3::deltaR2(eta[a2], phi[a2], eta[b2], phi[b2]);
    const double key = std::min(dr1, dr2);  // 기존 min(dR1,dR2) 기준
    if (key < bestKey) {
      bestKey = key;
      best = c;
    }
  }

  p1 = {combos[best][0], combos[best][1]};
  p2 = {combos[best][2], combos[best][3]};
}

}  // namespace

namespace raRun3 {
namespace CandidateBuilder {

Candidate buildFourMuon(const std::vector<Muon>& mus) {
  Candidate cand;
  if (mus.size() != 4) return cand;  // valid=false

  // pairing dR 은 muon reco 방향 사용 (기존 pat::Muon::eta()/phi()).
  double eta[4], phi[4];
  for (int i = 0; i < 4; ++i) {
    eta[i] = mus[i].recoEta;
    phi[i] = mus[i].recoPhi;
  }
  chooseByDR(eta, phi, cand.pair1, cand.pair2);
  cand.valid = true;
  return cand;
}

Candidate buildFourElectron(const std::vector<Electron>& eles) {
  Candidate cand;
  if (eles.size() != 4) return cand;

  // 1) reciprocal-GSF 짝 탐색 (기존 scanBoostedElectrons).
  //    addGsfIdx<0 = secondary GSF 없음 → 후보 아님.
  for (int i = 0; i < 4 && !cand.valid; ++i) {
    if (eles[i].addGsfIdx < 0) continue;
    for (int j = i + 1; j < 4; ++j) {
      if (eles[j].addGsfIdx < 0) continue;
      if (eles[i].addGsfIdx == eles[j].index &&
          eles[j].addGsfIdx == eles[i].index) {
        cand.pair1 = {i, j};
        // 나머지 두 위치가 pair2
        int rest[2], n = 0;
        for (int k = 0; k < 4; ++k)
          if (k != i && k != j) rest[n++] = k;
        cand.pair2 = {rest[0], rest[1]};
        cand.valid = true;
        break;
      }
    }
  }
  if (cand.valid) return cand;

  // 2) fallback: closest-dR (electron 은 polarP4 방향 = rawEta/rawPhi).
  double eta[4], phi[4];
  for (int i = 0; i < 4; ++i) {
    eta[i] = eles[i].rawEta;
    phi[i] = eles[i].rawPhi;
  }
  chooseByDR(eta, phi, cand.pair1, cand.pair2);
  cand.valid = true;
  return cand;
}

Candidate buildTwoETwoMu(const std::vector<Electron>& eles, const std::vector<Muon>& mus) {
  Candidate cand;
  if (eles.size() != 2 || mus.size() != 2) return cand;
  cand.pair1 = {0, 1};  // ee (electron vector 위치)
  cand.pair2 = {0, 1};  // μμ (muon vector 위치)
  cand.valid = true;
  return cand;
}

}  // namespace CandidateBuilder
}  // namespace raRun3
