#ifndef ResolveAnalysisRun3_CandidateBuilder_h
#define ResolveAnalysisRun3_CandidateBuilder_h
// 무엇: 채널별 lepton pairing (4e / 4mu / 2e2mu). 기존 PairingHelper 재현.
// 어떻게: 순수 함수. 입력 vector 안의 "위치 index(0..n-1)"로 두 pair 를 돌려준다.
//   - 4mu   : pair1,pair2 둘 다 muon vector 위치
//   - 4e    : pair1,pair2 둘 다 electron vector 위치 (reciprocal-GSF 우선, 실패시 dR)
//   - 2e2mu : pair1 = electron vector 위치(0,1), pair2 = muon vector 위치(0,1)
// 의존: Muon.h, Electron.h, Kinematics.h (deltaR2).
#include <vector>
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Muon.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Electron.h"

namespace raRun3 {
namespace CandidateBuilder {

// 한 pair = 입력 vector 안의 두 위치.
struct PairIdx {
  int first  = -1;
  int second = -1;
};

struct Candidate {
  PairIdx pair1;
  PairIdx pair2;
  bool    valid = false;  // 입력 개수가 맞지 않으면 false
};

// 4mu: 정확히 4개. closest-dR (min(dR1²,dR2²) 최소 조합).
Candidate buildFourMuon(const std::vector<Muon>& mus);

// 4e: 정확히 4개. reciprocal-GSF 짝 우선, 없으면 closest-dR.
Candidate buildFourElectron(const std::vector<Electron>& eles);

// 2e2mu: electron 2 + muon 2. ee 1쌍 + μμ 1쌍 (조합 유일).
Candidate buildTwoETwoMu(const std::vector<Electron>& eles, const std::vector<Muon>& mus);

}  // namespace CandidateBuilder
}  // namespace raRun3
#endif
