#ifndef ResolveAnalysisRun3_ObjectSelector_h
#define ResolveAnalysisRun3_ObjectSelector_h
// 무엇: muon object 선택 (acceptance/ID/neighbor/modified-iso/P).
// 어떻게: 순수 함수. Muon struct 와 Config 만 입출력. 히스토그램/파일 접근 없음.
// 의존: Muon.h, Config.h. (기존 checkIso / ResolvedMuCRanalyzer 재현)
#include <vector>
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Muon.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Config.h"

namespace raRun3 {
namespace ObjectSelector {

// acceptance: tuneP pt > min, |tuneP eta| < max.
bool passMuonAccept(const Muon& mu, const Config& cfg);

// ID: global high-pT OR tracker high-pT (기존 저장 bool 재사용, 두 ID 유지).
bool passMuonId(const Muon& mu);

// self 와 other(neighbor 후보) 가 기존 checkIso 조건을 만족하는가.
bool isMuonNeighbor(const Muon& self, const Muon& other, const Config& cfg);

// idPool 중 self 제외 neighbor 후보에서 TuneP pt 최고 하나의 innerPt 를
// raw trackIso 에서 뺀 modified iso 값(빼기 전이면 raw 그대로).
double modifiedMuonIso(const Muon& self, const std::vector<Muon>& idPool, const Config& cfg);

// accept+ID 통과 muon 을 idPool 로 삼아, modified iso/tuneP pt < 컷 인 P muon 반환.
std::vector<Muon> selectMuonsP(const std::vector<Muon>& all, const Config& cfg);

}  // namespace ObjectSelector
}  // namespace raRun3
#endif
