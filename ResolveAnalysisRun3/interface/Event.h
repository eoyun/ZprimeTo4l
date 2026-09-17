#ifndef ResolveAnalysisRun3_Event_h
#define ResolveAnalysisRun3_Event_h
// 무엇: 후단 event 하나. 스칼라 + muon/electron 가변길이 배열.
// 어떻게: 순수 데이터. NtupleReader(Plan 3)가 채운다.
// 의존: Muon.h, Electron.h.
#include <vector>
#include <cstdint>
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Muon.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Electron.h"

namespace raRun3 {

struct Event {
  unsigned int run  = 0;
  unsigned int lumi = 0;
  uint64_t     event = 0;
  double       genWeight = 1.;  // MC weight (data=1)
  std::vector<Muon>     muons;
  std::vector<Electron> electrons;
};

}  // namespace raRun3
#endif
