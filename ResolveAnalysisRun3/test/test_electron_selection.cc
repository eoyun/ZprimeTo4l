#include "ZprimeTo4l/ResolveAnalysisRun3/interface/ObjectSelector.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Config.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"
#include <vector>

using namespace raRun3;

// 테스트용 config: EB-EE gap 1.4442~1.566, |etaSC|<2.5, HEEP 마스크 0x7B0/0xFFF.
static Config makeCfg() {
  Config c;
  c.electron.gapLo = 1.4442;
  c.electron.gapHi = 1.566;
  c.electron.eeEtaMax = 2.5;
  c.electron.heepMaskLoose = 0x7B0;  // 1968
  c.electron.heepAllPass   = 0xFFF;  // 4095
  return c;
}

// 기본 electron: EB 안(|etaSC|=0.5), modified-HEEP 통과.
static Electron baseEle() {
  Electron e;
  e.charge = -1;
  e.selEt = 50.0;
  e.etaSC = 0.5;
  e.passModHeep = true;
  e.modHeepBitmap = 0xFFF;  // 전부 통과
  return e;
}

// P: accept + passModHeep. gap/범위밖/HEEP실패는 P 아님.
static void test_selectElectronsP() {
  Config cfg = makeCfg();
  Electron good = baseEle();                        // P
  Electron gap  = baseEle(); gap.etaSC = 1.5;       // EB-EE gap → 탈락
  Electron oob  = baseEle(); oob.etaSC = 2.6;       // |etaSC|>2.5 → 탈락
  Electron noId = baseEle(); noId.passModHeep = false;  // VID 실패 → P 아님
  std::vector<Electron> all = {good, gap, oob, noId};
  std::vector<Electron> passP = ObjectSelector::selectElectronsP(all, cfg);
  CHECK(passP.size() == 1);
}

// acceptance 단독 확인: EE 영역(1.566~2.5)은 통과, gap 은 탈락.
static void test_accept_ee_and_gap() {
  Config cfg = makeCfg();
  Electron ee = baseEle(); ee.etaSC = 2.0;   // EE 영역 → accept
  Electron gap = baseEle(); gap.etaSC = 1.5; // gap → reject
  CHECK(ObjectSelector::passElectronAccept(ee, cfg) == true);
  CHECK(ObjectSelector::passElectronAccept(gap, cfg) == false);
}

int main() {
  RUN(test_selectElectronsP);
  RUN(test_accept_ee_and_gap);
  REPORT();
}
