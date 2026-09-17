#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Config.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"
#include <string>
#include <stdexcept>

using raRun3::Config;
using raRun3::loadConfigFromString;

// 정상 JSON은 값이 정확히 파싱돼야 한다.
static void test_parse_ok() {
  const std::string js = R"({
    "objects": {
      "muon": {
        "tunePptMin": 20.0, "etaMax": 2.4, "modIsoRelMax": 0.1,
        "neighbor": { "drMin": 0.01, "drMax": 0.3, "dzMax": 0.2, "dxyBSMax": 0.1 }
      }
    },
    "massCuts": { "dileptonMassMin": 1.0, "signalMassMin": 200.0 }
  })";
  Config cfg = loadConfigFromString(js);
  CHECK_CLOSE(cfg.muon.tunePptMin, 20.0, 1e-9);
  CHECK_CLOSE(cfg.muon.etaMax, 2.4, 1e-9);
  CHECK_CLOSE(cfg.muon.modIsoRelMax, 0.1, 1e-9);
  CHECK_CLOSE(cfg.muon.neighbor.drMax, 0.3, 1e-9);
  CHECK_CLOSE(cfg.muon.neighbor.dxyBSMax, 0.1, 1e-9);
  CHECK_CLOSE(cfg.massCuts.signalMassMin, 200.0, 1e-9);
}

// 누락 키는 조용히 기본값으로 때우지 말고 예외를 던져야 한다(설계 원칙).
static void test_missing_key_throws() {
  const std::string js = R"({ "objects": { "muon": { "etaMax": 2.4 } } })";
  bool threw = false;
  try {
    Config cfg = loadConfigFromString(js);
    (void)cfg;
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);
}

int main() {
  RUN(test_parse_ok);
  RUN(test_missing_key_throws);
  REPORT();
}
