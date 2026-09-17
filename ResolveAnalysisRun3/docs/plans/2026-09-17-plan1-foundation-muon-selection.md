# Plan 1 — Foundation + Muon Object Selection

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run-3 resolved framework의 후단(post-analysis) 기초 — 공유 data struct, JSON config 로더, muon P-object 선택과 modified-isolation 재계산 — 을 CMSSW/데이터 없이 단위 테스트로 검증 가능하게 구현한다.

**Architecture:** 후단은 CMSSW에 의존하지 않는 순수 C++ 라이브러리(`src/`)로 만든다. 물리 로직 함수는 히스토그램·파일·CMSSW 타입을 건드리지 않고 plain struct만 입출력한다. 각 함수는 합성 struct로 단위 테스트한다. ntuple I/O(TTreeReader)와 cmsRun ntuplizer는 이후 Plan에서 이 라이브러리를 재사용한다.

**Tech Stack:** C++17, CMSSW SCRAM 빌드(`scram b`, arch `el9_amd64_gcc12`), ROOT(수학용 `rootmath`), nlohmann/json(scram tool `json`, `#include <nlohmann/json.hpp>`), 자작 최소 테스트 하네스(의존성 0, 전부 읽어 이해 가능).

**전제:** 작업 셸에서 `cd $CMSSW_BASE/src && cmsenv` 를 먼저 실행해야 `scram`/`root`/PATH가 잡힌다. 설계 근거: `ResolveAnalysisRun3/docs/2026-09-17-resolved-run3-framework-design.md`.

---

## File Structure

이 Plan에서 만드는 파일과 책임:

- `ResolveAnalysisRun3/BuildFile.xml` — 패키지 루트. `src/*.cc`를 라이브러리로 빌드하고 의존성(root, json) export.
- `ResolveAnalysisRun3/interface/Kinematics.h` — `P4` struct + 불변질량 함수. 물리 계산의 최소 단위.
- `ResolveAnalysisRun3/interface/Muon.h` — 후단 muon struct (ntuple에서 읽어올 값의 후단 표현).
- `ResolveAnalysisRun3/interface/Event.h` — 후단 event struct (스칼라 + muon/electron 배열).
- `ResolveAnalysisRun3/interface/Electron.h` — electron struct 최소 정의(Plan 2에서 확장). Event가 참조.
- `ResolveAnalysisRun3/interface/Config.h` — config struct (muon/neighbor/mass 섹션).
- `ResolveAnalysisRun3/src/Config.cc` — JSON 파일/문자열 → `Config` 로더. 누락 키는 예외(기본값 금지).
- `ResolveAnalysisRun3/interface/ObjectSelector.h` — muon 선택 함수 선언.
- `ResolveAnalysisRun3/src/ObjectSelector.cc` — muon accept/ID/neighbor/modified-iso/P 선택 구현.
- `ResolveAnalysisRun3/test/TestMain.h` — 최소 테스트 하네스(CHECK 매크로).
- `ResolveAnalysisRun3/test/test_smoke.cc` — 빌드·실행 루프 확인용.
- `ResolveAnalysisRun3/test/test_kinematics.cc` — 불변질량 테스트.
- `ResolveAnalysisRun3/test/test_config.cc` — config 로더 테스트.
- `ResolveAnalysisRun3/test/test_muon_selection.cc` — neighbor/iso/P 선택 테스트.
- `ResolveAnalysisRun3/test/BuildFile.xml` — 테스트 실행파일 빌드 선언.

**테스트 실행 규칙(모든 Task 공통):** `scram b`로 빌드 후, 테스트 바이너리는
`$CMSSW_BASE/test/$SCRAM_ARCH/<name>` 에 생성된다. 전체 경로로 실행한다.

---

## Task 0: 패키지 골격 + 테스트 하네스 (빌드·실행 루프 확인)

**Files:**
- Create: `ResolveAnalysisRun3/BuildFile.xml`
- Create: `ResolveAnalysisRun3/test/TestMain.h`
- Create: `ResolveAnalysisRun3/test/test_smoke.cc`
- Create: `ResolveAnalysisRun3/test/BuildFile.xml`

- [ ] **Step 1: 패키지 루트 BuildFile 작성**

`ResolveAnalysisRun3/BuildFile.xml`:
```xml
<!-- src/*.cc 를 libZprimeTo4l_ResolveAnalysisRun3 로 빌드하고, -->
<!-- 이 라이브러리를 쓰는 쪽(test, bin, plugins)에 의존성을 물려준다. -->
<use name="root"/>       <!-- TMath 등 -->
<use name="rootmath"/>   <!-- ROOT 수학 -->
<use name="json"/>       <!-- nlohmann/json 3.10.2 : #include <nlohmann/json.hpp> -->
<export>
  <lib name="1"/>
</export>
```

- [ ] **Step 2: 최소 테스트 하네스 작성**

`ResolveAnalysisRun3/test/TestMain.h`:
```cpp
#ifndef ResolveAnalysisRun3_TestMain_h
#define ResolveAnalysisRun3_TestMain_h
// 의존성 없는 초경량 테스트 하네스. 프레임워크 대신 전부 읽어 이해할 수 있게 자작.
// CHECK(cond)          : 조건이 거짓이면 파일:라인과 함께 실패 기록
// CHECK_CLOSE(a,b,tol) : |a-b|>tol 이면 실패 기록 (부동소수 비교)
// RUN(fn)              : 테스트 함수 실행 (이름 출력)
// REPORT()             : 실패 개수에 따라 종료코드 반환 (0=성공, 1=실패)
#include <cstdio>
#include <cmath>

static int g_failures = 0;

#define CHECK(cond)                                                            \
  do {                                                                        \
    if (!(cond)) {                                                            \
      std::printf("  [FAIL] %s:%d  CHECK(%s)\n", __FILE__, __LINE__, #cond);  \
      ++g_failures;                                                           \
    }                                                                        \
  } while (0)

#define CHECK_CLOSE(a, b, tol)                                                 \
  do {                                                                        \
    const double _d = std::fabs((double)(a) - (double)(b));                   \
    if (_d > (double)(tol)) {                                                 \
      std::printf("  [FAIL] %s:%d  |%g - %g| = %g > %g\n", __FILE__,          \
                  __LINE__, (double)(a), (double)(b), _d, (double)(tol));     \
      ++g_failures;                                                           \
    }                                                                        \
  } while (0)

#define RUN(fn)                                                                \
  do {                                                                        \
    std::printf("[RUN ] %s\n", #fn);                                          \
    fn();                                                                     \
  } while (0)

#define REPORT()                                                              \
  do {                                                                        \
    if (g_failures) {                                                         \
      std::printf("FAILED: %d check(s)\n", g_failures);                       \
      return 1;                                                               \
    }                                                                        \
    std::printf("ALL PASS\n");                                               \
    return 0;                                                                 \
  } while (0)

#endif
```

- [ ] **Step 3: 스모크 테스트 작성 (일부러 실패시켜 하네스 동작 확인)**

`ResolveAnalysisRun3/test/test_smoke.cc`:
```cpp
#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"

static void test_harness_detects_failure() {
  CHECK(1 + 1 == 2);   // 통과해야 함
  CHECK(1 + 1 == 3);   // 일부러 실패 — Step 4에서 실패를 확인
}

int main() {
  RUN(test_harness_detects_failure);
  REPORT();
}
```

- [ ] **Step 4: 테스트 BuildFile 작성**

`ResolveAnalysisRun3/test/BuildFile.xml`:
```xml
<!-- 각 테스트를 독립 실행파일로 빌드. 패키지 라이브러리(위 export)를 링크한다. -->
<bin file="test_smoke.cc" name="raRun3_test_smoke">
  <use name="ZprimeTo4l/ResolveAnalysisRun3"/>
</bin>
```

- [ ] **Step 5: 빌드하고 스모크 테스트가 "실패"를 잡는지 확인**

Run:
```bash
cd $CMSSW_BASE/src && scram b -j4 2>&1 | tail -5
$CMSSW_BASE/test/$SCRAM_ARCH/raRun3_test_smoke; echo "exit=$?"
```
Expected: 빌드 성공. 실행 출력에 `[FAIL] .../test_smoke.cc:... CHECK(1 + 1 == 3)` 와 `FAILED: 1 check(s)`, `exit=1`.
→ 하네스가 실패를 정확히 잡는다는 증거.

- [ ] **Step 6: 스모크 테스트를 통과 상태로 고치고 재확인**

`test_smoke.cc`의 `CHECK(1 + 1 == 3);` 줄을 삭제한다.

Run:
```bash
cd $CMSSW_BASE/src && scram b -j4 2>&1 | tail -3
$CMSSW_BASE/test/$SCRAM_ARCH/raRun3_test_smoke; echo "exit=$?"
```
Expected: `ALL PASS`, `exit=0`.

- [ ] **Step 7: 커밋**

```bash
cd $CMSSW_BASE/src/ZprimeTo4l
git add ResolveAnalysisRun3/BuildFile.xml ResolveAnalysisRun3/test/
git commit -m "raRun3: package skeleton + minimal test harness"
```

---

## Task 1: Kinematics (P4 + 불변질량)

**Files:**
- Create: `ResolveAnalysisRun3/interface/Kinematics.h`
- Test: `ResolveAnalysisRun3/test/test_kinematics.cc`
- Modify: `ResolveAnalysisRun3/test/BuildFile.xml`

- [ ] **Step 1: 실패하는 테스트 작성**

`ResolveAnalysisRun3/test/test_kinematics.cc`:
```cpp
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Kinematics.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"

using raRun3::P4;
using raRun3::invariantMass;

// 정지-질량이 같은 두 입자를 정반대 방향으로 두면 불변질량은 해석적으로 계산된다.
static void test_back_to_back_dimuon() {
  const double mMu = 0.1056583745;
  // pt=45, eta=0, phi=0 과 phi=pi (정반대). 두 입자 E = sqrt(45^2 + mMu^2).
  P4 a{45.0, 0.0, 0.0, mMu};
  P4 b{45.0, 0.0, M_PI, mMu};
  // 합: px=0, py=0, pz=0, E=2*sqrt(45^2+mMu^2) → M = 2*sqrt(45^2+mMu^2)
  const double eOne = std::sqrt(45.0 * 45.0 + mMu * mMu);
  CHECK_CLOSE(invariantMass(a, b), 2.0 * eOne, 1e-6);
}

// 일반 pair 하나를 손 계산값과 비교 (회귀 방지용 고정값).
// 손 계산: px/py/pz/E 로 전개하면 M ≈ 79.094 GeV (아래 tol 0.05 안).
static void test_known_pair() {
  P4 a{40.0, 0.5, 0.0, 0.1056583745};
  P4 b{35.0, -0.4, 2.5, 0.1056583745};
  CHECK_CLOSE(invariantMass(a, b), 79.094, 0.05);
}

int main() {
  RUN(test_back_to_back_dimuon);
  RUN(test_known_pair);
  REPORT();
}
```

- [ ] **Step 2: 테스트 BuildFile에 항목 추가**

`ResolveAnalysisRun3/test/BuildFile.xml` 에 추가:
```xml
<bin file="test_kinematics.cc" name="raRun3_test_kinematics">
  <use name="ZprimeTo4l/ResolveAnalysisRun3"/>
</bin>
```

- [ ] **Step 3: 최소 구현 작성**

`ResolveAnalysisRun3/interface/Kinematics.h`:
```cpp
#ifndef ResolveAnalysisRun3_Kinematics_h
#define ResolveAnalysisRun3_Kinematics_h
// 무엇: 4-운동량(P4)과 두 입자 불변질량 계산.
// 어떻게: (pt,eta,phi,mass) → (px,py,pz,E) 로 바꿔 더한 뒤 M=sqrt(E^2-|p|^2).
// 의존: 표준 <cmath> 만. CMSSW/ROOT 타입에 의존하지 않아 어디서든 테스트 가능.
#include <cmath>

namespace raRun3 {

// 하나의 lepton kinematics 표현. mass 는 입자 정지질량(예: muon 0.1057).
struct P4 {
  double pt;    // [GeV] 횡운동량 (muon: TuneP, electron: ecalTrk 보정값)
  double eta;   // pseudorapidity
  double phi;   // [rad] 방위각
  double mass;  // [GeV] 정지질량
};

// P4 → 데카르트 성분 (내부 헬퍼)
inline double px(const P4& v) { return v.pt * std::cos(v.phi); }
inline double py(const P4& v) { return v.pt * std::sin(v.phi); }
inline double pz(const P4& v) { return v.pt * std::sinh(v.eta); }
inline double energy(const P4& v) {
  return std::sqrt(px(v) * px(v) + py(v) * py(v) + pz(v) * pz(v) + v.mass * v.mass);
}

// 두 입자의 불변질량 M = sqrt( (E1+E2)^2 - |p1+p2|^2 ).
inline double invariantMass(const P4& a, const P4& b) {
  const double e = energy(a) + energy(b);
  const double sx = px(a) + px(b);
  const double sy = py(a) + py(b);
  const double sz = pz(a) + pz(b);
  const double m2 = e * e - (sx * sx + sy * sy + sz * sz);
  return m2 > 0.0 ? std::sqrt(m2) : 0.0;  // 수치오차로 음수면 0
}

}  // namespace raRun3
#endif
```

- [ ] **Step 4: 빌드 후 실행, 전체 통과 확인**

Run:
```bash
cd $CMSSW_BASE/src && scram b -j4 2>&1 | tail -3
$CMSSW_BASE/test/$SCRAM_ARCH/raRun3_test_kinematics; echo "exit=$?"
```
Expected: `ALL PASS`, `exit=0`. 만약 `test_known_pair`가 `[FAIL]`이면서 출력값이
79.094 와 0.05 이상 다르면 `invariantMass` 공식(px/py/pz/E 전개)을 재검토한다.

- [ ] **Step 5: 커밋**

```bash
cd $CMSSW_BASE/src/ZprimeTo4l
git add ResolveAnalysisRun3/interface/Kinematics.h ResolveAnalysisRun3/test/test_kinematics.cc ResolveAnalysisRun3/test/BuildFile.xml
git commit -m "raRun3: P4 + invariant mass with tests"
```

---

## Task 2: Data struct (Muon/Electron/Event) + Config 로더

**Files:**
- Create: `ResolveAnalysisRun3/interface/Muon.h`
- Create: `ResolveAnalysisRun3/interface/Electron.h`
- Create: `ResolveAnalysisRun3/interface/Event.h`
- Create: `ResolveAnalysisRun3/interface/Config.h`
- Create: `ResolveAnalysisRun3/src/Config.cc`
- Test: `ResolveAnalysisRun3/test/test_config.cc`
- Modify: `ResolveAnalysisRun3/test/BuildFile.xml`

- [ ] **Step 1: 실패하는 config 테스트 작성**

`ResolveAnalysisRun3/test/test_config.cc`:
```cpp
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
```

- [ ] **Step 2: 테스트 BuildFile에 항목 추가 (json 필요)**

`ResolveAnalysisRun3/test/BuildFile.xml` 에 추가:
```xml
<bin file="test_config.cc" name="raRun3_test_config">
  <use name="ZprimeTo4l/ResolveAnalysisRun3"/>
  <use name="json"/>
</bin>
```

- [ ] **Step 3: Muon/Electron/Event struct 작성**

`ResolveAnalysisRun3/interface/Muon.h`:
```cpp
#ifndef ResolveAnalysisRun3_Muon_h
#define ResolveAnalysisRun3_Muon_h
// 무엇: 후단에서 다루는 muon 하나의 표현. ntuple branch 값이 그대로 들어온다.
// 어떻게: 순수 데이터(로직 없음). 단위/의미는 각 필드 주석 참고.
// 의존: 없음. (kinematics 는 P4 로 변환해 Kinematics.h 사용)
namespace raRun3 {

struct Muon {
  int    charge = 0;

  // --- kinematics (Particle Flow 아님: TuneP 사용) ---
  double corrTunePpt = 0.;  // [GeV] 보정 후 TuneP pt (selection/mass 기준)
  double rawTunePpt  = 0.;  // [GeV] 보정 전 TuneP pt (보정 검증용)
  double tunePeta    = 0.;  // TuneP 방향 eta (mass p4)
  double tunePphi    = 0.;  // TuneP 방향 phi (mass p4)

  // --- ID 결과 (CMSSW muon::isHighPtMuon / isTrackerHighPtMuon) ---
  bool   isHighPt        = false;  // global high-pT ID
  bool   isTrackerHighPt = false;  // tracker high-pT ID

  // --- modified isolation 재료 ---
  double trackIso = 0.;   // [GeV] raw track isolation
  double innerPt  = 0.;   // [GeV] inner-track pt (neighbor subtraction에 사용)

  // --- neighbor 기하 (checkIso 재현용) ---
  // self 로 쓸 때: recoEta/recoPhi/recoVz. neighbor 로 쓸 때: innerEta/innerPhi/innerVz/innerDxyBS.
  double recoEta    = 0.;  // muon reco 방향 eta (deltaR 기준: self)
  double recoPhi    = 0.;  // muon reco 방향 phi
  double recoVz     = 0.;  // [cm] muon vertex z (dz 비교: self)
  double innerEta   = 0.;  // inner-track eta (deltaR 기준: neighbor)
  double innerPhi   = 0.;  // inner-track phi
  double innerVz    = 0.;  // [cm] inner-track vertex z (dz 비교: neighbor)
  double innerDxyBS = 0.;  // [cm] inner-track signed dxy wrt beamspot (neighbor)

  // --- ID 입력변수 (nominal Plan1 에선 미사용, 향후 threshold 스캔용) ---
  double relPtErr        = 0.;
  int    trkLayers       = 0;
  int    pixelHits       = 0;
  int    matchedStations = 0;
  double dxy             = 0.;  // [cm] wrt PV
  double dz              = 0.;  // [cm] wrt PV
};

}  // namespace raRun3
#endif
```

`ResolveAnalysisRun3/interface/Electron.h`:
```cpp
#ifndef ResolveAnalysisRun3_Electron_h
#define ResolveAnalysisRun3_Electron_h
// 무엇: 후단 electron 최소 표현. Plan 2에서 modified-HEEP 필드로 확장한다.
// 어떻게: 순수 데이터. kinematics 는 ECAL-driven+GSF (Particle Flow 아님).
// 의존: 없음.
namespace raRun3 {

struct Electron {
  int    charge = 0;
  double selEt  = 0.;   // [GeV] selection용 ET = pat::Electron::et()
  double etaSC  = 0.;   // supercluster eta (acceptance 기준)
  // 보정 전/후 mass kinematics
  double rawPt = 0., rawEta = 0., rawPhi = 0., rawEnergy = 0.;   // polarP4(), energy()
  double corrPt = 0., corrEta = 0., corrPhi = 0., corrM = 0.;    // ecalTrkEnergyPostCorr 적용
  bool   passModHeep = false;  // 기존 VID 결과 (재현 기준)
};

}  // namespace raRun3
#endif
```

`ResolveAnalysisRun3/interface/Event.h`:
```cpp
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
```

- [ ] **Step 4: Config struct + JSON 로더 작성**

`ResolveAnalysisRun3/interface/Config.h`:
```cpp
#ifndef ResolveAnalysisRun3_Config_h
#define ResolveAnalysisRun3_Config_h
// 무엇: 후단의 모든 물리 컷 값을 담는 struct + JSON 로더 선언.
// 어떻게: selection.json → Config. 누락 키는 예외(기본값으로 때우지 않음).
// 의존: 로더 구현(Config.cc)만 nlohmann/json 에 의존. 헤더는 순수 struct.
#include <string>

namespace raRun3 {

struct NeighborCuts {
  double drMin    = 0.;  // deltaR 하한 (self-매칭/중복 track 배제)
  double drMax    = 0.;  // deltaR 상한
  double dzMax    = 0.;  // [cm] |self.recoVz - neighbor.innerVz| 상한
  double dxyBSMax = 0.;  // [cm] neighbor.innerDxyBS 상한 (signed, 기존 checkIso 재현)
};

struct MuonCuts {
  double       tunePptMin = 0.;  // [GeV]
  double       etaMax     = 0.;  // |eta| 상한
  double       modIsoRelMax = 0.;  // modified iso / TuneP pt 상한
  NeighborCuts neighbor;
};

struct MassCuts {
  double dileptonMassMin = 0.;  // [GeV] 각 pair mll 하한
  double signalMassMin   = 0.;  // [GeV] m4l SR 하한
};

struct Config {
  MuonCuts muon;
  MassCuts massCuts;
};

// 파일 경로 / 문자열에서 로드. 필수 키 누락 시 std::runtime_error.
Config loadConfig(const std::string& path);
Config loadConfigFromString(const std::string& jsonText);

}  // namespace raRun3
#endif
```

`ResolveAnalysisRun3/src/Config.cc`:
```cpp
// 무엇: JSON → Config 변환. 누락 키는 명확한 메시지와 함께 예외.
// 어떻게: nlohmann/json 으로 파싱하고 required<T>() 헬퍼로 키 존재를 강제.
// 의존: nlohmann/json (scram tool "json"), Config.h.
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Config.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <nlohmann/json.hpp>

using nlohmann::json;

namespace {

// 키가 없으면 어느 경로에서 무엇이 빠졌는지 알려주고 예외.
template <typename T>
T required(const json& node, const char* key, const std::string& where) {
  if (!node.contains(key))
    throw std::runtime_error("Config: missing key '" + std::string(key) +
                             "' under " + where);
  return node.at(key).get<T>();
}

raRun3::Config parse(const json& j) {
  raRun3::Config cfg;

  if (!j.contains("objects") || !j.at("objects").contains("muon"))
    throw std::runtime_error("Config: missing 'objects.muon'");
  const json& mu = j.at("objects").at("muon");
  cfg.muon.tunePptMin   = required<double>(mu, "tunePptMin", "objects.muon");
  cfg.muon.etaMax       = required<double>(mu, "etaMax", "objects.muon");
  cfg.muon.modIsoRelMax = required<double>(mu, "modIsoRelMax", "objects.muon");

  if (!mu.contains("neighbor"))
    throw std::runtime_error("Config: missing 'objects.muon.neighbor'");
  const json& nb = mu.at("neighbor");
  cfg.muon.neighbor.drMin    = required<double>(nb, "drMin", "objects.muon.neighbor");
  cfg.muon.neighbor.drMax    = required<double>(nb, "drMax", "objects.muon.neighbor");
  cfg.muon.neighbor.dzMax    = required<double>(nb, "dzMax", "objects.muon.neighbor");
  cfg.muon.neighbor.dxyBSMax = required<double>(nb, "dxyBSMax", "objects.muon.neighbor");

  if (!j.contains("massCuts"))
    throw std::runtime_error("Config: missing 'massCuts'");
  const json& mc = j.at("massCuts");
  cfg.massCuts.dileptonMassMin = required<double>(mc, "dileptonMassMin", "massCuts");
  cfg.massCuts.signalMassMin   = required<double>(mc, "signalMassMin", "massCuts");

  return cfg;
}

}  // namespace

namespace raRun3 {

Config loadConfigFromString(const std::string& jsonText) {
  return parse(json::parse(jsonText));
}

Config loadConfig(const std::string& path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Config: cannot open file '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  return loadConfigFromString(ss.str());
}

}  // namespace raRun3
```

- [ ] **Step 5: 빌드 후 config 테스트 통과 확인**

Run:
```bash
cd $CMSSW_BASE/src && scram b -j4 2>&1 | tail -3
$CMSSW_BASE/test/$SCRAM_ARCH/raRun3_test_config; echo "exit=$?"
```
Expected: `ALL PASS`, `exit=0`. (`test_missing_key_throws`가 예외를 확인)

- [ ] **Step 6: 커밋**

```bash
cd $CMSSW_BASE/src/ZprimeTo4l
git add ResolveAnalysisRun3/interface/Muon.h ResolveAnalysisRun3/interface/Electron.h \
        ResolveAnalysisRun3/interface/Event.h ResolveAnalysisRun3/interface/Config.h \
        ResolveAnalysisRun3/src/Config.cc ResolveAnalysisRun3/test/test_config.cc \
        ResolveAnalysisRun3/test/BuildFile.xml
git commit -m "raRun3: data structs + JSON config loader with tests"
```

---

## Task 3: Muon neighbor 판정 + modified isolation + P 선택

**Files:**
- Create: `ResolveAnalysisRun3/interface/ObjectSelector.h`
- Create: `ResolveAnalysisRun3/src/ObjectSelector.cc`
- Test: `ResolveAnalysisRun3/test/test_muon_selection.cc`
- Modify: `ResolveAnalysisRun3/test/BuildFile.xml`

기존 로직(재현 대상):
- `checkIso` (Analysis/src/MuonCorrectionHelper.cc): `0.0001 < dR2 < 0.09`,
  `|mu.vz - trk.vz| < 0.2`, `trk.dxy(bs) <= 0.1` (signed).
- modified iso (Analysis/plugins/ResolvedMuCRanalyzer.cc:571-605): ID-pass muon 중
  self 제외 checkIso 통과 후보를 TuneP pt 내림차순 정렬, 가장 높은 하나의 innerPt 를
  raw trackIso 에서 뺀 뒤 `iso/tunePpt < 0.1` 이면 isolated P.

- [ ] **Step 1: 실패하는 테스트 작성**

`ResolveAnalysisRun3/test/test_muon_selection.cc`:
```cpp
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
```

- [ ] **Step 2: 테스트 BuildFile에 항목 추가**

`ResolveAnalysisRun3/test/BuildFile.xml` 에 추가:
```xml
<bin file="test_muon_selection.cc" name="raRun3_test_muon_selection">
  <use name="ZprimeTo4l/ResolveAnalysisRun3"/>
</bin>
```

- [ ] **Step 3: ObjectSelector 헤더 작성**

`ResolveAnalysisRun3/interface/ObjectSelector.h`:
```cpp
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
```

- [ ] **Step 4: ObjectSelector 구현 작성**

`ResolveAnalysisRun3/src/ObjectSelector.cc`:
```cpp
// 무엇: muon 선택 로직 구현. 기존 checkIso / ResolvedMuCRanalyzer 재현.
// 어떻게: neighbor 후보를 TuneP pt 로 정렬해 최고 하나만 빼는 modified iso.
// 의존: ObjectSelector.h, <algorithm>, <cmath>.
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/ObjectSelector.h"

#include <algorithm>
#include <cmath>

namespace {

// deltaR^2 (phi 는 [-pi,pi] 로 감싼다).
double deltaR2(double eta1, double phi1, double eta2, double phi2) {
  double dphi = phi1 - phi2;
  while (dphi >  M_PI) dphi -= 2.0 * M_PI;
  while (dphi < -M_PI) dphi += 2.0 * M_PI;
  const double deta = eta1 - eta2;
  return deta * deta + dphi * dphi;
}

}  // namespace

namespace raRun3 {
namespace ObjectSelector {

bool passMuonAccept(const Muon& mu, const Config& cfg) {
  return mu.corrTunePpt > cfg.muon.tunePptMin &&
         std::abs(mu.tunePeta) < cfg.muon.etaMax;
}

bool passMuonId(const Muon& mu) {
  return mu.isHighPt || mu.isTrackerHighPt;
}

bool isMuonNeighbor(const Muon& self, const Muon& other, const Config& cfg) {
  // self reco 방향 vs neighbor inner-track 방향 (기존 checkIso 와 동일)
  const double dr2 = deltaR2(self.recoEta, self.recoPhi, other.innerEta, other.innerPhi);
  const double drMin2 = cfg.muon.neighbor.drMin * cfg.muon.neighbor.drMin;
  const double drMax2 = cfg.muon.neighbor.drMax * cfg.muon.neighbor.drMax;
  if (dr2 > drMax2 || dr2 < drMin2) return false;
  if (std::abs(self.recoVz - other.innerVz) > cfg.muon.neighbor.dzMax) return false;
  // 기존 checkIso 는 signed dxy 를 그대로 비교(> 면 탈락). 재현.
  if (other.innerDxyBS > cfg.muon.neighbor.dxyBSMax) return false;
  return true;
}

double modifiedMuonIso(const Muon& self, const std::vector<Muon>& idPool, const Config& cfg) {
  // self 와 다른(주소 아닌 값 비교 대신 동일 객체 판정을 위해 innerPt+방향으로 구분하지 않고,
  // 호출자가 idPool 에 self 를 포함해 넘겨도 되도록 checkIso 의 drMin(자기 자신 dR≈0 배제)이 self 를 자동 제외).
  double iso = self.trackIso;

  // neighbor 후보 수집
  std::vector<const Muon*> cands;
  for (const auto& other : idPool) {
    if (isMuonNeighbor(self, other, cfg))
      cands.push_back(&other);
  }
  if (cands.empty()) return iso;

  // TuneP pt 내림차순 정렬 후 최고 하나의 innerPt 를 뺀다.
  std::sort(cands.begin(), cands.end(),
            [](const Muon* a, const Muon* b) { return a->corrTunePpt > b->corrTunePpt; });
  iso -= cands.front()->innerPt;
  return iso;
}

std::vector<Muon> selectMuonsP(const std::vector<Muon>& all, const Config& cfg) {
  // 1) accept + ID 통과 pool
  std::vector<Muon> idPool;
  for (const auto& mu : all) {
    if (passMuonAccept(mu, cfg) && passMuonId(mu))
      idPool.push_back(mu);
  }
  // 2) 각 pool muon 의 modified iso 로 격리 판정
  std::vector<Muon> passP;
  for (const auto& mu : idPool) {
    const double iso = modifiedMuonIso(mu, idPool, cfg);
    if (iso / mu.corrTunePpt < cfg.muon.modIsoRelMax)
      passP.push_back(mu);
  }
  return passP;
}

}  // namespace ObjectSelector
}  // namespace raRun3
```

> **주의(자기 자신 제외):** `modifiedMuonIso`는 idPool 에 self 가 포함돼도 `isMuonNeighbor`의
> `drMin`(dR≈0 배제) 덕에 self 를 neighbor 로 세지 않는다. 단, self 와 방향·vz 가 우연히
> 완전히 같은 **다른** muon이 있으면 배제될 수 있는데, 실제 event 에서 서로 다른 muon 이
> dR<0.01 로 겹치는 경우는 기존 코드도 동일하게 배제하므로 재현 관점에서 일치한다.

- [ ] **Step 5: 빌드 후 muon 선택 테스트 통과 확인**

Run:
```bash
cd $CMSSW_BASE/src && scram b -j4 2>&1 | tail -3
$CMSSW_BASE/test/$SCRAM_ARCH/raRun3_test_muon_selection; echo "exit=$?"
```
Expected: `ALL PASS`, `exit=0`.

- [ ] **Step 6: 전체 테스트 재실행 (회귀 확인)**

Run:
```bash
for t in raRun3_test_smoke raRun3_test_kinematics raRun3_test_config raRun3_test_muon_selection; do
  echo "== $t =="; $CMSSW_BASE/test/$SCRAM_ARCH/$t; echo "exit=$?";
done
```
Expected: 네 개 모두 `ALL PASS`, `exit=0`.

- [ ] **Step 7: 커밋**

```bash
cd $CMSSW_BASE/src/ZprimeTo4l
git add ResolveAnalysisRun3/interface/ObjectSelector.h ResolveAnalysisRun3/src/ObjectSelector.cc \
        ResolveAnalysisRun3/test/test_muon_selection.cc ResolveAnalysisRun3/test/BuildFile.xml
git commit -m "raRun3: muon neighbor/modified-iso/P selection with tests"
```

---

## Plan 1 완료 조건

- `scram b` 성공, 네 테스트 실행파일 모두 `ALL PASS`.
- muon P 선택과 modified isolation 이 기존 `checkIso`/`ResolvedMuCRanalyzer` 로직을
  재현(neighbor 창, dz/dxy, 최고-TuneP innerPt subtraction, iso/pt 컷).
- 이후 모든 후단 모듈이 따를 패턴 확립: 순수 함수 + 합성 struct 단위 테스트 + config 주도 컷.

---

## 후속 Plan 로드맵 (각각 별도 문서로 작성)

- **Plan 2 — Electron object + Candidate pairing**
  modified-HEEP 재현(입력변수 저장분 활용) + electron/muon F(loose denominator) 판정
  + `CandidateBuilder`(4e reciprocal-GSF→dR, μ closest-dR, PairingHelper 재현).
- **Plan 3 — NtupleReader + runResolved main + nominal 4P SR**
  `Events.root`(합성/실제) 읽기, RegionSelector(trigger match·multiplicity·pair-ID·
  mll/m4l·SR/CR), HistogramWriter, `bin/runResolved.cc`. selection.json round-trip 검증.
- **Plan 4 — ResolvedNtuplizer (cmsRun EDAnalyzer)**
  MiniAOD → Events.root. correction 전/후 저장, 느슨한 skim(3-lepton 보존),
  metadata/정규화. 작은 ZZ sample 로 통합 검증(branch 존재·전후 보정 합리성).
- **Plan 5 — CR/FF 적용 + Systematics + BackgroundBuilder**
  3P1F/2P2F, isolated/nonisolated FF, weight-only + kinematic(migration) systematics,
  prompt subtraction/ template(clipping 정책). ZZ Z-peak 검증.
