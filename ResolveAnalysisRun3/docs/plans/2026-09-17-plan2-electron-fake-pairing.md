# Plan 2 — Electron object + Fake(F) 판정 + Candidate pairing

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:executing-plans (inline). Steps는 checkbox 로 추적.

**Goal:** Plan 1(muon P)에 이어, electron P/F·muon F(loose denominator)·flavor별 pairing(4e/4mu/2e2mu)을 순수 C++ + 단위테스트로 구현한다.

**Architecture:** Plan 1과 동일 — CMSSW/데이터 불필요, 합성 struct 로 TDD. 기존 `ResolvedEleCRanalyzer`/`ResolvedMuCRanalyzer`/`PairingHelper` 로직을 정확히 재현.

**Tech Stack:** Plan 1과 동일. 빌드/테스트는 [[scram-dircache-new-srcdir]] 패턴(el8 컨테이너, 새 파일이 아닌 기존 src/에 추가라 DirCache 재삭제 불필요).

**빌드/실행 프리픽스 (모든 Task 공통):**
```bash
/cvmfs/cms.cern.ch/common/cmssw-el8 --command-to-run bash -c \
 'source /cvmfs/cms.cern.ch/cmsset_default.sh; export SCRAM_ARCH=el8_amd64_gcc11;
  cd /u/user/eoyun/4l/25.10.20/CMSSW_13_0_13/src; eval `scram runtime -sh`;
  cd ZprimeTo4l/ResolveAnalysisRun3; scram b -j4; <TEST 바이너리 전체경로>'
```

---

## 기존 코드 재현 근거 (읽어서 확인함)

- **Electron P** (ResolvedEleCRanalyzer.cc:544-572): `|etaSC|<2.5`, EB-EE gap(1.4442~1.566) veto,
  `electronID("modifiedHeepElectronID")==true`.
- **Electron F** (동 574-570): P 아님 + 마스크 통과 `(userInt("modifiedHeepElectronID") | 0x7B0) == 0xFFF`
  (12 HEEP 컷 중 마스크 안 된 것만 요구).
- **Muon F** (ResolvedMuCRanalyzer.cc:692-714): `isTrackerMuon`, track!=null, `tunePpt≥ptThres`,
  `|eta|<2.4`, P collection 제외, `trkLayers>5`, `pixelHits>0`, `|dxy(PV)|<0.2`, `|dz(PV)|<0.5`,
  `matchedStations≥1`. iso 안 봄.
- **Muon acceptance eta** (동 554): `aMuon->eta()` = **reco eta** (Plan 1의 tunePeta 는 수정 필요).
- **Pairing** (PairingHelper.cc): 4e reciprocal-GSF→dR, μ closest-dR, 최종 pair1/pair2 dR 정렬.

---

## Task 1: muon acceptance eta 수정 + isTrackerMuon 필드

**Files:** Modify `interface/Muon.h`, `src/ObjectSelector.cc`, `test/test_muon_selection.cc`

- [ ] **Step 1: Muon struct 에 isTrackerMuon 추가**
`Muon.h` 의 ID 결과 블록에 추가:
```cpp
  bool   isTrackerMuon    = false;  // muon::isTrackerMuon (F loose denominator 용)
```

- [ ] **Step 2: 실패 테스트 추가** (`test_muon_selection.cc` 의 baseMuon 아래, main 위)
```cpp
// recoEta 로 acceptance 판정 (기존 aMuon->eta() 재현). tunePeta 와 다른 값을 줘서 구분.
static void test_accept_uses_recoEta() {
  Config cfg = makeCfg();
  Muon m = baseMuon();
  m.recoEta = 2.5;   // reco eta 는 창 밖(>2.4)
  m.tunePeta = 0.0;  // tunePeta 는 창 안 — 여기에 속으면 안 됨
  CHECK(ObjectSelector::passMuonAccept(m, cfg) == false);
}
```
그리고 `main()` 에 `RUN(test_accept_uses_recoEta);` 추가.

- [ ] **Step 3: passMuonAccept 를 recoEta 로 수정** (`ObjectSelector.cc`)
```cpp
bool passMuonAccept(const Muon& mu, const Config& cfg) {
  return mu.corrTunePpt > cfg.muon.tunePptMin &&
         std::abs(mu.recoEta) < cfg.muon.etaMax;   // 기존 aMuon->eta() 재현
}
```

- [ ] **Step 4: 빌드 + muon 테스트 통과**
Run: 프리픽스 + `$CMSSW_BASE/test/$SCRAM_ARCH/raRun3_test_muon_selection`
Expected: `ALL PASS`.

- [ ] **Step 5: 커밋**
```bash
git add ResolveAnalysisRun3/interface/Muon.h ResolveAnalysisRun3/src/ObjectSelector.cc ResolveAnalysisRun3/test/test_muon_selection.cc
git commit -m "raRun3: fix muon acceptance to use reco eta + add isTrackerMuon"
```

---

## Task 2: Electron struct 확장 + electron P 선택

**Files:** Modify `interface/Electron.h`, `interface/Config.h`, `src/Config.cc`;
Create in `ObjectSelector`; Test `test/test_electron_selection.cc`

- [ ] **Step 1: Electron struct 확장** (`Electron.h` 에 필드 추가)
```cpp
  int    modHeepBitmap = 0;   // userInt("modifiedHeepElectronID") — F 마스크 판정용
  int    addGsfIdx     = -1;  // reciprocal-GSF 상대 electron index (없으면 -1)
```

- [ ] **Step 2: Config 에 electron acceptance + fake 마스크 추가** (`Config.h`)
```cpp
struct ElectronCuts {
  double ebEtaMax = 0.;   // |etaSC| EB 상한 (1.4442)
  double gapLo    = 0.;   // EB-EE gap 하한 (1.4442)
  double gapHi    = 0.;   // EB-EE gap 상한 (1.566)
  double eeEtaMax = 0.;   // |etaSC| EE 상한 (2.5)
  int    heepMaskLoose = 0;  // 0x7B0 = 1968 (F denominator 마스크)
  int    heepAllPass   = 0;  // 0xFFF = 4095 (12 컷 전부 통과)
};
```
`Config` 에 `ElectronCuts electron;` 추가. `Config.cc` 에 파싱 추가(누락시 예외).

- [ ] **Step 3: 실패 테스트** (`test_electron_selection.cc`) — electron P: accept + passModHeep.
```cpp
// baseEle(): |etaSC|=0.5(EB), passModHeep=true
// 테스트: P 통과 / gap 이면 탈락 / |etaSC|>2.5 탈락 / passModHeep=false 탈락
```
(구현 Step 에서 실제 코드 확정)

- [ ] **Step 4~6: `selectElectronsP` 구현, 빌드/통과, 커밋** — Plan 1 패턴과 동일.

---

## Task 3: Electron F (bitmap 마스크)

**Files:** `ObjectSelector` + `test_electron_selection.cc`

- [ ] **Step 1: 실패 테스트** — F: P 아님 + accept + `(bitmap|maskLoose)==allPass`.
- [ ] **Step 2: `selectElectronsF` 구현**
```cpp
// P 아님(passModHeep==false) + acceptance + (modHeepBitmap | maskLoose) == allPass
```
- [ ] **Step 3~4: 빌드/통과, 커밋.**

---

## Task 4: Muon F (loose denominator)

**Files:** `interface/Muon.h`(필드 확인), `ObjectSelector`, `test_muon_selection.cc`

- [ ] **Step 1: 실패 테스트** — F: isTrackerMuon + tunePpt≥min + |recoEta|<max + P 아님
  + trkLayers>5 + pixelHits>0 + |dxy|<0.2 + |dz|<0.5 + matchedStations≥1.
- [ ] **Step 2: `selectMuonsF(all, passP, cfg)` 구현** (passP 를 받아 제외)
```cpp
// P collection(selectMuonsP 결과)에 없는 것 중 위 loose 조건 통과
```
- [ ] **Step 3~4: 빌드/통과, 커밋.**

---

## Task 5: CandidateBuilder — 4mu (closest-dR)

**Files:** Create `interface/CandidateBuilder.h`, `src/CandidateBuilder.cc`, `test/test_candidate.cc`

pairing 규칙(PairingHelper::pairByDR 재현): 3조합 `(01)(23)/(02)(13)/(03)(12)` 중
`min(dR1²,dR2²)` 최소 선택.

- [ ] Candidate struct: `struct Pair { int i, j; }; struct Candidate { Pair pair1, pair2; };`
- [ ] `buildFourMuon(const std::vector<Muon>&) -> Candidate` (정확히 4개 입력).
- [ ] 테스트: 4개 muon 배치로 올바른 조합 선택 검증.

---

## Task 6: CandidateBuilder — 4e (reciprocal-GSF → dR)

pairing 규칙(PairingHelper::pair4E 재현): reciprocal-GSF 짝(`a.addGsfIdx==b && b.addGsfIdx==a`)
먼저 찾고, 있으면 나머지 2개가 pair2. 없으면 dR fallback.

- [ ] `buildFourElectron(const std::vector<Electron>&) -> Candidate`.
- [ ] 테스트: reciprocal-GSF 케이스 / fallback 케이스 둘 다.

---

## Task 7: CandidateBuilder — 2e2mu

pairing: ee 1쌍 + μμ 1쌍 (각 flavor 2개씩이므로 조합 유일).

- [ ] `buildTwoETwoMu(eles, mus) -> Candidate`.
- [ ] 테스트: index 매핑 검증.

---

## Plan 2 완료 조건
- 모든 테스트 `ALL PASS`.
- electron P/F, muon F, 4e/4mu/2e2mu pairing 이 기존 코드 재현.
- 다음 Plan 3에서 NtupleReader + RegionSelector + main 으로 event loop 완성.
