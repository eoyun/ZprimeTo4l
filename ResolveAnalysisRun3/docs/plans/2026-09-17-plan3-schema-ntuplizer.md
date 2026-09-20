# Plan 3 — NtupleSchema (공유 계약) + ResolvedNtuplizer (cmsRun plugin)

> REQUIRED SUB-SKILL: superpowers:executing-plans (inline). Steps는 checkbox 로 추적.

**Goal:** MiniAOD → `Events.root` ntuple 을 생산하는 CMSSW plugin 을 만든다.
후단(Plan 1+2)의 `Muon`/`Electron`/`Event` struct 에 대응하는 branch 를,
공유 헤더 `NtupleSchema.h` 로 이름을 고정해 채운다. correction 전/후 + variation 저장,
느슨한 skim(3-lepton 보존).

**Architecture:** ntuplizer 는 `plugins/ResolvedNtuplizer.cc` **단일 CMSSW plugin**
(EDAnalyzer). src+interface 안 씀. branch 이름은 `interface/NtupleSchema.h` 한 곳에서만
정의해 writer(ntuplizer)와 reader(Plan 4)가 공유 → drift 방지.

**Tech Stack:** CMSSW EDAnalyzer, TFileService, TTree(std::vector 가변길이 branch),
기존 `MuonCorrectionHelper`/`ElectronSystematicsHelper`(correction), el8 컨테이너 빌드.

**전제:** 실행 검증엔 Run-3 MiniAOD 1개 + GT 필요(사용자 제공 예정). 그 전까지는
**빌드(컴파일)까지만 검증**한다.

---

## 기존 코드에서 각 값 추출 방법 (재현 근거)

- muon TuneP: `aMu->tunePMuonBestTrack()->pt()/eta()/phi()`; inner: `aMu->innerTrack()->pt()/eta()/phi()/vz()`; reco: `aMu->eta()/phi()/vz()`; ID: `muon::isHighPtMuon(*aMu,pv)`, `muon::isTrackerHighPtMuon`, `aMu->isTrackerMuon()`; iso: `aMu->trackIso()`; ID inputs: `innerTrack()->hitPattern().trackerLayersWithMeasurement()/numberOfValidPixelHits()`, `numberOfMatchedStations()`, `innerTrack()->dxy(pv)/dz(pv)`; dxyBS: `innerTrack()->dxy(bs.position())`; correction: `MuonCorrectionHelper::nominalMC/rochesterMC` (전/후 비교용 raw=tuneP pt, corr=보정 적용).
- electron: `aEle->et()`(selEt), `superCluster()->eta()`(etaSC), `polarP4()`(raw), `energy()`, `userFloat("ecalTrkEnergyPreCorr"/"ecalTrkEnergyPostCorr")`(corr), `electronID("modifiedHeepElectronID")`(passModHeep bool), `userInt("modifiedHeepElectronID")`(bitmap); addGsfIdx: `ValueMap<GsfTrackRef> addGsfTrkMap` 로 이웃 electron 매칭.
- trigger: `TriggerResults` decision, `TriggerObjectStandAlone` unpackPathNames.
- weight: `GenEventInfoProduct::weight()`(genWeight), `PileupSummaryInfo`(puTrue).

---

## Task 1: NtupleSchema.h (branch 이름 공유 계약)

**Files:** Create `interface/NtupleSchema.h`; Test `test/test_schema.cc`; Modify `test/BuildFile.xml`

- [ ] **Step 1: 실패 테스트** — schema 이름이 비어있지 않고 유일한지 최소 확인.
```cpp
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/NtupleSchema.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"
#include <string>
static void test_names_nonempty() {
  CHECK(std::string(raRun3::schema::mu::corrTunePpt).rfind("mu_",0)==0);
  CHECK(std::string(raRun3::schema::ele::selEt).rfind("ele_",0)==0);
  CHECK(raRun3::schema::kVersion >= 1);
}
int main(){ RUN(test_names_nonempty); REPORT(); }
```

- [ ] **Step 2: NtupleSchema.h 작성** — 모든 branch 이름 상수 + 버전. (구현 Step 에서 전체 확정;
  Muon/Electron struct 의 각 필드에 1:1 대응하는 이름, event/trigger 포함.)

- [ ] **Step 3~4: 빌드/통과, 커밋.**

---

## Task 2: ResolvedNtuplizer 골격 + event 스칼라 + muon branch

**Files:** Create `plugins/ResolvedNtuplizer.cc`, `plugins/BuildFile.xml`, `python/ResolvedNtuplizer_cfi.py`

- [ ] EDAnalyzer 골격: consumes(muon/electron/pv/beamspot/gen/pileup/triggerResults/triggerObjects),
  `beginJob` 에서 TTree + `NtupleSchema` branch 등록, `analyze` 에서 채우고 `Fill()`.
- [ ] event 스칼라: run/lumi/event/genWeight/puTrue.
- [ ] muon vector branch 채우기: charge, raw/corr TuneP, reco/inner 방향·vz, innerDxyBS,
  innerPt, trackIso, ID bool 3종, ID inputs, index. (correction 전/후: raw=tuneP, corr=helper).
- [ ] **검증(파일 없을 때)**: `scram b` 로 plugin 컴파일 성공.

---

## Task 3: electron branch

- [ ] electron vector branch: charge, selEt, etaSC, raw p4, corr p4(ecalTrkEnergyPostCorr),
  passModHeep, modHeepBitmap, addGsfIdx(addGsfTrkMap 로 매칭), index.
- [ ] 검증: 컴파일 성공.

---

## Task 4: trigger object + HLT decision + METfilter + metadata/정규화

- [ ] trigger object branch(pt/eta/phi/filter or path id), HLT decision(설정된 path bool),
  METfilter pass, metadata TTree(schemaVersion/sourceCommit/appliedSkim),
  정규화 TH1(Nevents/sumGenWeights/sumGenWeights2, skim 이전).
- [ ] 느슨한 skim: >=3 lepton(P∪F loose) 유지 등(정확 값은 사용자 조사).
- [ ] 검증: 컴파일 성공.

---

## Task 5: cfg.py + 실제 파일 검증

- [ ] `test/runResolvedNtuplizer_cfg.py`: source(사용자 dataset), GT, producer sequence, TFileService.
- [ ] 작은 Run-3 MiniAOD 1개로 실행 → branch 존재, `corr/raw` 합리성, 3-lepton event 보존 확인.

---

## Plan 3 완료 조건
- plugin 컴파일 성공. (파일 확보 후) 작은 sample 로 ntuple 생산 + 검증.
- `NtupleSchema.h` 가 writer/reader 공유 계약으로 확립 → Plan 4 reader 가 동일 이름으로 읽음.
