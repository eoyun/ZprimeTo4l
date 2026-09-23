# Resolved μ CR analyzer — Run-3 포팅 (템플릿) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `Analysis/plugins/ResolvedMuCRanalyzer.cc`(Run-2 UL)를 Run-3(2022 preEE)용 `ResolvedMuCRanalyzerRun3`로 포팅한다 — prefiring 제거, egamma VID-only, Run-3 GT/era/trigger/MET-filter, MC/Data 통합 cfg.

**Architecture:** MiniAOD monolithic EDAnalyzer를 파일 3종(.cc/cfi/cfg)으로 복제·수정. 물리 로직(selection·pairing·FF·SF systematic)은 원본과 100% 동일하게 유지하고, prefiring과 Run-2 전용 cfg/cfi 요소만 Run-3로 교체한다. helper(PairingHelper/MuonCorrectionHelper/RoccoR)는 무변경 재사용.

**Tech Stack:** CMSSW 13_0_13 (C++17 EDAnalyzer, cmsRun python cfg), correctionlib, ROOT, egammaPostRecoTools(VID), scram build (cmssw-el8 컨테이너).

**검증 방식 주의:** EDAnalyzer는 단위 테스트 하네스가 없다. 검증 = ① `scram b` 빌드 통과 ② `cmsRun` 실행 crash 없음 ③ prefiring 흔적 grep 0 ④ Z peak sanity. TDD의 "failing test"는 여기서 "빌드/실행 실패를 먼저 관찰 → 수정 → 통과 관찰"로 치환한다.

**참고 권위 레퍼런스:** `ResolveAnalysisRun3/test/26.09.21_22preEE/runResolvedNtuplizer_run3MC_cfg.py` (Run-3 2022preEE GT/era/trigger/MET-filter 값의 원천).

**Spec:** `docs/superpowers/specs/2026-09-23-resolved-mu-cr-run3-port-design.md`

---

## File Structure

| 파일 | 책임 | 조치 |
|---|---|---|
| `Analysis/plugins/ResolvedMuCRanalyzerRun3.cc` | Run-3 μ CR EDAnalyzer 본체 | Create (원본 복사 후 수정) |
| `Analysis/python/ResolvedMuCRanalyzerRun3_cfi.py` | 모듈 기본 PSet (Run-3 payload/trigger/filter) | Create |
| `Analysis/test/runResolvedMuCRRun3_cfg.py` | cmsRun 실행 설정 (MC/Data 통합) | Create |
| `Analysis/plugins/BuildFile.xml` | 빌드 — `*.cc` glob이라 **수정 불필요** (확인만) | 무변경 |
| 기존 Run-2 3파일 | 회귀 방지 | 무변경 |

---

## Task 1: `.cc` 복사 + 클래스명 변경 (prefiring 제외, 그대로 복제)

**Files:**
- Create: `Analysis/plugins/ResolvedMuCRanalyzerRun3.cc` (from `Analysis/plugins/ResolvedMuCRanalyzer.cc`)

- [ ] **Step 1: 원본을 새 파일로 복사**

```bash
cd /u/user/eoyun/4l/25.10.20/CMSSW_13_0_13/src/ZprimeTo4l
cp Analysis/plugins/ResolvedMuCRanalyzer.cc Analysis/plugins/ResolvedMuCRanalyzerRun3.cc
```

- [ ] **Step 2: 클래스명·생성자·소멸자·모듈등록을 `ResolvedMuCRanalyzerRun3`로 일괄 변경**

`ResolvedMuCRanalyzerRun3.cc`에서 다음을 치환(단어 경계 정확히):
- 클래스 선언: `class ResolvedMuCRanalyzer :` → `class ResolvedMuCRanalyzerRun3 :`
- 생성자 선언: `explicit ResolvedMuCRanalyzer(` → `explicit ResolvedMuCRanalyzerRun3(`
- 소멸자: `virtual ~ResolvedMuCRanalyzer()` → `virtual ~ResolvedMuCRanalyzerRun3()`
- 생성자 정의: `ResolvedMuCRanalyzer::ResolvedMuCRanalyzer(` → `ResolvedMuCRanalyzerRun3::ResolvedMuCRanalyzerRun3(`
- 모든 멤버함수 정의 스코프: `ResolvedMuCRanalyzer::` → `ResolvedMuCRanalyzerRun3::` (beginJob/analyze/endJob/lvecFromTuneP)
- 파일 맨 끝 등록: `DEFINE_FWK_MODULE(ResolvedMuCRanalyzer);` → `DEFINE_FWK_MODULE(ResolvedMuCRanalyzerRun3);`

명령 예(정확 치환, `::`와 선언/등록 모두 포함하도록 `\b` 경계 사용):
```bash
sed -i -E 's/\bResolvedMuCRanalyzer\b/ResolvedMuCRanalyzerRun3/g' Analysis/plugins/ResolvedMuCRanalyzerRun3.cc
```
> 이 sed는 `ResolvedMuCRanalyzer` 전체 토큰만 바꾸므로 안전(다른 심볼에 부분매치 없음).

- [ ] **Step 3: 치환 검증**

Run:
```bash
grep -c 'ResolvedMuCRanalyzerRun3' Analysis/plugins/ResolvedMuCRanalyzerRun3.cc
grep -c '\bResolvedMuCRanalyzer\b' Analysis/plugins/ResolvedMuCRanalyzerRun3.cc
tail -1 Analysis/plugins/ResolvedMuCRanalyzerRun3.cc
```
Expected: 첫 grep ≥ 8, 둘째 grep = 0, 마지막 줄 = `DEFINE_FWK_MODULE(ResolvedMuCRanalyzerRun3);`

- [ ] **Step 4: Commit**

```bash
git add Analysis/plugins/ResolvedMuCRanalyzerRun3.cc
git commit -m "raRun3: copy ResolvedMuCRanalyzer -> Run3 variant (rename only)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `.cc` prefiring 완전 제거

**Files:**
- Modify: `Analysis/plugins/ResolvedMuCRanalyzerRun3.cc`

> 주의: 원본 라인번호는 참고용. Task 1 후 라인이 밀리지 않았으므로 동일하지만, **문자열로 매칭**하여 편집한다.

- [ ] **Step 1: 멤버 토큰 3개 선언 삭제**

`ResolvedMuCRanalyzerRun3.cc`에서 아래 3줄(연속) 삭제:
```cpp
  const edm::EDGetTokenT<float> prefweight_token;
  const edm::EDGetTokenT<float> prefweightUp_token_;
  const edm::EDGetTokenT<float> prefweightDn_token_;
```

- [ ] **Step 2: 생성자 초기화 리스트 3줄 삭제**

아래 3줄(연속) 삭제:
```cpp
prefweight_token(consumes<float>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
prefweightUp_token_(consumes<float>(edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
prefweightDn_token_(consumes<float>(edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
```

- [ ] **Step 3: histogram 부킹 2줄 삭제**

아래 2줄 삭제:
```cpp
  histo1d_["4P0F_CR_llll_invM_prefireUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_prefireUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_prefireDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_prefireDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
```

- [ ] **Step 4: analyze() 지역 prefire 변수 3줄 삭제**

아래 3줄 삭제:
```cpp
  float prefireNo = 1.;
  float prefireUp = 1.;
  float prefireDn = 1.;
```

- [ ] **Step 5: analyze() prefire getByToken 블록 삭제**

아래 연속 블록 삭제:
```cpp
    edm::Handle<float> theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight);
    prefireNo = *theprefweight;

    edm::Handle<float> theprefweightUp;
    iEvent.getByToken(prefweightUp_token_, theprefweightUp);
    prefireUp = *theprefweightUp;

    edm::Handle<float> theprefweightDn;
    iEvent.getByToken(prefweightDn_token_, theprefweightDn);
    prefireDn = *theprefweightDn;

```

- [ ] **Step 6: weight 조합식에서 prefireNo 제거**

찾기:
```cpp
    aWeight = prefireNo*mcweight/std::abs(mcweight);
```
바꾸기:
```cpp
    aWeight = mcweight/std::abs(mcweight);
```

- [ ] **Step 7: 4P0F prefire systematic fill 2줄 삭제**

아래 2줄 삭제(주변 PU fill 라인은 유지):
```cpp
      histo1d_["4P0F_CR_llll_invM_prefireUp"]->Fill( m4l, aWeight*prefireUp/prefireNo );
      histo1d_["4P0F_CR_llll_invM_prefireDn"]->Fill( m4l, aWeight*prefireDn/prefireNo );
```

- [ ] **Step 8: prefiring 흔적 0 확인**

Run:
```bash
grep -in 'prefir' Analysis/plugins/ResolvedMuCRanalyzerRun3.cc
```
Expected: 출력 없음 (exit 1).

- [ ] **Step 9: Commit**

```bash
git add Analysis/plugins/ResolvedMuCRanalyzerRun3.cc
git commit -m "raRun3: remove prefiring from Run-3 mu CR analyzer

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: cfi 작성 (Run-3 2022preEE, payload placeholder+TODO)

**Files:**
- Create: `Analysis/python/ResolvedMuCRanalyzerRun3_cfi.py`

> 원본 `ResolvedMuCRanalyzer_cfi.py`의 파라미터 이름/타입을 그대로 따르되(analyzer가 그 이름으로 읽음), 값만 Run-3화. payload 경로는 repo에 Run-3 버전이 없으므로 Run-2 UL을 placeholder로 두고 `# TODO(Run-3 2022preEE)` 주석. 전체 파라미터 목록은 원본 cfi에서 확인하고 누락 없이 옮긴다.

- [ ] **Step 1: 원본 파라미터 목록 확인 (누락 방지)**

Run:
```bash
grep -nE '= *cms\.' Analysis/python/ResolvedMuCRanalyzer_cfi.py | head -60
```
Expected: `isMC, srcEle, srcMuon, srcMET, srcPv, genptc, generator, pileupSummary, beamSpot, METfilters, triggerResults, triggerObjects, METfilterList, trigList, rochesterPath, triggerSF, muonIdIsoSFpath, muonBoostIsoSFpath, muonRecoSFpath, FFpath, PUrwgt, PUname, muScaleBias, muSmearFactors, muSmearParams, ptThres, ptThresTrig, ffSystCL95` 등 전 파라미터 확인. **이 목록 전체가 새 cfi에 존재해야 한다.**

- [ ] **Step 2: 새 cfi 작성**

`Analysis/python/ResolvedMuCRanalyzerRun3_cfi.py`:
```python
import FWCore.ParameterSet.Config as cms

# ResolvedMuCRanalyzerRun3 — Run-3 (2022 preEE) 4mu resolved CR analyzer.
# 물리 로직은 Run-2 ResolvedMuCRanalyzer 와 동일. prefiring 제거, egamma VID-only.
# payload(SF/PU/RoccoR/FF)는 아직 Run-3 버전 부재 → Run-2 UL placeholder + TODO.
resolvedMuCRanalyzerRun3 = cms.EDAnalyzer('ResolvedMuCRanalyzerRun3',
    isMC            = cms.bool(True),
    srcEle          = cms.InputTag('slimmedElectrons'),
    srcMuon         = cms.InputTag('slimmedMuons'),
    srcMET          = cms.InputTag('slimmedMETs'),
    srcPv           = cms.InputTag('offlineSlimmedPrimaryVertices'),
    genptc          = cms.InputTag('prunedGenParticles'),
    generator       = cms.InputTag('generator'),
    pileupSummary   = cms.InputTag('slimmedAddPileupInfo'),
    beamSpot        = cms.InputTag('offlineBeamSpot'),
    METfilters      = cms.InputTag('TriggerResults', '', 'PAT'),
    triggerResults  = cms.InputTag('TriggerResults', '', 'HLT'),
    triggerObjects  = cms.InputTag('slimmedPatTrigger'),

    # Run-3 MET filter (HBHENoise/HBHENoiseIso 제거, Run-3 필터 추가)
    METfilterList = cms.vstring(
        'Flag_goodVertices',
        'Flag_globalSuperTightHalo2016Filter',
        'Flag_EcalDeadCellTriggerPrimitiveFilter',
        'Flag_BadPFMuonFilter',
        'Flag_BadPFMuonDzFilter',
        'Flag_hfNoisyHitsFilter',
        'Flag_eeBadScFilter',
        'Flag_ecalBadCalibFilter',
    ),
    # Run-3 single-mu / high-pT-mu triggers (Run-2 TkMu50 부재 -> Cascade/HighPtTk)
    trigList = cms.vstring(
        'HLT_Mu50_v*',
        'HLT_CascadeMu100_v*',
        'HLT_HighPtTkMu100_v*',
    ),

    # ---- payloads: Run-3 버전 부재 → Run-2 UL placeholder (TODO 교체) ----
    # TODO(Run-3 2022preEE): RoccoR Run-3 공식 파일 부재. 유지, muon momentum 마이그레이션은 후속.
    rochesterPath     = cms.FileInPath('ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt'),
    # TODO(Run-3 2022preEE): Muon POG Run-3 highPt SF json 으로 교체.
    triggerSF         = cms.FileInPath('ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2016_schemaV2.json'),
    muonIdIsoSFpath   = cms.FileInPath('ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2016_schemaV2.json'),
    muonBoostIsoSFpath= cms.FileInPath('ZprimeTo4l/Analysis/data/muIso20UL16.root'),
    muonRecoSFpath    = cms.FileInPath('ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2016_schemaV2.json'),
    # TODO(Run-3 2022preEE): FF 재측정 payload 로 교체.
    FFpath            = cms.FileInPath('ZprimeTo4l/Analysis/data/RMFF_run2.root'),
    # TODO(Run-3 2022preEE): LUM POG Run-3 PU json 으로 교체.
    PUrwgt            = cms.FileInPath('ZprimeTo4l/Analysis/data/puWeights20UL16.json'),
    PUname            = cms.string('Collisions16_UltraLegacy_goldenJSON'),

    # 고pT muon scale/smear (GE-method). TODO(Run-3): Run-3 재검증. 값은 Run-2 유지.
    muScaleBias   = cms.vdouble(),   # Step 3에서 원본 값으로 채움
    muSmearFactors= cms.vdouble(),   # Step 3에서 원본 값으로 채움
    muSmearParams = cms.vdouble(),   # Step 3에서 원본 값으로 채움

    ptThres     = cms.double(20.),
    ptThresTrig = cms.double(52.),
    ffSystCL95  = cms.double(1.96),
)
```

- [ ] **Step 3: `muScaleBias`/`muSmearFactors`/`muSmearParams`/`PUname` 실제 값을 원본에서 복사**

Run (원본 2016b 블록 값 확인):
```bash
grep -nE 'muScaleBias|muSmearFactors|muSmearParams|PUname|PUrwgt|ffSystCL95|ptThres' Analysis/python/ResolvedMuCRanalyzer_cfi.py | head -40
sed -n '/RoccoR2016bUL/,/RMFF/p' Analysis/python/ResolvedMuCRanalyzer_cfi.py
```
확인한 `muScaleBias`(18개), `muSmearFactors`(2개), `muSmearParams`(8개) vdouble 값과 `PUname`/`ffSystCL95`/`ptThres*` 값을 위 cfi의 빈 `cms.vdouble()`/placeholder에 그대로 채운다. (원본 2016b 변형 블록 기준. 값 자체는 detector-level이라 Run-2 값 유지가 의도.)

- [ ] **Step 4: python 구문 검증**

Run:
```bash
cd /u/user/eoyun/4l/25.10.20/CMSSW_13_0_13/src/ZprimeTo4l
python3 -c "import ast; ast.parse(open('Analysis/python/ResolvedMuCRanalyzerRun3_cfi.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 5: 파라미터 누락 대조**

Run:
```bash
diff <(grep -oE '^ *[a-zA-Z0-9_]+ *= *cms\.' Analysis/python/ResolvedMuCRanalyzer_cfi.py | grep -oE '[a-zA-Z0-9_]+' | grep -v cms | sort -u) \
     <(grep -oE '^ *[a-zA-Z0-9_]+ *= *cms\.' Analysis/python/ResolvedMuCRanalyzerRun3_cfi.py | grep -oE '[a-zA-Z0-9_]+' | grep -v cms | sort -u)
```
Expected: 차이 없음(빈 출력). 차이가 있으면 누락 파라미터를 새 cfi에 추가.

- [ ] **Step 6: Commit**

```bash
git add Analysis/python/ResolvedMuCRanalyzerRun3_cfi.py
git commit -m "raRun3: cfi for Run-3 mu CR analyzer (2022preEE, payload TODO)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: cfg 작성 (MC/Data 통합, Run-3 2022preEE)

**Files:**
- Create: `Analysis/test/runResolvedMuCRRun3_cfg.py`

- [ ] **Step 1: cfg 작성**

`Analysis/test/runResolvedMuCRRun3_cfg.py`:
```python
# =============================================================================
# ResolvedMuCRanalyzerRun3 cfg (Run 3, 2022 preEE) — MC/Data 통합
#
# isMC 스위치 하나로 GT/era/analyzer.isMC 분기. egamma 는 VID-only
# (runEnergyCorrections=False) — mu 채널은 electron 을 veto ID 로만 사용.
# prefiring 없음.
#
# 실행:  cmsRun runResolvedMuCRRun3_cfg.py
# 출력:  hists.root
# =============================================================================
import FWCore.ParameterSet.Config as cms

isMC = True   # Data 로 돌릴 땐 False 로 (아래 GT/METfilters 태그 분기)

# Run-3 era modifier 필수 (Run-3 GT 엔 PCastorRcd 없음 → run3_common 이 CASTOR 제외)
from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('ResolvedMuCRanalysisRun3', Run3)

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        # 로컬 테스트용 예시 (Run3Summer22 DY MiniAODSIM). crab/실행 시 override.
        '/store/mc/Run3Summer22MiniAODv4/DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/017afd05-e111-41dc-9802-86e708952417.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.TFileService = cms.Service('TFileService', fileName = cms.string('hists.root'))

process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')

# ---- GlobalTag (2022 preEE) ----
if isMC:
    process.GlobalTag.globaltag = cms.string('130X_mcRun3_2022_realistic_v5')
else:
    # TODO(Run-3 2022preEE Data): prompt/재처리 GT 확정치로 교체.
    process.GlobalTag.globaltag = cms.string('130X_dataRun3_PromptAnalysis_v1')

# ---- VID + egamma post-reco: modifiedHeepElectronID 생성 (electron veto 용) ----
# energy scale/smearing 보정 off → 후단 correctionlib. mu 채널은 energy 미사용.
runVIDmodules = [
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
    'ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_cff',
]
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False,
                       runVID=True,
                       eleIDModules=runVIDmodules,
                       phoIDModules=[],
                       era='2022-Prompt')

process.load('ZprimeTo4l.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi')
process.load('ZprimeTo4l.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi')
process.load('ZprimeTo4l.Analysis.ResolvedMuCRanalyzerRun3_cfi')
process.resolvedMuCRanalyzerRun3.isMC = cms.bool(isMC)
if not isMC:
    process.resolvedMuCRanalyzerRun3.METfilters = cms.InputTag('TriggerResults', '', 'RECO')

process.evtCounter = cms.EDAnalyzer('SimpleEventCounter')
process.evtCounter.isMC = cms.bool(isMC)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone()
process.hltFilter.throw = cms.bool(True)
process.hltFilter.HLTPaths = cms.vstring('HLT_Mu50_v*', 'HLT_CascadeMu100_v*', 'HLT_HighPtTkMu100_v*')
process.hltFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'HLT')

process.p = cms.Path(
    process.evtCounter+
    process.hltFilter+
    process.ModifiedHEEPIDVarValueMaps+
    process.ModifiedEcalRecHitIsolationScone+
    process.egammaPostRecoSeq+
    process.resolvedMuCRanalyzerRun3
)
```

- [ ] **Step 2: python 구문 검증**

Run:
```bash
cd /u/user/eoyun/4l/25.10.20/CMSSW_13_0_13/src/ZprimeTo4l
python3 -c "import ast; ast.parse(open('Analysis/test/runResolvedMuCRRun3_cfg.py').read()); print('OK')"
```
Expected: `OK`

- [ ] **Step 3: prefiring 흔적 0 확인 (cfg)**

Run:
```bash
grep -in 'prefir' Analysis/test/runResolvedMuCRRun3_cfg.py
```
Expected: 출력 없음.

- [ ] **Step 4: Commit**

```bash
git add Analysis/test/runResolvedMuCRRun3_cfg.py
git commit -m "raRun3: unified MC/Data cfg for Run-3 mu CR analyzer (2022preEE)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: 빌드

**Files:** (없음 — 빌드만)

- [ ] **Step 1: BuildFile glob 확인 (수정 불필요 재확인)**

Run:
```bash
grep 'library file' Analysis/plugins/BuildFile.xml
```
Expected: `<library file="*.cc" name="ZprimeTo4l_Analysis_plugins">` — 새 .cc 자동 포함.

- [ ] **Step 2: scram 빌드 (el8 컨테이너)**

Run (cmssw-el8 컨테이너 안에서, CMSSW 환경 세팅 후):
```bash
cd /u/user/eoyun/4l/25.10.20/CMSSW_13_0_13/src
scram b -j8 2>&1 | tail -40
```
Expected: `ResolvedMuCRanalyzerRun3` 컴파일 성공, 에러 없이 완료. 새 plugin이 새로 추가된 첫 파일이라 인식이 안 되면 `rm -f ../tmp/*/src/ZprimeTo4l/Analysis/plugins/DirCache.json` 후 재빌드(메모리 노트 참조).

- [ ] **Step 3: 컴파일 산출물 확인**

Run:
```bash
scram b 2>&1 | grep -i 'error' || echo "no errors"
```
Expected: `no errors`.

- [ ] **Step 4: Commit (빌드 산출물은 커밋 안 함 — 코드 변경 없으면 skip)**

빌드만이므로 커밋 없음. 빌드 실패 시 원인 파일 수정 후 해당 Task 재방문.

---

## Task 6: 실행 스모크 테스트 + prefiring 부재 확인

**Files:** (없음 — 실행 검증)

- [ ] **Step 1: 소수 이벤트 cmsRun 실행**

먼저 cfg의 `maxEvents`를 임시로 2000으로(또는 명령행 override) 낮춰 실행:
```bash
cd /u/user/eoyun/4l/25.10.20/CMSSW_13_0_13/src/ZprimeTo4l/Analysis/test
cmsRun runResolvedMuCRRun3_cfg.py maxEvents=2000 2>&1 | tail -30
```
> `maxEvents=2000` override가 안 먹으면 cfg에서 `input`을 잠시 2000으로 바꿔 실행 후 되돌린다. 접근 가능한 로컬 Run3Summer22 DY MiniAOD가 없으면 사용자에게 테스트 파일 경로를 요청한다.
Expected: crash 없이 종료, `Events total = ... passed = ...` 요약 출력.

- [ ] **Step 2: 출력 히스토그램 존재 + prefiring 히스토그램 부재 확인**

Run:
```bash
python3 - <<'PY'
import ROOT
f = ROOT.TFile.Open('hists.root')
d = f.Get('resolvedMuCRanalyzerRun3') or f
keys = [k.GetName() for k in (d.GetListOfKeys() if hasattr(d,'GetListOfKeys') else [])]
print('has 4P0F_CR_llll_invM:', '4P0F_CR_llll_invM' in keys)
print('prefire hists present:', any('prefire' in k for k in keys))
PY
```
Expected: `has 4P0F_CR_llll_invM: True`, `prefire hists present: False`.
> TFileService 출력이 `resolvedMuCRanalyzerRun3` 디렉토리로 안 들어가면 파일 최상위에서 키를 조회.

- [ ] **Step 3: Z peak sanity (DY 샘플)**

Run:
```bash
python3 - <<'PY'
import ROOT
f = ROOT.TFile.Open('hists.root')
for hn in ['4P0F_CR_ll1_invM','4P0F_CR_ll2_invM','3P0F_ll_invM']:
    h = f.FindObjectAny(hn)
    if h and h.GetEntries()>0:
        print(hn, 'peak bin center =', round(h.GetBinCenter(h.GetMaximumBin()),1), 'entries=', int(h.GetEntries()))
    else:
        print(hn, 'empty/absent (low stats OK for 2000 evt)')
PY
```
Expected: 통계가 충분하면 peak이 ~91 GeV 근처. 2000 이벤트로 4μ 4P0F 통계가 없으면 `3P0F_ll_invM` 또는 더 많은 이벤트로 재확인(정보용, 실패 기준 아님).

- [ ] **Step 4: 전체 회귀 확인 — 기존 Run-2 파일 무변경**

Run:
```bash
cd /u/user/eoyun/4l/25.10.20/CMSSW_13_0_13/src/ZprimeTo4l
git status --porcelain Analysis/plugins/ResolvedMuCRanalyzer.cc Analysis/python/ResolvedMuCRanalyzer_cfi.py Analysis/test/runResolvedMuCRMC_cfg.py
```
Expected: 빈 출력(3개 원본 무변경).

- [ ] **Step 5: 최종 확인 요약을 사용자에게 보고**

빌드/실행/prefiring-부재/Z-peak 결과를 사용자에게 보고하고, payload TODO 목록(muon SF/PU/RoccoR/FF Run-3 교체)을 다음 작업으로 명시.

---

## Self-Review 결과 (작성자 체크)

- **Spec coverage:** §2(.cc prefiring)→Task2, §3(cfi)→Task3, §4(cfg 통합)→Task4, §5(검증 빌드/실행/grep/Zpeak)→Task5·6, 파일배치 §1→Task1·File Structure. 누락 없음.
- **Placeholder scan:** 코드 블록은 전부 실제 내용. cfi의 `cms.vdouble()` 빈 값은 Task3 Step3에서 원본 값 복사로 채우는 명시적 단계 있음(플랜 placeholder 아님, 의도된 데이터 이관 단계).
- **Type/이름 일관성:** 모듈/클래스명 `ResolvedMuCRanalyzerRun3` 전 Task 일치. cfi 변수명 `resolvedMuCRanalyzerRun3` cfg에서 동일 사용. 파라미터명은 원본 cfi와 diff로 대조(Task3 Step5).
