# Resolved μ CR analyzer — Run-3 포팅 (템플릿) 설계

작성일: 2026-09-23
상태: **설계 확정(사용자 승인). 구현 착수 전 기준 문서.**
대상: `ResolvedMuCRanalyzer` → `ResolvedMuCRanalyzerRun3` (4μ resolved CR, **2022 preEE only**)
근거: 기존 `Analysis/` Run-2 UL 코드, `ResolveAnalysisRun3` Run-3 ntuplizer cfg (권위 레퍼런스)

---

## 0. 목적 / 범위

기존 MiniAOD monolithic EDAnalyzer `Analysis/plugins/ResolvedMuCRanalyzer.cc`(Run-2 UL)를
**Run-3(2022 preEE)용으로 직접 포팅**한다. 이 μ 채널을 **템플릿**으로 먼저 완성·검증한 뒤
Ele/EMu 채널을 같은 패턴으로 복제한다(별도 작업).

핵심 변경 3가지 + 채택된 추가 제안 1가지:
1. **Prefiring 완전 제거** (Run-3 무시 가능, 도구도 Run-2 전용)
2. **egammaPostRecoTools: energy 보정 off, VID만 유지** (μ 채널은 electron을 veto ID로만 사용 →
   energy 보정 불필요. correctionlib energy 배선은 Ele/EMu 포팅에서 등장)
3. **cfg/cfi Run-3화** (GT, era modifier, egamma era, trigger/MET filter, payload 경로)
4. **(E2) Data/MC cfg 통합** — `isMC` 플래그 하나로 단일 cfg 파일

### 제외(이번 범위 아님)
- Ele/EMu 채널 포팅
- FF 재측정 (Run-3 payload는 재측정 필요; 이번엔 파라미터만 유지)
- muon momentum 보정 마이그레이션(RoccoR→correctionlib) — 결정상 **유지, 나중에**
- jet-map veto / MET rejection 확장
- 다중 era 스위칭(E1) — 단일 era(2022preEE)라 불필요

---

## 1. 결과물 파일 (모두 `Analysis/` 패키지, `Run3` 접미사)

| 파일 | 기반 | 비고 |
|---|---|---|
| `Analysis/plugins/ResolvedMuCRanalyzerRun3.cc` | `ResolvedMuCRanalyzer.cc` 복사 | `plugins/BuildFile.xml`이 `*.cc` glob → 자동 빌드(수정 불필요) |
| `Analysis/python/ResolvedMuCRanalyzerRun3_cfi.py` | `ResolvedMuCRanalyzer_cfi.py` | 단일 era(2022preEE), payload placeholder+TODO |
| `Analysis/test/runResolvedMuCRRun3_cfg.py` | `runResolvedMuCRMC_cfg.py` + ntuplizer Run-3 cfg | **MC/Data 통합**, `isMC` 플래그 |

기존 helper(`PairingHelper`, `MuonCorrectionHelper`, `RoccoR`)는 **그대로 재사용**(수정 없음).
기존 Run-2 파일 3종은 손대지 않고 그대로 둔다.

EDAnalyzer 클래스 등록 이름: 새 클래스명 `ResolvedMuCRanalyzerRun3`,
`DEFINE_FWK_MODULE(ResolvedMuCRanalyzerRun3)`. cfi의 모듈 타입도 동일 문자열.

---

## 2. `.cc` 변경 — prefiring 완전 제거 (그 외 로직 100% 동일)

`ResolvedMuCRanalyzer.cc`(1400줄) 대비 삭제 지점(원본 라인 기준):

| 위치 | 원본 | 조치 |
|---|---|---|
| 멤버 선언 | L60–62 `prefweight_token / prefweightUp_token_ / prefweightDn_token_` | 삭제 |
| 생성자 초기화 | L124–126 위 토큰 3개 초기화 | 삭제 |
| analyze 지역변수 | L426–428 `float prefireNo/Up/Dn = 1.` | 삭제 |
| analyze getByToken | L431–441 prefiring weight 3개 읽기 | 삭제 |
| weight 조합 | L447 `aWeight = prefireNo*mcweight/std::abs(mcweight);` | → `aWeight = mcweight/std::abs(mcweight);` |
| 히스토그램 부킹 | L378–379 `4P0F_CR_llll_invM_prefireUp/prefireDn` | 삭제 |
| 히스토그램 채우기 | L889–890 prefireUp/Dn fill | 삭제 |

- 클래스명/생성자/소멸자/`DEFINE_FWK_MODULE` 를 `ResolvedMuCRanalyzerRun3`로 변경.
- 그 외(muon 선택, isolation, boosted, pairing `pair4M`, 4P0F/3P0F/2P1F/3P1F/2P2F, FF,
  muon scale/smear·ID/ISO/boostISO/trig/reco SF·PU systematic)는 **동일 유지**.
- muon momentum 보정(RoccoR `nominalMC/Data`, `altScale`, `smear`)·SF 로직 **그대로**.
  Run-3 payload가 아직 없다는 점은 cfi의 TODO로만 표시(코드는 payload-무관).

---

## 3. `cfi` 변경 (`ResolvedMuCRanalyzerRun3_cfi.py`)

Run-2 cfi는 4개 era 변형을 한 파일에 담았으나, Run-3는 **단일 base 모듈(2022preEE)** 로 단순화.

- 모듈 타입: `cms.EDAnalyzer('ResolvedMuCRanalyzerRun3', ...)`.
- **trigList (Run-3 μ)**: `HLT_Mu50_v*`, `HLT_CascadeMu100_v*`, `HLT_HighPtTkMu100_v*`
  (Run-2의 `HLT_TkMu50_v*`는 Run-3 부재 → 교체).
- **METfilterList (Run-3)**: `Flag_goodVertices`, `Flag_globalSuperTightHalo2016Filter`,
  `Flag_EcalDeadCellTriggerPrimitiveFilter`, `Flag_BadPFMuonFilter`, `Flag_BadPFMuonDzFilter`,
  `Flag_hfNoisyHitsFilter`, `Flag_eeBadScFilter`, `Flag_ecalBadCalibFilter`
  (Run-2의 HBHENoise/HBHENoiseIso 제거).
- **payload 경로 (repo에 Run-3 버전 없음 → Run-2 UL을 placeholder로 두고 `# TODO(Run-3 2022preEE)` 주석)**:
  - `rochesterPath` — RoccoR(Run-3 공식 부재; 유지, TODO)
  - `triggerSF / muonIdIsoSFpath / muonBoostIsoSFpath / muonRecoSFpath` — Run-3 Muon POG SF 필요, TODO
  - `PUrwgt` + `PUname` — Run-3 LUM POG PU json 필요, TODO
  - `FFpath` (`RMFF_run2.root`) — Run-3 FF 재측정 필요, TODO
- threshold(`ptThres`, `ptThresTrig`, `ffSystCL95`)·`muScaleBias/muSmearFactors/muSmearParams` 값은
  일단 Run-2 값 유지(고pT scale GE-method Run-3 재검증은 `# TODO` 주석만).

> payload는 코드가 `edm::FileInPath`로만 참조하므로, 파일 교체만으로 Run-3 반영 가능.
> 이번 목표는 "컴파일·실행되는 Run-3 골격 + 명확한 payload TODO"이다.

---

## 4. `cfg` 변경 (`runResolvedMuCRRun3_cfg.py`) — MC/Data 통합(E2)

ntuplizer Run-3 cfg(`26.09.21_22preEE/runResolvedNtuplizer_run3MC_cfg.py`)를 권위 레퍼런스로 정렬.

- 파일 상단 `isMC = True` 스위치(주석으로 Data 전환법 명시). 아래가 `isMC`에 분기:
  - `process.GlobalTag.globaltag`: MC `130X_mcRun3_2022_realistic_v5` /
    Data 2022preEE prompt GT(예: `130X_dataRun3_PromptAnalysis_v1` — 실제 값은 실행 시 확정, TODO 주석).
  - analyzer `isMC` 파라미터, `SimpleEventCounter.isMC`.
  - MET filter process 태그(`'PAT'` MC vs `'RECO'` Data) 확인 주석.
- **Run-3 era modifier 필수**: `from Configuration.Eras.Era_Run3_cff import Run3` →
  `process = cms.Process('ResolvedMuCRanalysisRun3', Run3)`.
- GlobalTag load: `Configuration.StandardSequences.FrontierConditions_GlobalTag_cff`
  (Run-2 cfg의 `_condDBv2_cff` 아님).
- **egamma VID-only**:
  ```python
  setupEgammaPostRecoSeq(process,
      runEnergyCorrections=False,   # energy 보정 off → 후단 correctionlib (μ 채널은 미사용)
      runVID=True,                  # electron veto용 modifiedHeepElectronID 생성
      eleIDModules=runVIDmodules, phoIDModules=[],
      era='2022-Prompt')
  ```
- **prefiring 제거**: prefiring producer 로드/경로 없음(analyzer도 토큰 삭제되어 정합).
- hltFilter HLTPaths: `HLT_Mu50_v*`, `HLT_CascadeMu100_v*`, `HLT_HighPtTkMu100_v*`.
- Path: `evtCounter + hltFilter + ModifiedHEEPIDVarValueMaps + ModifiedEcalRecHitIsolationScone
  + egammaPostRecoSeq + resolvedMuCRanalyzerRun3`
  (μ 채널은 electron veto만 하므로 `modifiedHEEPIDVarValueMaps2nd`/`seedGainEle` 불필요 —
  단, veto가 `modifiedHeepElectronID` embed 값을 읽으면 VID 시퀀스로 충분).
- 입력: Run3Summer22 MiniAODSIM 예시 파일(로컬 테스트용, 주석).

---

## 5. 검증 기준 (완료 판정)

1. **빌드**: `scram b`가 el8 컨테이너에서 새 plugin 포함해 통과(신규 파일 추가 시 필요하면
   `DirCache.json` 제거 후 재빌드).
2. **실행**: `cmsRun runResolvedMuCRRun3_cfg.py`가 Run3Summer22 DY MiniAODSIM 소수 이벤트에서
   crash 없이 `hists.root` 생성.
3. **prefiring 흔적 0**: 새 `.cc`/`cfi`/`cfg`에 `prefir`/`Prefir` 문자열 grep 결과 없음.
   `4P0F_CR_llll_invM_prefireUp/Dn` 히스토그램 미생성.
4. **물리 sanity**: `4P0F_CR_ll1_invM`/`ll2_invM`에서 Z peak(~91 GeV) 확인(DY 샘플).
5. **회귀 없음**: 기존 Run-2 3파일 무변경(git status).

---

## 6. 후속(다음 작업, 이번 아님)
- Ele/EMu 채널 동일 패턴 포팅 (여기서 egamma energy correctionlib 본격 배선 + `seedGainEle`).
- Run-3 payload 채우기: muon SF/PU/RoccoR, FF 재측정.
- E1(다중 era), E3(muon correctionlib), E4(GE-method 재검증), E5(UL data 정리).
