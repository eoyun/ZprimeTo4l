# Resolved 4-lepton Run-3 Framework — 설계 문서

작성일: 2026-09-17
상태: **설계 확정(사용자 승인). 구현 착수 전 기준 문서.**
대상: resolved `4e`, `2e2mu`, `4mu` (Run 3)
근거 문서: `ResolveAnalysisRun3/resolved_analysis_llm_handoff.md`, 기존 `Analysis/` 코드

---

## 0. 이 문서의 목적

기존 CMS resolved 분석(한 파일에 selection·pairing·weight·histogram·FF가 엉켜 있는
CMSSW analyzer)을 **2단계 구조**로 재설계한 Run-3 전용 framework의 확정 설계다.
사용자가 모든 코드를 직접 이해·검증하는 것이 최우선 목표이므로, 구현 시
**모든 파일에 한 줄 단위 주석 + 모듈 상단 3줄 요약(무엇을/어떻게/무엇에 의존)**을 단다.

### 확정된 핵심 결정

| # | 결정 | 선택 |
|---|---|---|
| 1 | 출력 구조 | **B: ntuple 생산 + 별도 후단 분석** |
| 2 | 후단 backend | **TTreeReader 명시적 event 루프 (C++ 컴파일 실행파일)** |
| 3 | Fake Factor | **적용 우선(기존 payload 재사용), 재측정은 나중 별도 모드** |
| 4 | ID 저장 | **기존 ID 결과(bool) 재현 + ID 입력 scalar 변수 저장** |
| 5 | Config 형식 | **JSON** (nlohmann/json, correctionlib 경유로 이미 가용) |
| 6 | 대상 | **Run 3 전용. Run-2 재현은 목표 아님.** |
| 7 | Correction | **앞단(cmsRun)에서 계산, 보정 전/후 값 모두 저장** |
| 8 | 확장성 | jet-map veto·MET rejection을 나중에 깨지 않고 추가 가능하게 (지금은 미포함) |
| 9 | Prefiring | **완전 제거** (Run-3에서 무시 가능, 도구도 Run-2 전용) |
| 10 | 검증 | **ZZ sample의 Z peak 등 data-driven** (Run-2 대조 아님, 별도 영역 추가 없음) |

---

## 1. 전체 아키텍처

```
[MiniAOD]
   │  cmsRun : ResolvedNtuplizer (EDAnalyzer 1개)
   │   · CMSSW 의존 계산만: modified-HEEP 입력, secondary-GSF association,
   │     muon TuneP/inner/global 정보, modified isolation 재료(raw iso + neighbor),
   │     HLT decision + trigger object 연결, weights(genWeight/PU),
   │     lepton correction 전/후 + up/dn variation
   │   · 느슨한 production skim만 (event 제거 최소화)
   ▼
[Events.root]   event당 1 entry, lepton은 가변길이 배열, 연결은 index 직렬화
   │  후단 C++ 실행파일 (bin/runResolved.cc)  ← test/selection.json
   │
   ├─ NtupleReader     : Events.root → 후단 struct (Event / Electron[] / Muon[])
   ├─ ObjectSelector   : P/F 판정, acceptance, modified iso 재계산, opposite-flavor veto
   ├─ CandidateBuilder : flavor별 pairing (4e reciprocal-GSF→dR, μ closest-dR)
   ├─ RegionSelector   : trigger match, multiplicity, pair-ID, mll/m4l, SR/CR 분류
   ├─ WeightProvider   : weight component 리스트(genWeight·PU·SF·FF) + variation 조합
   ├─ HistogramWriter  : 공통 selection 결과 + weight → 히스토그램
   └─ BackgroundBuilder: prompt subtraction, CR 조합, 최종 template (clipping은 여기서만)
   ▼
[histograms.root]
```

**설계 원칙**
- 각 모듈은 작은 struct + 함수로 시작(거대 class 계층 금지).
- **selection 함수는 히스토그램/파일을 건드리지 않는다** — "통과/불통과 + 후보"만 반환.
  한 함수만 떼어 읽어도 이해·단위검증 가능.
- 같은 물리 컷을 nominal/systematic 코드에 복제하지 않는다(한 곳 정의).
- 모든 물리 컷 값은 `selection.json`에서만 온다(하드코딩 금지, 누락 시 명시적 에러).

### 저장소 배치 (기존 `Analysis/`는 건드리지 않음)

```
ResolveAnalysisRun3/
├── plugins/    ResolvedNtuplizer.cc            (cmsRun EDAnalyzer)
│               BuildFile.xml
├── interface/  NtupleSchema.h                  (branch 정의 = 생산/후단 공유 계약)
│               Event.h  Electron.h  Muon.h     (후단 struct)
│               NtupleReader.h  ObjectSelector.h  CandidateBuilder.h
│               RegionSelector.h  WeightProvider.h  HistogramWriter.h
│               BackgroundBuilder.h
├── src/        위 헤더 구현
├── bin/        runResolved.cc                  (후단 main: 4l 분석)
│               measureFF.cc                    (나중: FF 재측정 모드)
│               BuildFile.xml
├── python/     production 설정 조각 (cfi)
├── test/       runResolvedNtuplizer_cfg.py     (cmsRun 설정)
│               selection.json                  (후단 컷)
├── data/       SF/FF payload (기존 Analysis/data 재사용 또는 복사)
└── docs/       본 설계 문서
```

---

## 2. Ntuple Schema (Events.root) — 되돌릴 수 없는 계약

### 2.1 저장 원칙
- `Events` TTree: event당 1 entry. lepton은 가변길이 `std::vector<...>` branch.
- lepton 간 연결(secondary-GSF, neighbor)은 **컬렉션 index**로 직렬화. 살아있는
  `pat::*Ref`나 메모리 주소를 ntuple에 남기지 않는다. null = `-1`.
- **selection용 값과 mass용 값 분리**, **correction 전/후 모두 저장**.
- precision: 대부분 `float`. `event`는 `ULong64_t`, index는 `int`.
- 각 branch의 단위/의미/null 표현을 `NtupleSchema.h` 주석에 명시.

### 2.2 Event 수준 (스칼라)
| branch | 의미 |
|---|---|
| `run, lumi, event` | event 식별(검증 대조) |
| `era, sampleId` | era/sample 태그 |
| `HLT_*` (bool 여러 개) | 필요한 각 trigger path decision |
| `passMETfilters` 또는 filter별 bool | event filter |
| `nPV, pvNdof, pvZ, pvRho` | PV 정보(muon ID가 PV 의존) |
| `beamSpotX/Y/Z` | neighbor dxy(BS) 재계산용 |
| `genWeight, puTrue` | MC weight, PU true interactions(후단서 PU nominal/up/dn 재평가) |

> prefiring branch 없음(결정 9). MET/jet branch 없음(결정 8, 확장 시 추가).

> **Kinematics source — Particle Flow 사용 안 함.** muon=TuneP(2.4절), electron=ECAL-driven
> + GSF 기반. 기존 `ResolvedEleCRanalyzer.cc` 확인 결과:
> selection/정렬 ET = `pat::Electron::et()` (ECAL transverse energy, PF pt 아님, line 541·583),
> acceptance η = `superCluster()->eta()` (line 547~551),
> mass p4 = `polarP4() * userFloat("ecalTrkEnergyPostCorr") / energy()` (line 731·972~975).
> → 아래 branch가 이 세 가지를 각각 보존한다.

### 2.3 Electron 배열 (index i)
| branch | 의미 |
|---|---|
| `ele_charge[i]` | 전하 |
| `ele_selEt[i]` | **selection용 ET = `et()`** (ECAL-driven+GSF, PF 아님) |
| `ele_eta, ele_phi, ele_etaSC[i]` | GSF reco 방향 + supercluster η(acceptance 기준) |
| `ele_rawPt/eta/phi/energy[i]` | **보정 전** GSF/ECAL reco p4 (`polarP4()`, `energy()`) |
| `ele_corrPt/eta/phi/m[i]` | **보정 후** mass p4 (`polarP4*ecalTrkEnergyPostCorr/energy`) |
| `ele_ecalTrkEnergyPreCorr, ele_ecalTrkEnergyPostCorr[i]` | scale factor 재구성용 |
| `ele_energyScaleUp/Dn, ele_energySmearUp/Dn[i]` | variation p4(앞단 저장) |
| `ele_passModHeep[i]` | 기존 VID/modified-HEEP 결과(재현 기준) |
| `ele_passModHeepLoose[i]` | loose denominator(F 판정) |
| `ele_HoverE, ele_sigmaIeta, ele_dEtaSeed, ele_dPhiIn, ele_dxy, ele_lostHits, ele_modTrkIso, ele_modEcalHcalIso, ele_rho[i]` | ID 입력변수(threshold 재판정, 결정 4) |
| `ele_addGsfIdx[i]` | secondary-GSF가 가리키는 electron index(없으면 -1) — 4e reciprocal pairing 재현 |
| `ele_modHeepBitmap[i]` | `isNonHeepEle` masked bit(F 판정 계약) |

### 2.4 Muon 배열 (index j)
| branch | 의미 |
|---|---|
| `mu_charge[j]` | 전하 |
| `mu_rawTunePpt[j]` | **보정 전** TuneP pt |
| `mu_corrTunePpt[j]` | **보정 후** TuneP pt (Rochester/scale/smear) |
| `mu_rochesterSF[j]` | 적용된 Rochester 계수 |
| `mu_tunePeta, mu_tunePphi[j]` | TuneP 방향 |
| `mu_innerPt[j]` | inner-track pt(neighbor subtraction 재료) |
| `mu_recoPt/eta/phi[j]` | reco(SF lookup + neighbor self 방향) |
| `mu_recoVz[j]` | [cm] muon vertex z (neighbor dz 비교: self) |
| `mu_innerEta/innerPhi/innerVz[j]` | inner-track 방향+vz (neighbor 로 쓰일 때 checkIso 재현) |
| `mu_innerDxyBS[j]` | [cm] inner-track signed dxy wrt beamspot (neighbor 컷, signed) |
| `mu_isHighPt, mu_isTrackerHighPt[j]` | 두 ID 결과(하나로 합치지 않음) |
| `mu_relPtErr, mu_trkLayers, mu_pixelHits, mu_matchedStations, mu_dxy, mu_dz[j]` | ID 입력변수 |
| `mu_trackIso[j]` | raw track iso |
| `mu_neighborIdx[j]` | 뺀 neighbor muon index(없으면 -1) |
| `mu_ptScaleUp/Dn, mu_ptSmearUp/Dn[j]` | variation |

> modified iso는 값 하나만 저장하지 않고 raw iso + 모든 muon inner pt + neighbor 관계를
> 저장해 후단에서 규칙 변경 시 재계산 가능하게 한다(handoff 경고 반영).

### 2.5 Trigger 배열 (index k)
| branch | 의미 |
|---|---|
| `trigObj_pt/eta/phi[k]` | trigger object kinematics |
| `trigObj_filterBits[k]` 또는 `trigObj_pathId[k]` | filter/path 연결(matching threshold 후단 변경 가능하게) |

### 2.6 Metadata (별도 TTree/TH1, job당)
`schemaVersion, sourceCommit, productionConfig, payloadIds, appliedSkim`
+ 정규화용 **skim 이전** `Nevents, sumGenWeights, sumGenWeights2`.

### 2.7 Optional reference (검증 대조용, on/off)
기존식 선택 결과 `channel, P/F, pair index, mll1/mll2, m4l, cut flags` — 대조할 때만 켠다.

### 2.8 확장성 장치
1. **Collection block 패턴**: 모든 컬렉션을 `{prefix}_n` + `{prefix}_<var>[]` 동일 규칙.
   `NtupleReader`는 prefix로 블록을 읽는 제네릭 함수 → jet 추가 = 블록 등록만.
2. **schemaVersion + optional branch**: 분석이 요구하는 branch가 없으면 명시적 에러
   (임의 기본값 금지). 새 블록은 버전으로 구분.
3. **Event-level scalar 자유 확장**: MET(`met_pt/phi`)·jet-veto flag는 나중에 스칼라 추가만.
   후단 `RegionSelector::passEventVeto`는 config `eventVeto`가 비면 no-op.

---

## 3. 후단 분석 흐름 (TTreeReader)

`bin/runResolved.cc` main이 config를 읽고 event 루프를 돈다. 한 바퀴 = 기존 `analyze()` 한 번.

```cpp
Config cfg = loadConfig("selection.json");    // 모든 물리 컷의 유일한 출처
NtupleReader reader(inputFiles, cfg);
HistogramWriter hw(outFile, cfg);

while (reader.next()) {
  Event ev = reader.event();                  // 스칼라 + electrons[] + muons[]

  // 1) event-level 전처리
  double w = WeightProvider::eventWeight(ev, cfg);       // genWeight 부호 · PU
  if (!RegionSelector::passHLT(ev, cfg))        continue;
  if (!RegionSelector::passMETfilters(ev, cfg)) continue;
  if (!RegionSelector::passPV(ev, cfg))         continue;
  // (확장) if (!RegionSelector::passEventVeto(ev, cfg)) continue;  // 지금 no-op

  // 2) object 판정 — 순수 함수 (히스토그램/파일 접근 없음)
  auto sel = ObjectSelector::select(ev, cfg);            // P/F 컬렉션, modified iso 재계산

  // 3) 채널 결정 (multiplicity + opposite-flavor veto)
  Channel ch = RegionSelector::classifyChannel(sel, cfg);
  if (ch == Channel::none) continue;

  // 4) trigger plateau + matching (채널별)
  if (!RegionSelector::passTrigger(ev, sel, ch, cfg)) continue;

  // 5) pairing (4e reciprocal-GSF→dR, μ closest-dR)
  auto cand = CandidateBuilder::build(sel, ch, cfg);

  // 6) pair-level muon 조건 + dilepton mass cut
  if (!RegionSelector::passPairMuonID(cand, sel, ch, cfg)) continue;
  if (!RegionSelector::passDileptonMass(cand, cfg))        continue;   // mll1,mll2 > min

  // 7) SR/CR 분류 → weight 조합 → fill (nominal + 모든 variation)
  Region reg = RegionSelector::classifyRegion(sel, cand, cfg);
  double wTot = w * WeightProvider::sfWeight(cand, sel, cfg)
                  * WeightProvider::ffWeight(reg, cand, sel, cfg);
  hw.fill(reg, cand, ev, wTot, cfg);
}
```

### Pairing 규칙 (기존 `PairingHelper` 재현)
- **4e**: reciprocal secondary-GSF 우선(`a.gsf==addGsf[b] && b.gsf==addGsf[a]`),
  실패 시 3조합 `(01)(23)/(02)(13)/(03)(12)` 중 `min(dR1,dR2)` 최소 선택.
- **μ**: closest-dR (`min(dR1,dR2)` 최소).
- pairing 후 muon pair-ID 실패해도 **다른 pairing으로 자동 재시도하지 않음**(selection 변경 금지).

---

## 4. Config 소유권 (같은 컷 두 곳 복제 금지)

- **`production_cfg.py`** (cmsRun): source, era, GT, producer, **느슨한 skim**, ntuple 출력.
- **`selection.json`** (후단): 모든 최종 물리 컷.

```jsonc
{
  "schemaVersion": 1,
  "dataset": {
    "era": "2023",
    "logicalPD": { "Muon": ["Muon0","Muon1"], "EGamma": ["EGamma0","EGamma1"] },
    "channelPD": { "4e": "EGamma", "2e2mu": "Muon", "4mu": "Muon" },
    "triggerOverlapRemoval": true
  },
  "objects": {
    "electron": { "etMin": 20.0, "ebEtaMax": 1.4442, "eeEtaMin": 1.566, "eeEtaMax": 2.5,
                  "modTrkIsoMax": 5.0, "dxyEB": 0.02, "dxyEE": 0.05, "lostHitsMax": 1 },
    "muon":     { "tunePptMin": 20.0, "etaMax": 2.4, "relPtErrMax": 0.3,
                  "trkLayersMin": 5, "pixelHitsMin": 0, "dxyMax": 0.2, "dzMax": 0.5,
                  "modIsoRelMax": 0.1,
                  "neighbor": { "drMin": 0.01, "drMax": 0.3, "dzMax": 0.2, "dxyBSMax": 0.1 } }
  },
  "trigger": {
    "doubleElectron": { "paths": ["<사용자 조사>"], "leadEtMin": 35.0, "subEtMin": 35.0 },
    "singleMuon":     { "paths": ["<사용자 조사>"], "leadTunePptMin": 52.0 },
    "matchDR": 0.1
  },
  "pairing":  { "electron": "reciprocalGsfThenDR", "muon": "closestDR" },
  "massCuts": { "dileptonMassMin": 1.0, "validationMassMin": 50.0, "signalMassMin": 200.0 },
  "fake":     { "pairDRBoundary": 0.3, "muonRelIsoMax": 0.5, "ffPayload": "data/RMFF.root" },
  "eventVeto": {}
}
```

- 값은 handoff 5~7절/기존 코드에서 확인한 nominal. era별 override 허용.
- **trigger path / GT / era별 payload / 실제 dataset은 사용자가 DAS에서 조사해 채운다**(임의 금지).
- `eventVeto` 비면 no-op → jet/MET 확장 전까지 영향 없음.

---

## 5. 재현할 object selection (AN nominal, 기존 코드 기준)

### Electron (modified HEEP)
ET>20; EB `|etaSC|<1.4442`, EE `1.566<|etaSC|<2.5`; ECAL-driven seed;
H/E EB `<1/E+0.05` EE `<5/E+0.05`; EE `sigma_ieta<0.03`;
track-cluster EB `|dEtaSeed|<0.004 OR |dUin5x5|<0.004` EE `|dEtaSeed|<0.006`;
`|dPhiIn|<0.06`; modified track iso `<5`; modified ECAL+HCALd1 iso
EB `<2+0.03ET+0.28rho`, EE(ET<50) `<2.5+0.28rho`, EE(ET≥50) `<2.5+0.03(ET-50)+0.28rho`;
lost hits `≤1`; `|dxy|` EB`<0.02` EE`<0.05`.
mass p4 = `polarP4*ecalTrkEnergyPostCorr/energy`.

### Muon
`pT(TuneP)>20`, `|eta|<2.4`; global-highPt OR tracker-highPt ID;
modified tracker iso / TuneP pT `<0.1`;
공통 track: TuneP relPtErr`<0.3`, layers`>5`, pixel`>0`, `|dxy|<0.2`, `|dz|<0.5`.
두 ID를 하나의 flag로 합치지 않음. modified iso = raw iso − (조건 만족 neighbor 중
TuneP pT 최대 하나의 inner pT). neighbor: `0.01<dR<0.3`, `|dz|<0.2`, `|dxyBS|<0.1`.

### Resolved SR (채널별)
| 항목 | 4e | 2e2mu | 4mu |
|---|---|---|---|
| Passing mult. | NeP=4 | NeP=2, NmuP=2 | NmuP=4 |
| Opp-flavor veto | highPt μ veto(iso 전 collection) | 정확 P mult. 분류 | modHEEP P e veto |
| Trigger | double e | single μ | single μ |
| Plateau | lead/sublead e ET>35 | lead μ TuneP pT>52 | 동일 |
| Matching | seeded/unseeded leg, dR<0.1 | lead μ dR<0.1 | 동일 |
| Pairing | reciprocal GSF→dR | ee, μμ | closest dR |
| μ-pair 조건 | — | pair 내 global-highPt P μ ≥1 | 각 pair ≥1 |
| Paired mass | 각 pair mll>1 | 동일 | 동일 |
| Final SR | m4l>200 | 동일 | 동일 |

SR에 OS/Z-window/equal-mass/dR>0.3/MET cut 추가 금지. event filter는 유지.

---

## 6. Background regions & Fake Factor

### Passing/Failing (F ≠ !P)
- **Electron F**: acceptance + loose modHEEP 만족 + P 아님. `isNonHeepEle` masked bitmap
  + `|dU|<0.01 OR |dEtaSeed|<0.01` 기준(ID bit 순서도 계약).
- **Muon F**: tracker muon, acceptance, stations≥1, layers>5, pixel>0, `|dxy|<0.2`,
  `|dz|<0.5`, P collection에 없음, modified iso/TuneP pT<0.5. ID 통과·iso 실패도 포함.
  `passID/passIso/passP/passF` 각각 기록.

### 영역
| 영역 | 조건 |
|---|---|
| 4P | SR selection |
| Isolated 3P1F | F와 선택 partner dR>0.3 (single-lepton FF) |
| Isolated 2P2F | 각 F와 partner dR>0.3 (두 FF 곱) |
| Nonisolated 2P2F | 두 F 같은 flavor, 선택된 동일 pair, dR<0.3 (별도 dilepton FF) |

dR=0.3 경계(`<`/`>`) 명시. CR muon pair 조건은 기존처럼 완화(`isGlobalMuon()≥1`) —
SR의 global-highPt 조건을 모든 CR에 복사하지 않음.

### FF (적용 모드)
- 기존 `RMFF.root` TF1을 `WeightProvider`가 읽어 3P1F/2P2F에 적용.
- `A = Data(3P1F·f) − ZZ(3P1F·f)`, `B = Data(2P2F·f1f2)` 를 각각 히스토그램 보존.
- clipping(`max(A−2B,0)+B` 등)은 `BackgroundBuilder`에서 최종 template 만들 때만,
  정책을 config에 명시. raw prediction + uncertainty 보존.
- 재측정은 나중 별도 모드(`bin/measureFF.cc`). 그때 3-lepton tag sample 필요 →
  **production skim이 지금부터 3-lepton event를 남겨야** 재측정 가능(아래 8절).

---

## 7. Systematics

Run-2 대조 없이 처음부터 migration-aware.
- **Weight-only**: PU up/dn, muon ID/ISO/RECO/TRIG SF, electron modHEEP/RECO/TRIG SF, FF up/dn.
- **Kinematic**: muon scale/smear, electron scale/smear — 앞단 저장 up/dn p4로
  pairing/mass cut까지 재평가(candidate/region migration 반영).
- Smearing 재현성: seed를 `event`+object index로 결정론적 고정. nominal/up/dn 독립 난수 금지.
- Ordinary/boosted muon iso SF 이중 적용 금지(분기 구분).

---

## 8. Production skim & 정규화

- **skim은 느슨하게**: 후단에서 평가할 SR/CR/variation 전체의 합집합 포함.
  FF 재측정 대비 **3-lepton event까지 남긴다**(schema는 4l용이어도 event는 보존).
  correction variation에 따른 threshold migration도 포함.
- 정확한 HLT/pT/multiplicity skim 값은 **사용자가 조사 후 확정**(미결정).
- 정규화: job당 skim 이전 `Nevents, sumGenWeights, sumGenWeights2` 기록.
  병합 시 각 job 합 정확히 한 번. 후단 cutflow마다 `Nraw, sumw, sumw2`, 음수 weight 유지.

---

## 9. 검증 전략 (Run-3 data-driven, 별도 영역 추가 없음)

1. **Z peak (ZZ sample)**: ZZ MC의 진짜 Z→ll pair로 dilepton mass가 91 GeV에 서는지.
   기존 영역 활용, 별도 tag 영역 추가하지 않음.
2. **Correction 전/후**: `corr/raw` pt 분포, Z mass resolution 개선 확인(2.3/2.4 branch).
3. **Cutflow sanity**: 단계별 event 수가 물리적으로 타당한지.
4. **Data/MC**: CR(3P1F/2P2F)에서 data vs (prompt ZZ + nonprompt FF prediction).
5. **Config round-trip**: `selection.json`의 mll만 바꿔 후단 결과가 변하고
   ntuple 재생산이 불필요함을 확인(아키텍처 B 이점 검증).

---

## 10. 범위 밖 (이번 미포함, 확장 지점만 확보)

- Merged-electron SR, cleaned-muon SR (handoff 범위 제외).
- Jet-map veto, MET-based rejection (schema `eventVeto`/collection-block로 확장 지점만).
- FF 재측정 (skim만 대비, 코드는 나중).
- Prefiring (완전 제거).
- Run-2 재현 (목표 아님).

---

## 11. 사용자가 채울 미결정 항목 (임의로 채우지 않음)

- 실제 dataset 경로(2022/2023, Muon/EGamma[0/1]), Global Tag, era.
- Trigger path/filter label (Run-3용), plateau 값 재확인.
- Production skim 정확한 HLT/pT/multiplicity threshold.
- FF payload: 기존 `RMFF.root` 재사용 여부 / era 적합성.
- SF/PU payload의 Run-3 버전.

---

## 12. 구현 순서 (승인 시)

1. `NtupleSchema.h` 확정 + `production skim` 정의. 각 후단 cut ↔ branch 매핑 표 작성.
2. `ResolvedNtuplizer` 구현 + 작은 sample 실행. P/F·trigger denominator·3-lepton 보존 확인.
3. `NtupleReader` + 후단 struct. correction 전/후 branch 읽기 검증.
4. `ObjectSelector` + `CandidateBuilder` (nominal 4P). 세 채널 pairing/mass 확인.
5. `RegionSelector` + `HistogramWriter` (nominal 4P SR/CR).
6. 3P1F/2P2F + isolated/nonisolated FF 적용.
7. `WeightProvider` systematics + 정규화 + cutflow.
8. `BackgroundBuilder` (prompt subtraction, template, clipping 정책).
9. ZZ sample 검증(9절). 대표 sample로 용량/시간 측정.
```
