# Resolved 4-lepton framework: LLM 전달용 설계 검토 문서

작성일: 2026-09-17  
상태: **검토용 제안. 구현 완료 또는 최종 설계 확정을 의미하지 않는다.**

## 0. 이 문서를 받은 LLM에게

목표는 기존 CMS 분석의 resolved channel을 바탕으로 다음 분석에 사용할 framework 구조를 결정하는 것이다. 먼저 이 문서와 기존 코드를 읽고 설계의 타당성, 누락된 정보, 구현 범위를 검토하라. 사용자가 이 설계를 그대로 채택할지 수정할지 판단할 수 있도록 설명하라. 검토 요청을 저장소 전면 개편 지시로 해석하지 말라.

다음 세 가지를 구분하라.

1. **기존 동작:** 지정한 master snapshot과 AN에서 확인한 내용.
2. **설계 제안:** ntuple 생산과 후단 분석을 분리하는 권장 구조.
3. **미결정 사항:** production skim, 후단 실행 방식, schema 세부사항 등 사용자가 검토할 항목.

코드를 구현하게 되면 필요한 최소 수정부터 진행하라. 누락된 branch/payload를 임의의 기본값으로 대체하거나 오류를 NaN 등으로 숨기지 말라. 물리 selection 변경을 단순 refactoring에 섞지 말고, weighted/unweighted 값을 함께 기록하라.

## 1. 범위와 근거

### 사용자 요구

- 우선 대상은 resolved `4e`, `2e2mu`, `4mu`이다.
- Merged-electron SR 및 cleaned-muon SR 구현은 이번 범위에서 제외한다.
- 사용자가 dataset, Global Tag 및 기본적인 환경/era 설정을 직접 조사한다. 이를 임의로 채우지 않는다.
- 분석 컷은 설정 파일로 분리한다. **Dilepton mass cut도 반드시 포함한다.**
- `cmsRun`에서 histogram을 직접 만들지, ntuple을 생산하고 별도 analyzer로 분석할지 비교해서 결정한다.
- 기존 알고리즘과 선택을 먼저 재현하고, 다음 분석을 위한 변경 효과는 구분해서 평가한다.

### 검토 기준

- Repository: [eoyun/ZprimeTo4l](https://github.com/eoyun/ZprimeTo4l)
- Branch: `master`
- 확인한 snapshot: `8c7d9c2ac352c376aad4dd166c9c6d8bcaff5d19`
- Snapshot commit 날짜: 2025-07-14
- 분석 노트: **AN-2024/032 v9, 2025-03-03**, “Search for a heavy resonance decaying into four-lepton final states for boosted and resolved regimes in proton-proton collisions at 13 TeV”.
- 핵심 AN 범위: Sections 4–5, Section 6.1, Tables 9–12.

이 문서의 수치는 Run-2 분석 재현 기준이다. Run-3 trigger plateau, correction payload, GT로 검증 없이 전용하지 않는다. 새 master를 사용하면 위 snapshot과의 차이를 먼저 확인한다. 이 검토에서 CMSSW build나 실제 event 처리는 수행하지 않았다.

## 2. 기존 cmsRun 출력: histogram 중심 + 일부 목적별 TTree

기존 resolved 분석은 MiniAOD를 직접 읽는 CMSSW C++ analyzer 안에서 selection, pairing, weight 계산, histogram filling을 수행한다.

예를 들어 `Analysis/test/runMuMC_cfg.py`는 다음 출력을 설정한다.

```python
process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string("hists_mu.root")
)
```

| 출력 | 기존 내용 |
|---|---|
| Histograms | Cutflow, SR/CR kinematics, FF-weighted 분포, systematic variations |
| `numerTree`, `denomTree` | Electron/muon fake-factor 측정용 변수 |
| `tree_3P1F`, `tree_2P2F` | Failing lepton과 CR 관련 변수 |
| `tree_2P2F_dr03` | 가까운 fake dilepton 관련 변수 |

위 TTree들은 `ResolvedEleCRanalyzer.cc`, `ResolvedMuCRanalyzer.cc`에 있다. 확인한 `ResolvedEMuCRanalyzer.cc`에는 TTree 생성 코드가 없다. 현재 출력은 **세 resolved channel 전체를 event/object 수준에서 재선택하는 공통 ntuple이 아니다.** 기존 CR tree만으로 모든 SR/CR selection과 pairing을 다시 실행할 수 있다고 가정하지 않는다.

## 3. 출력 구조 선택

| 검토 항목 | A. MiniAOD → cmsRun → histogram | B. MiniAOD → cmsRun → ntuple → 후단 분석 → histogram |
|---|---|---|
| 기존 코드 수정 범위 | 상대적으로 적음 | Producer, schema, reader가 필요 |
| mll/m4l 컷 변경 | 대체로 MiniAOD 재처리 | 저장 acceptance 안에서 후단만 재실행 |
| Pairing 변경 | MiniAOD 재처리 | 모든 필요한 object와 연결 정보가 있으면 가능 |
| ID/iso threshold 변경 | MiniAOD 재처리 | 입력 변수, 관련 후보, 의존 관계까지 저장했으면 가능 |
| 새 분포 또는 binning | 필요한 분포가 없으면 재처리 | 저장된 변수로 재생성 |
| Track/recHit 기반 변수 계산법 변경 | MiniAOD 재처리 | 원자료가 없으면 역시 MiniAOD 재처리 |
| 저장량의 주된 의존성 | Histogram 개수와 bin 수 | Event 수, object 수, branch 구성 |
| 검증 부담 | 기존 selection 유지 확인 | 기존 결과와의 동등성 + schema/skim 정보 손실 확인 |

**권장안은 B이다.** Selection을 조정할 때마다 MiniAOD를 다시 처리하는 의존성을 줄일 수 있기 때문이다. 단, ntuple에는 후단에서 필요한 값과 event가 실제로 남아 있어야 한다. 용량, 처리 시간, 속도 향상 수치는 아직 측정하지 않았다.

작은 대표 sample에서 두 방식을 비교할 때 `events/s`, peak memory, compressed bytes/event, total output size를 측정한다. 측정 없이 성능 수치를 단정하지 않는다.

## 4. 권장안의 단계별 역할

### 4.1 cmsRun: CMSSW 의존 변수 계산 + ntuple 생산

- MiniAOD 및 필요한 EventSetup 읽기.
- Electron modified-HEEP 입력 변수와 secondary-track association 계산.
- Muon ID 기준값, TuneP/inner/reco 정보, isolation 및 neighbor 관계 계산/추출.
- Correction 결과 및 후단 재계산에 필요한 입력 저장.
- HLT decision, trigger object/필터 연결, event filter 정보 저장.
- 명시적인 production skim만 적용.
- `Events` TTree, production cutflow, normalization metadata, 소수의 진단 histogram 출력.

Resolved-only라도 **ModifiedHEEP의 secondary-track 및 modified-isolation 계산은 필요하다.** Merged SR를 제거하는 것과 이 의존성을 제거하는 것은 다르다. FF 측정 영역의 원래 merged-object veto까지 재현할 경우에도 관련 dependency를 따로 확인한다.

### 4.2 후단: C++/ROOT 분석

- 저장된 변수로 object P/F 판정.
- 채널 multiplicity와 opposite-flavor veto.
- Trigger plateau 및 matching selection.
- Pairing과 pair-level muon 조건.
- **Dilepton mass 및 four-lepton mass cut.**
- SR/CR 분류, SF/FF 적용, systematic variation.
- Physics histogram, cutflow, 선택 event 목록/진단 출력.

후단은 ROOT `TTreeReader` 기반 C++ executable 또는 `RDataFrame` 기반 분석 중 선택한다. 일반 TTree를 읽기 위해 반드시 두 번째 `cmsRun EDAnalyzer`를 만들 필요는 없다. 두 backend를 동시에 구현할 필요도 없다.

### 4.3 정보 보존의 한계

**Ntuple 생산 단계에서 제거한 event/object는 후단 cfg를 바꿔도 복구할 수 없다.**

- 생산 시 `mll > 1 GeV`로 잘랐다면 후단에서 `mll > 0.5 GeV`로 완화할 수 없다.
- 생산 시 `4P`만 저장하면 `3P1F`, `2P2F`를 만들 수 없다.
- 생산 시 lepton 4개 이상만 요구하면 FF 측정용 3-lepton sample을 잃을 수 있다.
- 최종 네 lepton만 저장하면 추가-lepton veto 및 다른 pairing/candidate 정책을 재현하지 못할 수 있다.
- ID bool만 저장하면 ID 자체의 threshold 변경은 불가능하다.
- Modified isolation 값 하나만 저장하면 neighbor 선택 규칙 변경은 불가능하다.
- HLT를 생산 단계에서 잘랐다면 후단에서 저장하지 않은 trigger stream으로 확장할 수 없다.

따라서 **production selection은 후단에서 평가할 SR/CR/variation 전체의 합집합을 포함**해야 한다. 구체적인 pT, multiplicity, HLT skim 값은 미결정이다. Skim은 correction variation에 따른 threshold migration도 포함해야 한다.

## 5. 재현할 object selection

아래 표는 AN의 nominal 설명이다. 코드의 정확한 경계 비교와 차이는 10절에서 별도로 관리한다.

### Electron

| 항목 | 기준 |
|---|---|
| ID | Modified HEEP |
| ET | > 20 GeV |
| Acceptance | EB: abs(etaSC) < 1.4442; EE: 1.566 < abs(etaSC) < 2.5 |
| Seed | ECAL-driven |
| H/E | EB < 1/E + 0.05; EE < 5/E + 0.05 |
| Shower shape | E2x5/E5x5 요구 없음; EE sigma_ieta_ieta < 0.03 |
| Track-cluster matching | EB: abs(deltaEtaSeed)<0.004 OR abs(deltaUin5x5)<0.004; EE: abs(deltaEtaSeed)<0.006 |
| abs(deltaPhiIn) | < 0.06 |
| Modified track iso | < 5 GeV |
| Modified ECAL + HCAL depth-1 iso, EB | < 2 + 0.03 ET + 0.28 rho |
| 동일 iso, EE ET<50 | < 2.5 + 0.28 rho |
| 동일 iso, EE ET>=50 | < 2.5 + 0.03(ET-50) + 0.28 rho |
| Inner lost hits | <= 1 |
| abs(dxy) | EB < 0.02 cm; EE < 0.05 cm |

CMSSW 쪽은 기존 ModifiedHEEP producer/VID 결과를 재사용한다. 후단에서 ID를 재판정하려면 해당 scalar 입력과 적용한 ID 정의/version을 저장하고 기존 VID 결과와 대조한다. Secondary-track 재탐색이나 recHit 재계산까지 가능하다는 뜻은 아니다.

기존 mass 계산은 `polarP4() * ecalTrkEnergyPostCorr / energy()`를 사용한다. Selection용 ET와 mass용 corrected p4를 구분해서 저장한다.

### Muon

- `pT(TuneP)>20 GeV`, `abs(eta)<2.4`.
- Global high-pT OR tracker high-pT ID.
- Modified tracker isolation / TuneP pT < 0.1.
- ID는 기존 `muon::isHighPtMuon`, `muon::isTrackerHighPtMuon` 기준을 따른다. 두 ID를 하나의 global/tracker reconstruction flag로 대체하지 않는다.
- 공통 track 조건: TuneP relative pT error <0.3, tracker layers>5, valid pixel hits>0, abs(dxy)<0.2 cm, abs(dz)<0.5 cm. Station/hit 조건은 두 ID가 다르므로 기존 구현을 재사용한다.

기존 modified isolation은 조건을 만족하는 neighbor 후보 중 TuneP pT가 가장 높은 **하나**의 inner-track pT를 뺀다.

```text
isoModified = trackIso - selectedNeighbor.innerTrack.pt
passIso = isoModified / tuneP_pt < 0.1
```

Neighbor가 없으면 subtraction은 없다. AN의 neighbor 조건은 `0.01 < deltaR < 0.3`, `abs(deltaZ)<0.2 cm`, `abs(dxyBS)<0.1 cm`이다. Raw iso, subtraction, neighbor index, modified iso를 함께 저장한다. ID/P-F 정의를 바꾸면 neighbor 후보 집합도 바뀔 수 있으므로, 재계산 범위를 정하고 필요한 연결/track 변수를 보존한다.

## 6. Resolved SR selection 및 pairing

| 항목 | 4e | 2e2mu | 4mu |
|---|---|---|---|
| Passing multiplicity | NeP=4 | NeP=2, NmuP=2 | NmuP=4 |
| Opposite-flavor veto | High-pT ID muon veto; 기존 코드는 isolation 전 collection | 정확한 P multiplicity로 분류 | Modified-HEEP P electron veto |
| Trigger | Double electron | Single muon | Single muon |
| Plateau | Leading/subleading e ET>35 GeV; 2018은 >28 | Leading passing mu TuneP pT>52 GeV | 동일 |
| Matching | Leading e seeded / subleading e unseeded leg, 각각 deltaR<0.1 | Leading passing mu deltaR<0.1 | 동일 |
| Pairing | Reciprocal secondary GSF association 우선, 이후 closest deltaR | ee 및 mumu | Closest deltaR |
| Muon pair 조건 | 해당 없음 | Pair 안에 global high-pT P muon >=1 | 각 pair마다 >=1 |
| Paired mass | 두 pair 각각 mll>1 GeV | 동일 | 동일 |
| Final SR mass | m4l>200 GeV | 동일 | 동일 |

Passing multiplicity와 loose multiplicity는 다르다. 기존 SR에 없던 `Nloose==4` 또는 추가 F veto를 임의로 추가하지 않는다. PV와 event noise/filter 조건도 유지한다. MET threshold가 없는 SR에서도 기존 event filters는 적용한다.

Pairing 세부 규칙:

1. 4e의 reciprocal association: `a.gsfTrack == b.additionalGsfTrack && b.gsfTrack == a.additionalGsfTrack`.
2. Fallback의 세 pairing: `(01)(23)`, `(02)(13)`, `(03)(12)`.
3. `min(deltaR(pair1), deltaR(pair2))`가 가장 작은 pairing 선택.
4. Pairing 후 muon pair ID 조건 검사. 실패 시 다른 pairing으로 자동 재시도하면 기존 selection 변경이다.

`sum(deltaR)`, `abs(m1-m2)`, Z-mass proximity 최소화로 대체하지 않는다. Tie와 object 순서는 재현 검증 항목이다. SR에 OS, Z window/veto, equal-mass, deltaR>0.3, MET threshold를 추가하지 않는다.

## 7. Dilepton mass cut: 위치와 설정 계약

### 7.1 기존 위치

Pairing 뒤 histogram filling 직전에 다음 조건이 반복된다.

```cpp
if (m4l > 50. && lvecA1.M() > 1. && lvecA2.M() > 1.) {
    // histogram filling
}
```

| 코드 | 4P 조건 위치, 기준 snapshot |
|---|---:|
| ResolvedEleCRanalyzer.cc | 1045행 |
| ResolvedEMuCRanalyzer.cc | 897행 |
| ResolvedMuCRanalyzer.cc | 962행 |

3P1F/2P2F 관련 분기에도 반복된다. **선택된 두 pair 각각의 mass cut이며 가능한 모든 pair에 거는 cut은 아니다.** `m4l>50`은 기존 histogram 수집 범위, AN 최종 SR은 `m4l>200`이다.

### 7.2 설정 파일의 소유권

- `production_cfg.py`: CMSSW source, era, producer, production skim, ntuple 출력.
- 후단 `selection.json` 등: 최종 P/F, trigger plateau, pairing 정책, mll/m4l, SR/CR 컷.
- JSON은 제안 형식이며 YAML/Python 등으로 바꿀 수 있다. **같은 physics cut을 두 파일에 독립적으로 복제하지 않는다.**
- 아래 config key는 신규 설계 이름이다. 기존 코드가 자동으로 읽는다고 가정하지 않는다.

후단 설정 예시 (부분 schema):

```json
{
  "massCuts": {
    "dileptonMassMin": 1.0,
    "validationMassMin": 50.0,
    "signalMassMin": 200.0
  },
  "triggerSelection": {
    "doubleElectronEtMin": 35.0,
    "singleMuonPtMin": 52.0,
    "matchDR": 0.1
  },
  "fakeSelection": {
    "pairDRBoundary": 0.3,
    "muonRelativeIsoMax": 0.5
  }
}
```

2018 재현 시 `doubleElectronEtMin=28.0`. Trigger menu/filter label 및 기타 object cut은 별도 schema 항목으로 명시적으로 공급한다. Electron/muon acceptance, isolation threshold, ID working point도 hard-code하지 않고 정의의 출처와 함께 설정한다.

공통 C++ predicate의 의미:

```cpp
const bool passDileptonMass =
    mll1 > cfg.massCuts.dileptonMassMin &&
    mll2 > cfg.massCuts.dileptonMassMin;

const bool passValidationMass =
    passDileptonMass && m4l > cfg.massCuts.validationMassMin;

const bool passSignalMass =
    passDileptonMass && m4l > cfg.massCuts.signalMassMin;
```

이는 **mass 판정만** 의미하며, 모든 selection을 통과한 SR flag는 아니다. `passValidationMass`는 SR과 겹치는 inclusive 조건이다. 배타적인 low-mass 검증 영역이 필요하면 `passValidationMass && !passSignalMass`로 별도 정의한다. 4P/3P1F/2P2F의 mass 판정은 공통 함수로 묶는다.

Histogram 직접 출력 방식을 선택한다면 동일 설정을 CMSSW analyzer가 읽는다:

```python
massCuts = cms.PSet(
    dileptonMassMin=cms.double(1.),
    validationMassMin=cms.double(50.),
    signalMassMin=cms.double(200.),
)
```

```cpp
const auto cuts = iConfig.getParameter<edm::ParameterSet>("massCuts");
dileptonMassMin_ = cuts.getParameter<double>("dileptonMassMin");
validationMassMin_ = cuts.getParameter<double>("validationMassMin");
signalMassMin_ = cuts.getParameter<double>("signalMassMin");
```

위 코드는 parameter를 analyzer PSet에 연결하고 C++ member/description까지 정의해야 동작하는 연결 예시다. Ntuple 권장안에서는 이 최종 mass selection을 생산 단계의 event 제거 조건으로 사용하지 않는다.

## 8. Background regions와 fake factors

### 8.1 Passing/Failing

F는 단순 `!P`가 아니라 loose denominator 안의 non-P object이다.

- Electron F: acceptance와 loose modified-HEEP 조건을 만족하고 P가 아님. 기존 `ElectronSystematicsHelper::isNonHeepEle`의 masked bitmap 및 `abs(deltaU)<0.01 OR abs(deltaEtaSeed)<0.01`를 기준으로 확인한다. ID bit 순서도 계약의 일부다.
- Muon F, 기존 helper: tracker muon, acceptance, matched stations>=1, tracker layers>5, valid pixel hits>0, abs(dxy)<0.2 cm, abs(dz)<0.5 cm, P collection에 없음, modified iso/TuneP pT<0.5.
- Muon F에는 ID 통과/iso 실패인 object도 들어갈 수 있다. `passID`, `passIso`, `passP`, `passF`를 각각 기록한다.

### 8.2 영역

| 영역 | 조건 | 적용 |
|---|---|---|
| 4P | SR selection | Data와 prompt ZZ + nonprompt prediction 비교 |
| Isolated 3P1F | F와 선택된 partner의 deltaR>0.3 | Single-lepton FF |
| Isolated 2P2F | 각 F와 선택된 partner의 deltaR>0.3 | 두 FF의 곱 |
| Nonisolated 2P2F | 두 F가 같은 flavor, 선택된 동일 pair, deltaR<0.3 | 별도 dilepton FF |

DeltaR=0.3의 `<`/`>` 경계 처리는 명시한다. 이를 SR의 deltaR cut으로 오해하지 않는다. Mixed-flavor CR은 trigger anchor가 되는 passing muon 조건 때문에 모든 P/F flavor 조합을 기계적으로 대칭 추가할 수 없다. 실제 조합별 기존 branch를 확인한다.

CR muon pair 조건은 기존 코드에서 `isGlobalMuon()` >=1로 완화되는 부분이 있다. SR의 global high-pT P 조건을 모든 CR에 복사하지 않는다.

### 8.3 FF 정의와 계산

- Single-lepton transfer factor `f=Npass/Nfail`를 사용한다. 이미 transfer factor인 값을 다시 `f/(1-f)`로 바꾸지 않는다.
- 측정 CR의 prompt contamination subtraction과 적용 CR의 subtraction을 구분한다.
- Nonisolated pair는 `f1*f2`가 아니라 별도 `Fll`을 사용한다.
- AN v9의 nonisolated pair-pT fit: ee는 exponential, mumu는 constant. Payload/계수는 외부 입력.
- FF를 기존 payload로 사용할지 새로 측정할지 미결정이다. 측정까지 한다면 3-lepton 및 서로 다른 flavor의 trigger tag sample을 production에 포함해야 한다. Z-tag OS/mass window 등의 측정용 컷은 SR 컷과 구분한다.

Isolated prediction의 signed 결합 구조:

```text
A = Data(3P1F weighted by f) - ZZ(3P1F weighted by f)
B = Data(2P2F weighted by f1*f2)
N_nonprompt_raw = A - B
```

기존 `runRMFF.cc`에는 bin별 `max(A-2B,0)+B` 처리가 있다. 단순 `A-B`와 음수 bin에서 다르다. 새 framework는 A, B, signed prediction과 uncertainty를 보존하고, 통계 입력용 clipping 정책은 마지막에 명시한다. 다른 채널의 후처리는 각 macro를 확인해 동일하다고 가정하지 않는다.

## 9. Ntuple schema와 공통 코드 구조 제안

### 9.1 저장 원칙

- `Events`: event당 한 entry, electron/muon collection은 가변 길이 배열.
- Electron/muon 연결은 collection index로 직렬화한다. 메모리 주소나 살아 있는 `pat::*Ref`를 일반 ntuple에 그대로 의존시키지 않는다.
- CMSSW producer 내부에서는 기존 ref 기반 helper를 활용하고, 출력 시 index 관계로 변환한다.
- 후단에서 재계산 가능한 값과 CMSSW 생산에 고정된 값을 구분한다.
- 저장 값의 단위, p4 정의, ID bit 의미, null association 표현을 schema에 적는다.

| 범주 | 필요한 정보 |
|---|---|
| Event | run/lumi/event, era/sample 식별, HLT decisions, event filter flags, PV 관련 정보 |
| Electron | original index, charge, selection ET, reco direction, raw/nominal/varied p4, etaSC, ID 입력/bitmap, modified iso, secondary GSF 관계 |
| Muon | original index, charge, reco/TuneP/inner-track 정보, global/tracker ID 및 입력, raw/modified iso, neighbor 정보, nominal/varied p4 |
| Trigger | Object kinematics, path/filter association 또는 재분석 범위를 만족하는 모든 object별 matching 정보 |
| Weights | genWeight, PU 및 필요한 calibration inputs, named correction components와 variations |
| Optional reference | 기존 선택의 channel/P-F/pair indices, masses, cut flags; 재현 비교용 |
| Metadata | Schema version, source commit, production config, payload 식별, 적용한 skim |

Trigger matching threshold를 나중에 바꾸려면 `matched=true/false` 한 값만으로는 부족하다. Object/filter label과 거리 또는 원래 trigger object 정보를 저장해야 한다. Muon station 조건이나 isolation neighbor 집합까지 바꾸려면 관련 입력을 추가로 저장해야 한다. “모든 cut 변경 가능”을 무조건 보장하지 않는다.

### 9.2 정규화

Job별로 production skim 이전의 `Nevents`, `sumGenWeights`, `sumGenWeights2`를 기록한다. 저장된 tree entry의 합으로 전체 sample 정규화를 대체하지 않는다. 병합 시 각 job의 합을 정확히 한 번 더하고 재시도 파일의 중복 합산을 피한다.

후단은 cutflow마다 `Nraw`, `sumw`, `sumw2`를 저장한다. 음수 MC weight를 유지한다. Event-level SF/FF를 하나의 total weight로만 덮어쓰지 않는다.

### 9.3 코드 책임 분리

| 모듈 이름, 제안 | 역할 |
|---|---|
| `ResolvedNtuplizer` | CMSSW 입력/계산/생산 skim/serialization |
| `NtupleReader` | Schema를 읽어 후단 event/object 구조에 연결 |
| `ObjectSelector` | P/F, acceptance, opposite-flavor veto collection |
| `CandidateBuilder` | Flavor별 pairing, indices, kinematics |
| `RegionSelector` | Trigger, multiplicity, pair ID, mass, SR/CR 판정 |
| `WeightProvider` | SF/FF와 variation 조합 |
| `HistogramWriter` | 공통 selection 결과와 weight를 출력 |
| `BackgroundBuilder` | Prompt subtraction, CR 조합, 최종 templates |

모듈마다 거대한 class hierarchy를 만들 필요는 없다. 작은 struct와 함수로 시작한다. Selection 함수는 histogram filling이나 파일 접근을 하지 않게 한다. 같은 cut을 nominal과 systematic code에 따로 복사하지 않는다.

## 10. Systematic 및 기존 코드와의 차이

### Systematic 권장 처리

- Weight-only: selection은 유지하고 해당 weight component를 변경.
- Kinematic: 저장한 variation으로 selection, threshold, candidate/region migration을 재평가할 수 있게 설계.
- 기존 코드는 nominal selection을 고정하고 shifted mass를 채우는 부분이 있으므로, 먼저 기존 방식 재현 결과를 만들고 migration-aware 결과와 비교.
- Smearing 재현성은 event/object 및 variation에 대해 정의한다. nominal/up/down에서 독립적인 난수를 무계획하게 사용하지 않는다.
- Muon nominal은 기존 TuneP 및 해당 MC smearing 흐름을 확인한다. Helper에 Rochester 함수가 있다는 이유로 추가 적용하지 않는다.
- Ordinary/boosted muon isolation SF는 적용 분기를 구분해 이중 적용하지 않는다.

### 물리 변경으로 관리할 항목

| 항목 | 확인된 차이 또는 주의점 |
|---|---|
| m4l | AN SR>200 GeV, 기존 histogram filling>50 GeV |
| Muon pT/eta 경계 | AN 표기는 >20, <2.4; 코드의 `<20`/`>2.4` reject는 경계 포함. Leading 52 GeV도 확인 |
| Neighbor dxy | AN은 abs(dxyBS), 기존 `checkIso`는 signed dxy>0.1 비교 |
| Muon F | AN의 두 ID 실패 설명과 기존 P-collection 제외 방식이 다름 |
| Direction/p4 | Pairing의 reco direction과 TuneP/corrected mass p4를 한 변수로 통합하면 결과가 바뀔 수 있음 |
| Fake deltaR | 선택된 F partner 기준과 global minimum same-flavor deltaR를 혼동하지 않음 |
| Negative bins | Raw prediction 보존과 최종 template clipping을 구분 |

차이를 발견했다고 즉시 수정하지 말고 기존 동작, 제안하는 변경, 선택 event/yield 변화를 기록한다.

## 11. 구현 전 결정할 사항

| 결정 | 권장 방향 | 아직 확인할 것 |
|---|---|---|
| 출력 architecture | Ntuple + 후단 C++/ROOT | 저장량/처리 시간과 관리 부담 |
| Production skim | SR/CR/variation 전체를 포함하는 느슨한 조건 | 정확한 HLT, pT, multiplicity 조건 |
| FF 측정 범위 | 적용과 측정을 별도 실행 모드로 구분 | 기존 payload 사용인지 새 측정인지 |
| 후단 backend | TTreeReader 또는 RDataFrame 중 하나 | Object/pair 조작과 유지보수 편의 |
| Config 형식 | 후단 독립 JSON 등 | 기존 환경과 의존성 |
| ID 수정 범위 | 기본은 기존 ID 재현, 필요한 입력 보존 | Threshold만 변경할지 변수 계산법도 바꿀지 |
| Systematics | 기존 재현 먼저, migration 비교 추가 | 최종 prescription |
| Schema | Event당 한 entry + object arrays | Branch별 필수/선택과 정밀도 |

설계 검토 LLM은 각 항목에 `채택 권장 / 수정 권장 / 확인 필요`를 붙이고 이유를 구체적으로 설명하라. 다음 단계가 막히는 결정만 사용자에게 질문하고, 구현 세부사항 전체를 질문 목록으로 넘기지 말라.

## 12. 구현이 승인될 경우의 순서와 완료 조건

1. 기준 snapshot 및 실제 사용할 CMSSW 환경 확인. 기존 resolved module의 의존성 정리.
2. Schema와 production skim 확정. 각 후단 cut이 어떤 branch를 쓰는지 mapping 작성.
3. Ntuple producer의 작은 event sample 실행. P/F와 trigger denominator event 보존 확인.
4. 후단 nominal 4P selection 구현. 세 채널의 event IDs, object indices, pairing, masses를 기존 코드와 대조.
5. 3P1F/2P2F와 isolated/nonisolated FF 적용 추가. Trigger anchor와 pair-ID 차이 확인.
6. Systematics, normalization, cutflow, background combination 추가.
7. 대표 sample으로 저장량과 처리 시간 측정 후 전체 생산 구조 결정.

검증 항목:

- 같은 입력에서 event ID 및 선택 pair 비교; float 비교 tolerance는 명시적으로 정한다.
- 1, 20, 28, 35, 50, 52, 200 GeV와 deltaR=0.3의 경계.
- Reciprocal electron pairing, tracker-only muon pair 실패, extra P/F object.
- 4P/3P1F/2P2F의 누락·중복 및 flavor 분류.
- Production skim으로 분석 대상 event와 variation migration이 손실되지 않음.
- Mll cfg 값만 바꾸어 후단 결과가 변하고 MiniAOD 재처리가 필요하지 않음.
- 병합 전후 normalization, Nraw/sumw/sumw2 일치.
- 없는 branch/payload 또는 schema 불일치는 원인을 식별할 수 있게 명시적으로 실패.

LLM의 검토 응답에는 추천 architecture, 유지/변경할 기존 코드, 확정할 결정, 첫 구현 범위, 실제 검증한 것과 검증하지 못한 것을 포함하라.

## 13. 기준 코드 링크

- [ResolvedEleCRanalyzer.cc](https://github.com/eoyun/ZprimeTo4l/blob/8c7d9c2ac352c376aad4dd166c9c6d8bcaff5d19/Analysis/plugins/ResolvedEleCRanalyzer.cc)
- [ResolvedEMuCRanalyzer.cc](https://github.com/eoyun/ZprimeTo4l/blob/8c7d9c2ac352c376aad4dd166c9c6d8bcaff5d19/Analysis/plugins/ResolvedEMuCRanalyzer.cc)
- [ResolvedMuCRanalyzer.cc](https://github.com/eoyun/ZprimeTo4l/blob/8c7d9c2ac352c376aad4dd166c9c6d8bcaff5d19/Analysis/plugins/ResolvedMuCRanalyzer.cc)
- [PairingHelper.cc](https://github.com/eoyun/ZprimeTo4l/blob/8c7d9c2ac352c376aad4dd166c9c6d8bcaff5d19/Analysis/src/PairingHelper.cc)
- [MuonCorrectionHelper.cc](https://github.com/eoyun/ZprimeTo4l/blob/8c7d9c2ac352c376aad4dd166c9c6d8bcaff5d19/Analysis/src/MuonCorrectionHelper.cc)
- [ElectronSystematicsHelper.cc](https://github.com/eoyun/ZprimeTo4l/blob/8c7d9c2ac352c376aad4dd166c9c6d8bcaff5d19/Analysis/src/ElectronSystematicsHelper.cc)
- [runMuMC_cfg.py](https://github.com/eoyun/ZprimeTo4l/blob/8c7d9c2ac352c376aad4dd166c9c6d8bcaff5d19/Analysis/test/runMuMC_cfg.py)
- [runRMFF.cc](https://github.com/eoyun/ZprimeTo4l/blob/8c7d9c2ac352c376aad4dd166c9c6d8bcaff5d19/Analysis/test/runRMFF.cc)

AN 원문이 필요하면 사용자에게 해당 PDF를 함께 전달받아 확인한다. 문서에 인용된 외부 calibration/ID 자료를 직접 읽지 않았다면 검증했다고 주장하지 않는다.
