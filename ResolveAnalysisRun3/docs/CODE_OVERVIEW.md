# ResolveAnalysisRun3 — 코드 개요 (파일별 용도 + 스켈레톤)

각 파일의 **용도**와 **큰 스켈레톤**을 주석으로 정리한 지도. 세부 구현은 실제 파일 참고.
전체 구조: `[MiniAOD] → (plugins) ResolvedNtuplizer → Events.root → (src+interface) 후단 → hist`.

디렉토리 역할:
- `interface/` : 헤더(선언·struct·공유 계약). 후단 로직의 "무엇".
- `src/`       : `interface/`의 구현. "어떻게".
- `plugins/`   : cmsRun EDAnalyzer(ntuplizer). MiniAOD 의존.
- `python/`    : ntuplizer cfi(기본 PSet).
- `test/`      : 합성 struct 단위 테스트(CMSSW/데이터 불필요).

---

## interface/NtupleSchema.h  ★공유 계약
**용도:** Events.root branch **이름을 한 곳에서만** 정의. writer(ntuplizer)와 reader가 같은 상수를
써서 이름이 어긋날 수 없게 함. schema 바꾸면 `kVersion` 증가.
```cpp
namespace raRun3::schema {
  constexpr int kVersion = 1;
  namespace ev  { /* run, lumi, event, genWeight, puTrue, hltFired, passMETfilters */ }
  namespace mu  { /* mu_index, mu_corrTunePpt, mu_rawTunePpt, ... mu_dz  (24개) */ }
  namespace ele { /* ele_index, ele_selEt, ele_rawPt.., ele_corrPt.., addGsfIdx (15개) */ }
  namespace trig{ /* trigObj_pt/eta/phi/filterBits */ }
}
// 확장: jet 추가 시 namespace jet 만 더하면 됨.
```

## interface/Kinematics.h
**용도:** 4-운동량과 각도 계산의 최소 단위. CMSSW/ROOT 타입 의존 없음 → 어디서든 테스트.
```cpp
namespace raRun3 {
  struct P4 { double pt, eta, phi, mass; };          // lepton kinematics 1개
  double px/py/pz/energy(const P4&);                 // 데카르트 변환(내부)
  double deltaPhi(phi1,phi2);                         // [-pi,pi] 감쌈
  double deltaR2(eta1,phi1,eta2,phi2);               // pairing/neighbor 공용
  double invariantMass(const P4& a, const P4& b);    // M = sqrt((E1+E2)^2-|p1+p2|^2)
}
```

## interface/Muon.h
**용도:** 후단 muon 1개의 표현(= ntuple에서 읽어온 값). 순수 데이터, 로직 없음.
```cpp
struct Muon {
  int index, charge;                                 // index=원본 collection 위치(P/F·pair용)
  double corrTunePpt, rawTunePpt, tunePeta, tunePphi;// kinematics: TuneP(=PF 아님)
  bool isHighPt, isTrackerHighPt, isTrackerMuon;     // ID 결과(재현 기준)
  double trackIso, innerPt;                          // modified iso 재료
  double recoEta/Phi/Vz, innerEta/Phi/Vz, innerDxyBS;// neighbor 기하(checkIso 재현)
  double relPtErr; int trkLayers, pixelHits, matchedStations; double dxy, dz; // ID 입력/F용
};
```

## interface/Electron.h
**용도:** 후단 electron 1개. kinematics는 ECAL-driven+GSF(PF 아님).
```cpp
struct Electron {
  int index, charge;
  double selEt;                                      // 선택용 ET = et()
  double etaSC;                                      // supercluster eta(acceptance)
  double rawPt/Eta/Phi/Energy;                       // 보정 전(polarP4, energy)
  double corrPt/Eta/Phi/M;                           // 보정 후(ecalTrkEnergyPostCorr)
  bool passModHeep; int modHeepBitmap;               // VID 결과 + bitmap(F 마스크용)
  int addGsfIdx;                                     // reciprocal-GSF 상대 index(-1=없음)
};
```

## interface/Event.h
**용도:** 후단 event 1개. 스칼라 + muon/electron 가변길이 배열. NtupleReader(Plan 4)가 채움.
```cpp
struct Event {
  unsigned run, lumi; uint64_t event; double genWeight;
  std::vector<Muon> muons;
  std::vector<Electron> electrons;
};
```

## interface/Config.h  +  src/Config.cc
**용도:** 후단 모든 물리 컷을 담는 struct + JSON 로더. **누락 키는 예외**(기본값으로 안 때움).
```cpp
// Config.h
struct NeighborCuts { drMin, drMax, dzMax, dxyBSMax; };
struct MuonFakeCuts { trkLayersMin, pixelHitsMin, matchedStationsMin, dxyMax, dzMax; };
struct MuonCuts    { tunePptMin, etaMax, modIsoRelMax; NeighborCuts neighbor; MuonFakeCuts fake; };
struct ElectronCuts{ gapLo, gapHi, eeEtaMax, heepMaskLoose, heepAllPass; };
struct MassCuts    { dileptonMassMin, signalMassMin; };
struct Config      { MuonCuts muon; ElectronCuts electron; MassCuts massCuts; };
Config loadConfig(path);  Config loadConfigFromString(json);
// Config.cc: required<T>(node,key,where) 로 키 강제 → 없으면 runtime_error("missing ...")
```

## interface/ObjectSelector.h  +  src/ObjectSelector.cc
**용도:** object P/F 판정. **순수 함수**(히스토그램/파일 접근 없음). 기존 analyzer 재현.
```cpp
namespace ObjectSelector {
  // --- electron ---
  bool passElectronAccept(e,cfg);                    // |etaSC|<eeEtaMax & gap 밖
  vector<Electron> selectElectronsP(all,cfg);        // accept + passModHeep
  vector<Electron> selectElectronsF(all,cfg);        // P아님 + accept + (bitmap|maskLoose)==allPass
  // --- muon ---
  bool passMuonAccept(mu,cfg);                        // corrTunePpt>=min & |recoEta|<=max (경계포함)
  bool passMuonId(mu);                               // isHighPt || isTrackerHighPt
  bool isMuonNeighbor(self,other,cfg);               // checkIso 재현(0.01<dR<0.3,|dz|<0.2,dxyBS<=0.1)
  double modifiedMuonIso(self,idPool,cfg);           // raw iso - (최고 TuneP neighbor innerPt)
  vector<Muon> selectMuonsP(all,cfg);                // accept+ID pool → iso/pt<컷
  vector<Muon> selectMuonsF(all,passP,cfg);          // tracker+accept+P제외+loose track 조건
}
```

## interface/CandidateBuilder.h  +  src/CandidateBuilder.cc
**용도:** 채널별 pairing(4e/4mu/2e2mu). 기존 PairingHelper 재현. PairIdx=입력 vector 위치.
```cpp
namespace CandidateBuilder {
  struct PairIdx { int first, second; };
  struct Candidate { PairIdx pair1, pair2; bool valid; };
  Candidate buildFourMuon(mus);        // closest-dR (3조합 중 min(dR1,dR2) 최소)
  Candidate buildFourElectron(eles);   // reciprocal-GSF 우선 → 실패시 dR
  Candidate buildTwoETwoMu(eles,mus);  // ee 1쌍 + μμ 1쌍
  // (내부) chooseByDR(eta[4],phi[4],p1,p2): 3조합 dR 비교
}
```

## plugins/ResolvedNtuplizer.cc  (cmsRun EDAnalyzer)
**용도:** MiniAOD → Events.root. acceptance 수준(ID 무관) object 저장, correction 전/후,
느슨한 skim, skim 이전 정규화. 물리 최종컷은 안 함.
```cpp
class ResolvedNtuplizer : one::EDAnalyzer<one::SharedResources> {
  beginJob():  TTree "Events" + 모든 branch(NtupleSchema 이름) 등록;
               norm TH1(Nevents/sumw/sumw2), meta TH1(kVersion, muonCorrApplied=0)
  analyze():   clear
               weight(genWeight,puTrue) → norm 에 skim 이전 합산
               run/lumi/event
               HLT(trigList OR) → hltFired ; METfilter → passMETfilters
               trigger objects(원하는 path 매칭) 저장
               PV 없으면 return
               muon 루프:  tunePpt>=20 & |eta|<2.4 (ID무관) → 모든 mu_* branch
                          (corrTunePpt=raw, meta muonCorrApplied=0 로 미적용 표시)
               electron 루프: |etaSC|<2.5 & gap밖 (ID무관) → 모든 ele_* branch
                          (corr = ecalTrkEnergyPostCorr 진짜값; addGsfIdx 계산)
               skim: (nMu+nEle) < minLeptons 면 return
               tree_->Fill()
}
DEFINE_FWK_MODULE(ResolvedNtuplizer);
```

## python/ResolvedNtuplizer_cfi.py
**용도:** ntuplizer 기본 PSet(InputTag·skim 기준). cfg에서 override.
```python
resolvedNtuplizer = cms.EDAnalyzer('ResolvedNtuplizer',
  isMC, srcMuon, srcEle, srcPv, beamSpot,
  addGsfTrk=modifiedHEEPIDVarValueMaps2nd:eleAddGsfTrk,
  generator, pileupSummary, triggerResults, triggerObjects, METfilters,
  trigList=[], METfilterList=[],                 # 사용자 조사
  muPtMinStore, muEtaMaxStore, eleEtaMaxStore, eleGapLo, eleGapHi, minLeptonsSkim)
```

## test/runResolvedNtuplizer_cfg.py  (예시 cfg — 틀)
**용도:** ntuplizer 실행 설정. [USER] GT/source/era/trigList/METfilterList만 채우면 됨.
```python
process = cms.Process('raRun3Ntuple')
# [USER] source, GlobalTag, era
# producer 체인(필수): ModifiedHEEPIDVarValueMaps + ModifiedEcalRecHitIsolationScone
#                     + egammaPostRecoSeq(ecalTrkEnergyPostCorr, modifiedHeepElectronID)
#                     + modifiedHEEPIDVarValueMaps2nd(addGsfTrk)
# process.p = Path(producers... + resolvedNtuplizer)
# 출력: Events.root
```

## test/TestMain.h  +  test/test_*.cc
**용도:** 의존성 0의 자작 테스트 하네스 + 각 모듈 단위 테스트(합성 struct).
```cpp
// TestMain.h : CHECK(cond) / CHECK_CLOSE(a,b,tol) / RUN(fn) / REPORT()
// test_smoke.cc            : 하네스 동작
// test_kinematics.cc       : 불변질량(back-to-back, 회귀값 79.094)
// test_config.cc           : JSON 파싱 + 누락키 예외
// test_muon_selection.cc   : accept(recoEta,경계) / neighbor / modIso / P / F
// test_electron_selection.cc: accept / P / F(bitmap 마스크)
// test_candidate.cc        : 4mu dR / 4e reciprocal·fallback / 2e2mu
// test_schema.cc           : branch 이름 prefix·유일성(50개)
```

## BuildFile.xml (3개)
**용도:** scram 빌드 규칙.
```
ResolveAnalysisRun3/BuildFile.xml        : src/*.cc → libZprimeTo4lResolveAnalysisRun3 (root,json export)
ResolveAnalysisRun3/plugins/BuildFile.xml: ResolvedNtuplizer → edm plugin
ResolveAnalysisRun3/test/BuildFile.xml   : test_*.cc → 실행파일 (패키지 lib 링크)
```

---

## 데이터 흐름 요약
```
MiniAOD
  │  plugins/ResolvedNtuplizer  (NtupleSchema 이름으로 branch 채움)
  ▼
Events.root  ── mu_*/ele_*/trig_*/ev_* + norm/meta
  │  (Plan 4) NtupleReader: NtupleSchema 이름으로 읽어 Event/Muon/Electron struct 로
  ▼
ObjectSelector(P/F) → CandidateBuilder(pairing) → [Plan 4] RegionSelector(채널/trigger/mass/SR·CR)
  │                                                     ↑ Config(selection.json)
  ▼
[Plan 4] HistogramWriter → histograms.root
```
아직 없음(다음 Plan): NtupleReader, RegionSelector, WeightProvider, HistogramWriter, runResolved main.
