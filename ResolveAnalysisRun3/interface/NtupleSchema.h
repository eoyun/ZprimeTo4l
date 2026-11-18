#ifndef ResolveAnalysisRun3_NtupleSchema_h
#define ResolveAnalysisRun3_NtupleSchema_h
// 무엇: Events.root 의 flat branch 이름을 한 곳에서만 정의하는 공유 계약.
// 어떻게: 모든 branch 는 std::vector<float/int> "flat" 배열 (muon_pt, electron_pt ...).
//   외부 코드(condor/uproot/RDataFrame)가 우리 C++ struct 없이 바로 읽을 수 있다.
//   writer(ResolvedNtuplizer)와 reader 가 같은 이름 상수를 써서 어긋날 수 없다.
// 규칙: object 컬렉션은 "{object}_{var}" flat 배열. ROOT vector branch 는 크기를 스스로 담음.
//   새 object(jet)는 namespace 추가만. schema 변경 시 kVersion 증가.
//
// 코드 가독성용으로 namespace(mu/ele/...)로 묶지만 branch 문자열은 전부 flat.

namespace raRun3 {
namespace schema {

constexpr int kVersion = 7;  // v3: electron scale/smearing 후단(correctionlib) 적용 →
                             //     ntuplizer 는 raw + 보정입력(r9,seedGain)만 저장(corr 제거)
                             // v4: 표준 HEEP V70 pass/bitmap 추가 (modified HEEP 대조용)
                             // v5: hard-process final-state gen lepton(e/mu) 저장 (gen 매칭은 후처리)

// v6: ModifiedHEEPRun3ID pass and per-cut bitmap added.
// v7: modified ECAL iso values read from 2nd-instance producers keyed on the post-reco
//     slimmedElectrons (v6 and earlier: electron_modEcalHcalIso was always -999 because the
//     ValueMap was keyed on the input-MiniAOD electrons); electron_modEcalHcalIsoRun3 added.

// ---- event 스칼라 ----
namespace ev {
constexpr const char* run            = "run";
constexpr const char* lumi           = "lumi";
constexpr const char* event          = "event";
constexpr const char* genWeight      = "genWeight";      // MC weight (data=1)
constexpr const char* puTrue         = "puTrue";         // in-time true interactions (MC)
constexpr const char* nPV            = "nPV";            // primary vertex 개수
constexpr const char* rho            = "rho";            // fixedGridRhoFastjetAll (iso 보정용)
constexpr const char* hltFired       = "hltFired";       // wanted path OR (0/1)
constexpr const char* passMETfilters = "passMETfilters"; // 요구 filter 전부 통과 (0/1)
}  // namespace ev

// ---- muon flat 배열 (prefix muon_) ----
namespace mu {
constexpr const char* index           = "muon_index";        // 원본 collection index
constexpr const char* charge          = "muon_charge";
constexpr const char* pt              = "muon_pt";           // TuneP pt (보정 후; 현재 raw=corr)
constexpr const char* ptRaw           = "muon_ptRaw";        // TuneP pt (보정 전)
constexpr const char* eta             = "muon_eta";          // TuneP eta
constexpr const char* phi             = "muon_phi";          // TuneP phi
constexpr const char* recoPt          = "muon_recoPt";      // reco 방향(neighbor self)
constexpr const char* recoEta         = "muon_recoEta";      // reco 방향(neighbor self)
constexpr const char* recoPhi         = "muon_recoPhi";
constexpr const char* recoVz          = "muon_recoVz";
constexpr const char* innerPt         = "muon_innerPt";      // inner-track pt (iso subtraction)
constexpr const char* innerEta        = "muon_innerEta";
constexpr const char* innerPhi        = "muon_innerPhi";
constexpr const char* innerVz         = "muon_innerVz";
constexpr const char* innerDxyBS      = "muon_innerDxyBS";
constexpr const char* isHighPt        = "muon_isHighPt";
constexpr const char* isTrackerHighPt = "muon_isTrackerHighPt";
constexpr const char* isTrackerMuon   = "muon_isTrackerMuon";
constexpr const char* trackIso        = "muon_trackIso";     // raw track iso
constexpr const char* relPtErr        = "muon_relPtErr";
constexpr const char* trkLayers       = "muon_trkLayers";
constexpr const char* pixelHits       = "muon_pixelHits";
constexpr const char* matchedStations = "muon_matchedStations";
constexpr const char* dxy             = "muon_dxy";          // wrt PV
constexpr const char* dz              = "muon_dz";           // wrt PV
}  // namespace mu

// ---- electron flat 배열 (prefix electron_) ----
// kinematics 는 MiniAOD slimmedElectron 의 reco 값(= EGM scale/smearing "적용 전").
// 실제 scale(data)/smear(MC) 보정은 후단에서 correctionlib 로 적용한다(deferred).
// 그래서 여기엔 raw kinematics + 보정 입력(etaSC, r9, seedGain)만 저장하고, corrected 는 저장 안 함.
namespace ele {
constexpr const char* index       = "electron_index";
constexpr const char* charge      = "electron_charge";
constexpr const char* et          = "electron_et";          // selection ET = et() (ECAL, PF 아님)
constexpr const char* ptRaw       = "electron_ptRaw";       // polarP4() pt (보정 전 = correctionlib 입력 pt)
constexpr const char* etaRaw      = "electron_etaRaw";
constexpr const char* phiRaw      = "electron_phiRaw";
constexpr const char* energyRaw   = "electron_energyRaw";   // energy()
constexpr const char* etaSC       = "electron_etaSC";       // supercluster eta (acceptance + 보정 입력 ScEta)
constexpr const char* r9          = "electron_r9";          // full5x5_r9 (보정 입력 + shower shape)
constexpr const char* seedGain    = "electron_seedGain";    // seed 결정체 gain(12/6/1) (보정 입력)
// ID 결과
constexpr const char* passModHeep   = "electron_passModHeep";   // electronID(...)
constexpr const char* modHeepBitmap = "electron_modHeepBitmap"; // userInt(...) 12-bit
constexpr const char* passHeepV70   = "electron_passHeepV70";   // electronID("heepElectronID-HEEPV70")
constexpr const char* heepV70Bitmap = "electron_heepV70Bitmap"; // userInt("heepElectronID-HEEPV70") 12-bit
constexpr const char* passModifiedHEEPRun3ID = "electron_passModifiedHEEPRun3ID";
constexpr const char* bitmapModifiedHEEPRun3ID = "electron_bitmapModifiedHEEPRun3ID"; // 12-bit
constexpr const char* addGsfIdx     = "electron_addGsfIdx";     // reciprocal-GSF 상대 index(-1=없음)
// --- modified-HEEP study: 표준 HEEP 입력 (pat::Electron 메서드) ---
constexpr const char* hOverE        = "electron_hOverE";           // hadronicOverEm  [bit6]
constexpr const char* sigmaIeta     = "electron_sigmaIeta";        // full5x5_sigmaIetaIeta [bit4]
constexpr const char* dEtaInSeed    = "electron_dEtaInSeed";       // [bit2]
constexpr const char* dPhiIn        = "electron_dPhiIn";           // [bit3]
constexpr const char* e2x5          = "electron_e2x5";             // full5x5_e2x5Max [bit5]
constexpr const char* e5x5          = "electron_e5x5";             // full5x5_e5x5    [bit5]
constexpr const char* dr03TkSumPt   = "electron_dr03TkSumPt";      // 표준 track iso(비교용)
constexpr const char* dr03EcalIso   = "electron_dr03EcalRecHitSumEt";
constexpr const char* dr03HcalIso   = "electron_dr03HcalTowerSumEt";
constexpr const char* missingHits   = "electron_missingInnerHits"; // [bit10]
constexpr const char* dxy           = "electron_dxy";              // [bit9]
constexpr const char* dz            = "electron_dz";
constexpr const char* ecalDriven    = "electron_ecalDriven";       // [bit11]
// --- modified-HEEP study: modified 변수 (ModifiedHEEP ValueMap; 없으면 -999 sentinel) ---
constexpr const char* modTrkIso        = "electron_modTrkIso";        // eleTrkPtIso [bit7]
constexpr const char* modEcalHcalIso   = "electron_modEcalHcalIso";   // EcalRecHitIso [bit8]
constexpr const char* modEcalHcalIsoRun3 = "electron_modEcalHcalIsoRun3"; // Run3 EcalRecHitIso [Run3 bit8]
constexpr const char* union5x5covIeIe  = "electron_union5x5covIeIe";
constexpr const char* union5x5covIeIp  = "electron_union5x5covIeIp";
constexpr const char* union5x5covIpIp  = "electron_union5x5covIpIp";
constexpr const char* union5x5dEtaIn   = "electron_union5x5dEtaIn";
constexpr const char* union5x5dPhiIn   = "electron_union5x5dPhiIn";
constexpr const char* union5x5Energy   = "electron_union5x5Energy";
constexpr const char* dEtaInSeed2nd    = "electron_dEtaInSeed2nd";
constexpr const char* dPhiInSC2nd      = "electron_dPhiInSC2nd";
constexpr const char* dPerpIn          = "electron_dPerpIn";
constexpr const char* alphaTrack       = "electron_alphaTrack";
constexpr const char* alphaCalo        = "electron_alphaCalo";
constexpr const char* normDParaIn      = "electron_normDParaIn";
}  // namespace ele

// ---- trigger object flat 배열 ----
namespace trig {
constexpr const char* pt         = "trigObj_pt";
constexpr const char* eta        = "trigObj_eta";
constexpr const char* phi        = "trigObj_phi";
constexpr const char* filterBits = "trigObj_filterBits";
}  // namespace trig

// ---- hard-process final-state gen lepton flat 배열 (prefix gen_) ----
// 선택: fromHardProcessFinalState() && (|pdgId|==11 || |pdgId|==13). MC 만 채움.
// reco<->gen 매칭은 저장 안 함(후처리). ΔR_gen pairing 은 momIdx 로(같은 A 형제 = 같은 momIdx).
namespace gen {
constexpr const char* pt     = "gen_pt";
constexpr const char* eta    = "gen_eta";
constexpr const char* phi    = "gen_phi";
constexpr const char* pdgId  = "gen_pdgId";   // 11=e-, -11=e+, 13=mu-, -13=mu+ (부호=charge)
constexpr const char* momIdx = "gen_momIdx";  // mother 의 prunedGenParticles index (형제 그룹핑 키; 없으면 -1)
}  // namespace gen

// modified ValueMap 을 못 읽었을 때(잘못된 InputTag 등) 저장할 명시적 sentinel.
constexpr float kMissing = -999.f;

}  // namespace schema
}  // namespace raRun3
#endif
