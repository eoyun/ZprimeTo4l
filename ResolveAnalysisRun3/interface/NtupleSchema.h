#ifndef ResolveAnalysisRun3_NtupleSchema_h
#define ResolveAnalysisRun3_NtupleSchema_h
// 무엇: Events.root 의 branch 이름을 한 곳에서만 정의하는 공유 계약.
// 어떻게: writer(ResolvedNtuplizer)와 reader(NtupleReader)가 모두 이 이름 상수를 써서
//   tree->Branch(name, ...) / TTreeReaderArray(reader, name) 하므로 이름이 어긋날 수 없다.
// 의존: 없음(문자열 상수). 각 이름은 후단 Muon/Electron/Event struct 필드에 1:1 대응.
//
// 규칙(확장성):
//   - 모든 object 컬렉션은 "{prefix}_{field}" 가변길이 배열. ROOT std::vector branch 는
//     크기를 스스로 담으므로 별도 n_ branch 는 두지 않는다.
//   - 새 컬렉션(jet 등)을 추가하려면 새 namespace(jet)와 이름 상수만 더하면 된다.
//   - schema 를 바꾸면 kVersion 을 올린다.

namespace raRun3 {
namespace schema {

// schema 버전. branch 구성을 바꾸면 증가시킨다.
constexpr int kVersion = 1;

// ---- event 스칼라 ----
namespace ev {
constexpr const char* run       = "run";
constexpr const char* lumi      = "lumi";
constexpr const char* event     = "event";
constexpr const char* genWeight = "genWeight";  // MC weight (data=1)
constexpr const char* puTrue    = "puTrue";     // in-time true interactions (MC)
}  // namespace ev

// ---- muon 가변길이 배열 (후단 Muon struct 대응) ----
namespace mu {
constexpr const char* index           = "mu_index";
constexpr const char* charge          = "mu_charge";
constexpr const char* corrTunePpt      = "mu_corrTunePpt";
constexpr const char* rawTunePpt       = "mu_rawTunePpt";
constexpr const char* tunePeta         = "mu_tunePeta";
constexpr const char* tunePphi         = "mu_tunePphi";
constexpr const char* isHighPt         = "mu_isHighPt";
constexpr const char* isTrackerHighPt  = "mu_isTrackerHighPt";
constexpr const char* isTrackerMuon    = "mu_isTrackerMuon";
constexpr const char* trackIso         = "mu_trackIso";
constexpr const char* innerPt          = "mu_innerPt";
constexpr const char* recoEta          = "mu_recoEta";
constexpr const char* recoPhi          = "mu_recoPhi";
constexpr const char* recoVz           = "mu_recoVz";
constexpr const char* innerEta         = "mu_innerEta";
constexpr const char* innerPhi         = "mu_innerPhi";
constexpr const char* innerVz          = "mu_innerVz";
constexpr const char* innerDxyBS       = "mu_innerDxyBS";
constexpr const char* relPtErr         = "mu_relPtErr";
constexpr const char* trkLayers        = "mu_trkLayers";
constexpr const char* pixelHits        = "mu_pixelHits";
constexpr const char* matchedStations  = "mu_matchedStations";
constexpr const char* dxy              = "mu_dxy";
constexpr const char* dz               = "mu_dz";
}  // namespace mu

// ---- electron 가변길이 배열 (후단 Electron struct 대응) ----
namespace ele {
constexpr const char* index         = "ele_index";
constexpr const char* charge        = "ele_charge";
constexpr const char* selEt         = "ele_selEt";       // et() (ECAL-driven, PF 아님)
constexpr const char* etaSC         = "ele_etaSC";
constexpr const char* rawPt         = "ele_rawPt";       // polarP4()
constexpr const char* rawEta        = "ele_rawEta";
constexpr const char* rawPhi        = "ele_rawPhi";
constexpr const char* rawEnergy     = "ele_rawEnergy";   // energy()
constexpr const char* corrPt        = "ele_corrPt";      // ecalTrkEnergyPostCorr 적용
constexpr const char* corrEta       = "ele_corrEta";
constexpr const char* corrPhi       = "ele_corrPhi";
constexpr const char* corrM         = "ele_corrM";
constexpr const char* passModHeep   = "ele_passModHeep"; // electronID(...)
constexpr const char* modHeepBitmap = "ele_modHeepBitmap"; // userInt(...)
constexpr const char* addGsfIdx     = "ele_addGsfIdx";   // reciprocal-GSF 상대 index(-1=없음)
}  // namespace ele

// ---- trigger object 가변길이 배열 ----
namespace trig {
constexpr const char* pt         = "trigObj_pt";
constexpr const char* eta        = "trigObj_eta";
constexpr const char* phi        = "trigObj_phi";
constexpr const char* filterBits = "trigObj_filterBits";  // path/filter 연결 비트
}  // namespace trig

}  // namespace schema
}  // namespace raRun3
#endif
