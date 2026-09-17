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

// muon F(loose denominator) 전용 컷 (기존 ResolvedMuCRanalyzer nonHighPtMuons 재현).
struct MuonFakeCuts {
  int    trkLayersMin       = 0;  // trkLayers > min (원본 >5 → min=5)
  int    pixelHitsMin       = 0;  // pixelHits > min (원본 >0 → min=0)
  int    matchedStationsMin = 0;  // matchedStations >= min (원본 >=1 → min=1)
  double dxyMax             = 0.; // |dxy(PV)| < max (0.2)
  double dzMax              = 0.; // |dz(PV)|  < max (0.5)
};

struct MuonCuts {
  double       tunePptMin = 0.;  // [GeV]
  double       etaMax     = 0.;  // |eta| 상한
  double       modIsoRelMax = 0.;  // modified iso / TuneP pt 상한
  NeighborCuts neighbor;
  MuonFakeCuts fake;
};

struct ElectronCuts {
  double gapLo    = 0.;   // EB-EE gap 하한 (|etaSC| 1.4442)
  double gapHi    = 0.;   // EB-EE gap 상한 (|etaSC| 1.566)
  double eeEtaMax = 0.;   // |etaSC| 상한 (2.5)
  int    heepMaskLoose = 0;  // 0x7B0 = 1968 : F denominator 마스크
  int    heepAllPass   = 0;  // 0xFFF = 4095 : 12 컷 전부 통과 값
};

struct MassCuts {
  double dileptonMassMin = 0.;  // [GeV] 각 pair mll 하한
  double signalMassMin   = 0.;  // [GeV] m4l SR 하한
};

struct Config {
  MuonCuts     muon;
  ElectronCuts electron;
  MassCuts     massCuts;
};

// 파일 경로 / 문자열에서 로드. 필수 키 누락 시 std::runtime_error.
Config loadConfig(const std::string& path);
Config loadConfigFromString(const std::string& jsonText);

}  // namespace raRun3
#endif
