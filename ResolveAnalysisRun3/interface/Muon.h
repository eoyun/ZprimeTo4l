#ifndef ResolveAnalysisRun3_Muon_h
#define ResolveAnalysisRun3_Muon_h
// 무엇: 후단에서 다루는 muon 하나의 표현. ntuple branch 값이 그대로 들어온다.
// 어떻게: 순수 데이터(로직 없음). 단위/의미는 각 필드 주석 참고.
// 의존: 없음. (kinematics 는 P4 로 변환해 Kinematics.h 사용)
namespace raRun3 {

struct Muon {
  int    charge = 0;

  // --- kinematics (Particle Flow 아님: TuneP 사용) ---
  double corrTunePpt = 0.;  // [GeV] 보정 후 TuneP pt (selection/mass 기준)
  double rawTunePpt  = 0.;  // [GeV] 보정 전 TuneP pt (보정 검증용)
  double tunePeta    = 0.;  // TuneP 방향 eta (mass p4)
  double tunePphi    = 0.;  // TuneP 방향 phi (mass p4)

  // --- ID 결과 (CMSSW muon::isHighPtMuon / isTrackerHighPtMuon) ---
  bool   isHighPt        = false;  // global high-pT ID
  bool   isTrackerHighPt = false;  // tracker high-pT ID
  bool   isTrackerMuon   = false;  // muon::isTrackerMuon (F loose denominator 용)

  // --- modified isolation 재료 ---
  double trackIso = 0.;   // [GeV] raw track isolation
  double innerPt  = 0.;   // [GeV] inner-track pt (neighbor subtraction에 사용)

  // --- neighbor 기하 (checkIso 재현용) ---
  // self 로 쓸 때: recoEta/recoPhi/recoVz. neighbor 로 쓸 때: innerEta/innerPhi/innerVz/innerDxyBS.
  double recoEta    = 0.;  // muon reco 방향 eta (deltaR 기준: self)
  double recoPhi    = 0.;  // muon reco 방향 phi
  double recoVz     = 0.;  // [cm] muon vertex z (dz 비교: self)
  double innerEta   = 0.;  // inner-track eta (deltaR 기준: neighbor)
  double innerPhi   = 0.;  // inner-track phi
  double innerVz    = 0.;  // [cm] inner-track vertex z (dz 비교: neighbor)
  double innerDxyBS = 0.;  // [cm] inner-track signed dxy wrt beamspot (neighbor)

  // --- ID 입력변수 (nominal Plan1 에선 미사용, 향후 threshold 스캔용) ---
  double relPtErr        = 0.;
  int    trkLayers       = 0;
  int    pixelHits       = 0;
  int    matchedStations = 0;
  double dxy             = 0.;  // [cm] wrt PV
  double dz              = 0.;  // [cm] wrt PV
};

}  // namespace raRun3
#endif
