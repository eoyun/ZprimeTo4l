#ifndef ResolveAnalysisRun3_Electron_h
#define ResolveAnalysisRun3_Electron_h
// 무엇: 후단 electron 최소 표현. Plan 2에서 modified-HEEP 필드로 확장한다.
// 어떻게: 순수 데이터. kinematics 는 ECAL-driven+GSF (Particle Flow 아님).
// 의존: 없음.
namespace raRun3 {

struct Electron {
  int    index  = -1;  // ntuple 원본 collection index (P/F 교차 제외, pair index 용)
  int    charge = 0;
  double selEt  = 0.;   // [GeV] selection용 ET = pat::Electron::et()
  double etaSC  = 0.;   // supercluster eta (acceptance 기준)
  // 보정 전/후 mass kinematics
  double rawPt = 0., rawEta = 0., rawPhi = 0., rawEnergy = 0.;   // polarP4(), energy()
  double corrPt = 0., corrEta = 0., corrPhi = 0., corrM = 0.;    // ecalTrkEnergyPostCorr 적용
  bool   passModHeep = false;  // 기존 VID 결과 electronID("modifiedHeepElectronID")
  int    modHeepBitmap = 0;    // userInt("modifiedHeepElectronID") — F 마스크 판정용
  int    addGsfIdx     = -1;   // reciprocal-GSF 상대 electron index (없으면 -1)
};

}  // namespace raRun3
#endif
