#ifndef ResolveAnalysisRun3_Kinematics_h
#define ResolveAnalysisRun3_Kinematics_h
// 무엇: 4-운동량(P4)과 두 입자 불변질량 계산.
// 어떻게: (pt,eta,phi,mass) → (px,py,pz,E) 로 바꿔 더한 뒤 M=sqrt(E^2-|p|^2).
// 의존: 표준 <cmath> 만. CMSSW/ROOT 타입에 의존하지 않아 어디서든 테스트 가능.
#include <cmath>

namespace raRun3 {

// 하나의 lepton kinematics 표현. mass 는 입자 정지질량(예: muon 0.1057).
struct P4 {
  double pt;    // [GeV] 횡운동량 (muon: TuneP, electron: ecalTrk 보정값)
  double eta;   // pseudorapidity
  double phi;   // [rad] 방위각
  double mass;  // [GeV] 정지질량
};

// P4 → 데카르트 성분 (내부 헬퍼)
inline double px(const P4& v) { return v.pt * std::cos(v.phi); }
inline double py(const P4& v) { return v.pt * std::sin(v.phi); }
inline double pz(const P4& v) { return v.pt * std::sinh(v.eta); }
inline double energy(const P4& v) {
  return std::sqrt(px(v) * px(v) + py(v) * py(v) + pz(v) * pz(v) + v.mass * v.mass);
}

// phi 차이를 [-pi, pi] 로 감싼다.
inline double deltaPhi(double phi1, double phi2) {
  double d = phi1 - phi2;
  while (d >  M_PI) d -= 2.0 * M_PI;
  while (d < -M_PI) d += 2.0 * M_PI;
  return d;
}

// deltaR^2 = deta^2 + dphi^2 (pairing/neighbor 공용).
inline double deltaR2(double eta1, double phi1, double eta2, double phi2) {
  const double deta = eta1 - eta2;
  const double dphi = deltaPhi(phi1, phi2);
  return deta * deta + dphi * dphi;
}

// 두 입자의 불변질량 M = sqrt( (E1+E2)^2 - |p1+p2|^2 ).
inline double invariantMass(const P4& a, const P4& b) {
  const double e = energy(a) + energy(b);
  const double sx = px(a) + px(b);
  const double sy = py(a) + py(b);
  const double sz = pz(a) + pz(b);
  const double m2 = e * e - (sx * sx + sy * sy + sz * sz);
  return m2 > 0.0 ? std::sqrt(m2) : 0.0;  // 수치오차로 음수면 0
}

}  // namespace raRun3
#endif
