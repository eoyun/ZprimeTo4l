#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Kinematics.h"
#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"

using raRun3::P4;
using raRun3::invariantMass;

// 정지-질량이 같은 두 입자를 정반대 방향으로 두면 불변질량은 해석적으로 계산된다.
static void test_back_to_back_dimuon() {
  const double mMu = 0.1056583745;
  // pt=45, eta=0, phi=0 과 phi=pi (정반대). 두 입자 E = sqrt(45^2 + mMu^2).
  P4 a{45.0, 0.0, 0.0, mMu};
  P4 b{45.0, 0.0, M_PI, mMu};
  // 합: px=0, py=0, pz=0, E=2*sqrt(45^2+mMu^2) → M = 2*sqrt(45^2+mMu^2)
  const double eOne = std::sqrt(45.0 * 45.0 + mMu * mMu);
  CHECK_CLOSE(invariantMass(a, b), 2.0 * eOne, 1e-6);
}

// 일반 pair 하나를 손 계산값과 비교 (회귀 방지용 고정값).
// 손 계산: px/py/pz/E 로 전개하면 M ≈ 79.094 GeV (아래 tol 0.05 안).
static void test_known_pair() {
  P4 a{40.0, 0.5, 0.0, 0.1056583745};
  P4 b{35.0, -0.4, 2.5, 0.1056583745};
  CHECK_CLOSE(invariantMass(a, b), 79.094, 0.05);
}

int main() {
  RUN(test_back_to_back_dimuon);
  RUN(test_known_pair);
  REPORT();
}
