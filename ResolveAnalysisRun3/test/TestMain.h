#ifndef ResolveAnalysisRun3_TestMain_h
#define ResolveAnalysisRun3_TestMain_h
// 의존성 없는 초경량 테스트 하네스. 프레임워크 대신 전부 읽어 이해할 수 있게 자작.
// CHECK(cond)          : 조건이 거짓이면 파일:라인과 함께 실패 기록
// CHECK_CLOSE(a,b,tol) : |a-b|>tol 이면 실패 기록 (부동소수 비교)
// RUN(fn)              : 테스트 함수 실행 (이름 출력)
// REPORT()             : 실패 개수에 따라 종료코드 반환 (0=성공, 1=실패)
#include <cstdio>
#include <cmath>

static int g_failures = 0;

#define CHECK(cond)                                                            \
  do {                                                                        \
    if (!(cond)) {                                                            \
      std::printf("  [FAIL] %s:%d  CHECK(%s)\n", __FILE__, __LINE__, #cond);  \
      ++g_failures;                                                           \
    }                                                                        \
  } while (0)

#define CHECK_CLOSE(a, b, tol)                                                 \
  do {                                                                        \
    const double _d = std::fabs((double)(a) - (double)(b));                   \
    if (_d > (double)(tol)) {                                                 \
      std::printf("  [FAIL] %s:%d  |%g - %g| = %g > %g\n", __FILE__,          \
                  __LINE__, (double)(a), (double)(b), _d, (double)(tol));     \
      ++g_failures;                                                           \
    }                                                                        \
  } while (0)

#define RUN(fn)                                                                \
  do {                                                                        \
    std::printf("[RUN ] %s\n", #fn);                                          \
    fn();                                                                     \
  } while (0)

#define REPORT()                                                              \
  do {                                                                        \
    if (g_failures) {                                                         \
      std::printf("FAILED: %d check(s)\n", g_failures);                       \
      return 1;                                                               \
    }                                                                        \
    std::printf("ALL PASS\n");                                               \
    return 0;                                                                 \
  } while (0)

#endif
