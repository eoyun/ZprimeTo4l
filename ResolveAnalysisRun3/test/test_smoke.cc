#include "ZprimeTo4l/ResolveAnalysisRun3/test/TestMain.h"

static void test_harness_basic() {
  CHECK(1 + 1 == 2);   // 통과해야 함
}

int main() {
  RUN(test_harness_basic);
  REPORT();
}
