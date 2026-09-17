// 무엇: muon 선택 로직 구현. 기존 checkIso / ResolvedMuCRanalyzer 재현.
// 어떻게: neighbor 후보를 TuneP pt 로 정렬해 최고 하나만 빼는 modified iso.
// 의존: ObjectSelector.h, <algorithm>, <cmath>.
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/ObjectSelector.h"

#include <algorithm>
#include <cmath>

namespace {

// deltaR^2 (phi 는 [-pi,pi] 로 감싼다).
double deltaR2(double eta1, double phi1, double eta2, double phi2) {
  double dphi = phi1 - phi2;
  while (dphi >  M_PI) dphi -= 2.0 * M_PI;
  while (dphi < -M_PI) dphi += 2.0 * M_PI;
  const double deta = eta1 - eta2;
  return deta * deta + dphi * dphi;
}

}  // namespace

namespace raRun3 {
namespace ObjectSelector {

bool passMuonAccept(const Muon& mu, const Config& cfg) {
  return mu.corrTunePpt > cfg.muon.tunePptMin &&
         std::abs(mu.tunePeta) < cfg.muon.etaMax;
}

bool passMuonId(const Muon& mu) {
  return mu.isHighPt || mu.isTrackerHighPt;
}

bool isMuonNeighbor(const Muon& self, const Muon& other, const Config& cfg) {
  // self reco 방향 vs neighbor inner-track 방향 (기존 checkIso 와 동일)
  const double dr2 = deltaR2(self.recoEta, self.recoPhi, other.innerEta, other.innerPhi);
  const double drMin2 = cfg.muon.neighbor.drMin * cfg.muon.neighbor.drMin;
  const double drMax2 = cfg.muon.neighbor.drMax * cfg.muon.neighbor.drMax;
  if (dr2 > drMax2 || dr2 < drMin2) return false;
  if (std::abs(self.recoVz - other.innerVz) > cfg.muon.neighbor.dzMax) return false;
  // 기존 checkIso 는 signed dxy 를 그대로 비교(> 면 탈락). 재현.
  if (other.innerDxyBS > cfg.muon.neighbor.dxyBSMax) return false;
  return true;
}

double modifiedMuonIso(const Muon& self, const std::vector<Muon>& idPool, const Config& cfg) {
  // self 와 다른 muon 만 neighbor 후보. self 를 idPool 에 포함해 넘겨도
  // isMuonNeighbor 의 drMin(dR≈0 배제) 이 self 자신을 걸러낸다.
  double iso = self.trackIso;

  // neighbor 후보 수집
  std::vector<const Muon*> cands;
  for (const auto& other : idPool) {
    if (isMuonNeighbor(self, other, cfg))
      cands.push_back(&other);
  }
  if (cands.empty()) return iso;

  // TuneP pt 내림차순 정렬 후 최고 하나의 innerPt 를 뺀다.
  std::sort(cands.begin(), cands.end(),
            [](const Muon* a, const Muon* b) { return a->corrTunePpt > b->corrTunePpt; });
  iso -= cands.front()->innerPt;
  return iso;
}

std::vector<Muon> selectMuonsP(const std::vector<Muon>& all, const Config& cfg) {
  // 1) accept + ID 통과 pool
  std::vector<Muon> idPool;
  for (const auto& mu : all) {
    if (passMuonAccept(mu, cfg) && passMuonId(mu))
      idPool.push_back(mu);
  }
  // 2) 각 pool muon 의 modified iso 로 격리 판정
  std::vector<Muon> passP;
  for (const auto& mu : idPool) {
    const double iso = modifiedMuonIso(mu, idPool, cfg);
    if (iso / mu.corrTunePpt < cfg.muon.modIsoRelMax)
      passP.push_back(mu);
  }
  return passP;
}

}  // namespace ObjectSelector
}  // namespace raRun3
