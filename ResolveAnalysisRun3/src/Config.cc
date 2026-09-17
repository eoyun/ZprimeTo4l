// 무엇: JSON → Config 변환. 누락 키는 명확한 메시지와 함께 예외.
// 어떻게: nlohmann/json 으로 파싱하고 required<T>() 헬퍼로 키 존재를 강제.
// 의존: nlohmann/json (scram tool "json"), Config.h.
#include "ZprimeTo4l/ResolveAnalysisRun3/interface/Config.h"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <nlohmann/json.hpp>

using nlohmann::json;

namespace {

// 키가 없으면 어느 경로에서 무엇이 빠졌는지 알려주고 예외.
template <typename T>
T required(const json& node, const char* key, const std::string& where) {
  if (!node.contains(key))
    throw std::runtime_error("Config: missing key '" + std::string(key) +
                             "' under " + where);
  return node.at(key).get<T>();
}

raRun3::Config parse(const json& j) {
  raRun3::Config cfg;

  if (!j.contains("objects") || !j.at("objects").contains("muon"))
    throw std::runtime_error("Config: missing 'objects.muon'");
  const json& mu = j.at("objects").at("muon");
  cfg.muon.tunePptMin   = required<double>(mu, "tunePptMin", "objects.muon");
  cfg.muon.etaMax       = required<double>(mu, "etaMax", "objects.muon");
  cfg.muon.modIsoRelMax = required<double>(mu, "modIsoRelMax", "objects.muon");

  if (!mu.contains("neighbor"))
    throw std::runtime_error("Config: missing 'objects.muon.neighbor'");
  const json& nb = mu.at("neighbor");
  cfg.muon.neighbor.drMin    = required<double>(nb, "drMin", "objects.muon.neighbor");
  cfg.muon.neighbor.drMax    = required<double>(nb, "drMax", "objects.muon.neighbor");
  cfg.muon.neighbor.dzMax    = required<double>(nb, "dzMax", "objects.muon.neighbor");
  cfg.muon.neighbor.dxyBSMax = required<double>(nb, "dxyBSMax", "objects.muon.neighbor");

  if (!mu.contains("fake"))
    throw std::runtime_error("Config: missing 'objects.muon.fake'");
  const json& mf = mu.at("fake");
  cfg.muon.fake.trkLayersMin       = required<int>(mf, "trkLayersMin", "objects.muon.fake");
  cfg.muon.fake.pixelHitsMin       = required<int>(mf, "pixelHitsMin", "objects.muon.fake");
  cfg.muon.fake.matchedStationsMin = required<int>(mf, "matchedStationsMin", "objects.muon.fake");
  cfg.muon.fake.dxyMax             = required<double>(mf, "dxyMax", "objects.muon.fake");
  cfg.muon.fake.dzMax              = required<double>(mf, "dzMax", "objects.muon.fake");

  if (!j.at("objects").contains("electron"))
    throw std::runtime_error("Config: missing 'objects.electron'");
  const json& el = j.at("objects").at("electron");
  cfg.electron.gapLo        = required<double>(el, "gapLo", "objects.electron");
  cfg.electron.gapHi        = required<double>(el, "gapHi", "objects.electron");
  cfg.electron.eeEtaMax     = required<double>(el, "eeEtaMax", "objects.electron");
  cfg.electron.heepMaskLoose = required<int>(el, "heepMaskLoose", "objects.electron");
  cfg.electron.heepAllPass   = required<int>(el, "heepAllPass", "objects.electron");

  if (!j.contains("massCuts"))
    throw std::runtime_error("Config: missing 'massCuts'");
  const json& mc = j.at("massCuts");
  cfg.massCuts.dileptonMassMin = required<double>(mc, "dileptonMassMin", "massCuts");
  cfg.massCuts.signalMassMin   = required<double>(mc, "signalMassMin", "massCuts");

  return cfg;
}

}  // namespace

namespace raRun3 {

Config loadConfigFromString(const std::string& jsonText) {
  return parse(json::parse(jsonText));
}

Config loadConfig(const std::string& path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Config: cannot open file '" + path + "'");
  std::stringstream ss;
  ss << in.rdbuf();
  return loadConfigFromString(ss.str());
}

}  // namespace raRun3
