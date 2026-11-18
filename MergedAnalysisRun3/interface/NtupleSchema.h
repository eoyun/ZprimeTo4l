#ifndef ZprimeTo4l_MergedAnalysisRun3_NtupleSchema_h
#define ZprimeTo4l_MergedAnalysisRun3_NtupleSchema_h

namespace mergedRun3::schema {
inline constexpr int version = 1;
inline constexpr float missing = -999.f;
// All arrays in a block have identical length; index is the MiniAOD source index.
// pt/eta/phi/mass are the unmodified pat object p4. TuneP is stored separately.
inline constexpr const char* muonFloats[] = {
    "pt", "eta", "phi", "mass", "energy", "vz",
    "tuneP_pt", "tuneP_eta", "tuneP_phi", "tuneP_relPtError",
    "inner_pt", "inner_eta", "inner_phi", "inner_vz", "inner_dxyBS", "dxy", "dz",
    "outer_pt", "outer_eta", "outer_phi", "global_pt", "global_eta", "global_phi",
    "best_pt", "best_eta", "best_phi", "normChi2_global", "normChi2_inner", "validFraction",
    "trackiso", "modifiedtrackiso", "modifiedIsoRun2", "relIsoRun2", "pfrelIso04",
    "em", "emS9", "emS25", "emMax", "had", "hadS9", "hadMax", "ho", "hoS9",
    "segCompatibility", "caloCompatibility", "chi2LocalPosition", "trkKink", "glbKink"};
inline constexpr const char* muonInts[] = {
    "index", "charge", "isHighPt", "isTrackerHighPt", "passMergedBaseline", "passResolved",
    "passResolvedAcceptance", "isGlobal", "isTracker", "isPF", "isLoose", "isMedium", "isTight",
    "trackIsoValid", "modifiedIsoValid", "pfIsoValid", "energyValid", "qualityValid", "matchesValid",
    "tuneP_valid", "inner_valid", "outer_valid", "global_valid", "best_valid",
    "isoNeighborStudy", "isoNeighborRun2", "numOfMatchedStations", "numOfChambers", "nMatches",
    "nValidMuonHits", "nValidPixelHits", "trkLayers", "pixelLayers"};
inline constexpr const char* electronFloats[] = {
    "pt", "eta", "phi", "mass", "energy", "et", "etaSC", "r9", "hOverE", "sigmaIetaIeta",
    "dEtaInSeed", "dPhiIn", "e2x5", "e5x5", "trackIso", "ecalIso", "hcalIso", "dxy", "dz"};
inline constexpr const char* electronInts[] = {
    "index", "charge", "gsfTrackValid", "ecalDriven", "missingHits", "passRun2Pt",
    "modifiedHeepRun3Valid", "passModifiedHeepRun3"};
}  // namespace mergedRun3::schema
#endif
