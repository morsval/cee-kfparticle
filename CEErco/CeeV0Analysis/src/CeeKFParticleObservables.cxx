#include "CeeKFParticleObservables.h"
#include <cmath>
#include <cstdlib>
#include <stdexcept>
namespace CeeV0Analysis {
namespace {
double NaN() { return std::numeric_limits<double>::quiet_NaN(); }
double Ratio(double x,double y) { return std::isfinite(x) && std::isfinite(y) && y>0. ? x/y : NaN(); }
double SafeSqrt(double x) { return std::isfinite(x) && x>=0. ? std::sqrt(x) : NaN(); }
double Rapidity(double e,double pz) { return e>std::fabs(pz) ? 0.5*std::log((e+pz)/(e-pz)) : NaN(); }
double MotherDca(const CeeV0Lite::Candidate& c,const AnalysisEvent& e) {
  if (c.type!=CeeV0Lite::kLambda && c.type!=CeeV0Lite::kAntiLambda && c.type!=CeeV0Lite::kK0Short) return NaN();
  return Ratio(Norm(Cross(c.decayVertexCm-e.input.primaryVertex.positionCm,c.momentumGeV)),Norm(c.momentumGeV));
}
}
int DaughterIndex(const CeeV0Lite::Candidate& c,int pdg) {
  int found=-1;
  for (int i=0;i<2;++i) if (std::abs(c.daughterPdg[i])==std::abs(pdg)) {
    if (found>=0) return -1;
    found=i;
  }
  return found;
}
double DcaProton(const CeeV0Lite::Candidate& c) { int i=DaughterIndex(c,2212); return i<0 ? NaN() : c.daughterDcaToPvCm[i]; }
double DcaPion(const CeeV0Lite::Candidate& c) { int i=DaughterIndex(c,211); return i<0 ? NaN() : c.daughterDcaToPvCm[i]; }
double Chi2PrimaryProton(const CeeV0Lite::Candidate& c) { int i=DaughterIndex(c,2212); return i<0 ? NaN() : c.daughterChi2FromPv[i]; }
double Chi2PrimaryPion(const CeeV0Lite::Candidate& c) { int i=DaughterIndex(c,211); return i<0 ? NaN() : c.daughterChi2FromPv[i]; }
const std::vector<CandidateObservable>& CandidateObservables() {
  static const std::vector<CandidateObservable> values = {
    {"candidateType", "1", "identity", "Candidate identity / original input track ID", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.type; }},
    {"trackId0", "1", "identity", "Candidate identity / original input track ID", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterId[0]; }},
    {"trackId1", "1", "identity", "Candidate identity / original input track ID", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterId[1]; }},
    {"daughterPdg0", "1", "identity", "Candidate identity / original input track ID", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterPdg[0]; }},
    {"daughterPdg1", "1", "identity", "Candidate identity / original input track ID", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterPdg[1]; }},
    {"massGeV", "GeV/c^2", "original", "Original mother invariant mass", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.massGeV; }},
    {"sigmaMassGeV", "GeV/c^2", "original", "Propagated mass uncertainty", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.sigmaMassGeV; }},
    {"daughterDcaCm", "cm", "original", "Daughter-to-daughter DCA in common-vertex fit", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterDcaCm; }},
    {"geometricChi2", "1", "original", "Common-vertex fit raw chi2", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.geometricChi2; }},
    {"geometricNdf", "1", "original", "Common-vertex fit NDF", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.geometricNdf; }},
    {"geometricChi2Ndf", "1", "original", "Common-vertex chi2 / NDF", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.geometricChi2Ndf; }},
    {"topologicalChi2", "1", "original", "Legacy geometric + PV compatibility chi2 (not 8D refit)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.topologicalChi2; }},
    {"topologicalNdf", "1", "original", "Legacy geometric NDF + 2", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.topologicalNdf; }},
    {"topologicalChi2Ndf", "1", "original", "Legacy combined chi2 / NDF", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.topologicalChi2Ndf; }},
    {"primaryVertexConstraintChi2", "1", "original", "Original mother PV compatibility chi2", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.primaryVertexConstraintChi2; }},
    {"primaryVertexConstraintNdf", "1", "original", "Original mother PV compatibility NDF", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.primaryVertexConstraintNdf; }},
    {"primaryVertexConstraintChi2Ndf", "1", "original", "Original mother PV compatibility chi2 / NDF", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.primaryVertexConstraintChi2Ndf; }},
    {"decayLengthCm", "cm", "original", "Unsigned straight distance |SV-PV|", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.decayLengthCm; }},
    {"sigmaDecayLengthCm", "cm", "original", "Uncertainty of original straight decay length", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.sigmaDecayLengthCm; }},
    {"decayLengthSignificance", "1", "original", "Original unsigned L / sigmaL", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.decayLengthSignificance; }},
    {"pointingCosine", "1", "original", "cos(angle of SV-PV and mother momentum)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.pointingCosine; }},
    {"DcaProton", "cm", "original", "3D DCA of proton/antiproton input state to PV", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return DcaProton(c); }},
    {"DcaPion", "cm", "original", "3D DCA of the unique pion input state to PV", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return DcaPion(c); }},
    {"Chi2PrimaryProton", "1", "original", "Proton/antiproton raw chi2 to PV (not sqrt)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return Chi2PrimaryProton(c); }},
    {"Chi2PrimaryPion", "1", "original", "Unique pion raw chi2 to PV (not sqrt)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return Chi2PrimaryPion(c); }},
    {"motherDcaToPvCm", "cm", "original", "Neutral original mother line DCA to PV; derived, no refit", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return MotherDca(c,e); }},
    {"decayTimeNs", "ns", "original", "Estimated decay time for TOF fixed-point loop", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.decayTimeNs; }},
    {"daughterDcaToPv0Cm", "cm", "original", "Input daughter state 3D DCA to PV", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterDcaToPvCm[0]; }},
    {"daughterChi2FromPv0", "1", "original", "Input daughter raw chi2 to PV", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterChi2FromPv[0]; }},
    {"daughterPvSignificance0", "1", "original", "sqrt(daughter chi2 to PV)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterPvSignificance[0]; }},
    {"daughter0PxGeV", "GeV/c", "original", "Daughter momentum AFTER common-vertex update", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterAtVertex[0].momentumGeV[0]; }},
    {"daughter0PyGeV", "GeV/c", "original", "Daughter momentum AFTER common-vertex update", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterAtVertex[0].momentumGeV[1]; }},
    {"daughter0PzGeV", "GeV/c", "original", "Daughter momentum AFTER common-vertex update", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterAtVertex[0].momentumGeV[2]; }},
    {"daughterDcaToPv1Cm", "cm", "original", "Input daughter state 3D DCA to PV", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterDcaToPvCm[1]; }},
    {"daughterChi2FromPv1", "1", "original", "Input daughter raw chi2 to PV", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterChi2FromPv[1]; }},
    {"daughterPvSignificance1", "1", "original", "sqrt(daughter chi2 to PV)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterPvSignificance[1]; }},
    {"daughter1PxGeV", "GeV/c", "original", "Daughter momentum AFTER common-vertex update", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterAtVertex[1].momentumGeV[0]; }},
    {"daughter1PyGeV", "GeV/c", "original", "Daughter momentum AFTER common-vertex update", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterAtVertex[1].momentumGeV[1]; }},
    {"daughter1PzGeV", "GeV/c", "original", "Daughter momentum AFTER common-vertex update", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.daughterAtVertex[1].momentumGeV[2]; }},
    {"svXCm", "cm", "original", "Unconstrained secondary vertex", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.decayVertexCm[0]; }},
    {"pvXCm", "cm", "event", "Input primary vertex", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return e.input.primaryVertex.positionCm[0]; }},
    {"pxGeV", "GeV/c", "original", "Original mother momentum", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.momentumGeV[0]; }},
    {"svYCm", "cm", "original", "Unconstrained secondary vertex", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.decayVertexCm[1]; }},
    {"pvYCm", "cm", "event", "Input primary vertex", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return e.input.primaryVertex.positionCm[1]; }},
    {"pyGeV", "GeV/c", "original", "Original mother momentum", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.momentumGeV[1]; }},
    {"svZCm", "cm", "original", "Unconstrained secondary vertex", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.decayVertexCm[2]; }},
    {"pvZCm", "cm", "event", "Input primary vertex", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return e.input.primaryVertex.positionCm[2]; }},
    {"pzGeV", "GeV/c", "original", "Original mother momentum", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.momentumGeV[2]; }},
    {"energyGeV", "GeV", "original", "Original mother energy", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.mother.energyGeV; }},
    {"pGeV", "GeV/c", "original", "Original mother momentum magnitude", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return Norm(c.momentumGeV); }},
    {"ptGeV", "GeV/c", "original", "Transverse momentum relative to beam z", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return std::hypot(c.momentumGeV[0],c.momentumGeV[1]); }},
    {"rapidity", "1", "original", "Rapidity along beam z", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return Rapidity(c.mother.energyGeV,c.momentumGeV[2]); }},
    {"motherCovarianceValid", "1", "status", "7D covariance valid flag", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.mother.covarianceValid; }},
    {"motherCovarianceStatus", "1", "status", "7D covariance status enum", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.mother.covarianceStatus; }},
    {"motherJacobianComplete", "1", "status", "All numerical derivative refits succeeded", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.motherJacobianComplete; }},
    {"covarianceDependentCutsEvaluated", "1", "status", "Covariance-dependent legacy cuts actually evaluated", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.covarianceDependentCutsEvaluated; }},
    {"tofIterations", "1", "status", "TOF outer iteration count", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.tofIterations; }},
    {"tofConverged", "1", "status", "TOF outer convergence flag", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.tofConverged; }},
    {"motherTopoValid", "1", "status", "8D PV fit valid flag", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.motherTopo.valid; }},
    {"motherTopoStatus", "1", "status", "ProductionFitStatus enum", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.motherTopo.status; }},
    {"motherTopoCorrelationModel", "1", "status", "PV correlation model enum", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.motherTopo.correlationModel; }},
    {"motherTopoAtProductionVertex", "1", "status", "Stored 8D state endpoint (normally PV)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.motherTopo.atProductionVertex; }},
    {"motherTopoIterations", "1", "status", "8D constraint iteration count", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.motherTopo.iterations; }},
    {"motherTopoCharge", "1", "status", "Mother charge in units e", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.motherTopo.charge; }},
    {"topoFlightValid", "1", "status", "Derived 8D flight observables valid", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return c.topoFlight.valid; }},
    {"motherTopoByTesla", "T", "PV refit", "Field carried by constrained state", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.magneticFieldByTesla : NaN()); }},
    {"motherTopoPvChi2", "1", "PV refit", "8D PV constraint raw chi2", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.pvConstraintChi2 : NaN()); }},
    {"motherTopoPvNdf", "1", "PV refit", "8D PV constraint NDF", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.pvConstraintNdf : NaN()); }},
    {"motherTopoPvChi2Ndf", "1", "PV refit", "8D PV constraint chi2 / NDF", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? Ratio(c.motherTopo.pvConstraintChi2,c.motherTopo.pvConstraintNdf) : NaN()); }},
    {"motherTopoChi2", "1", "PV refit", "8D geometric + PV constraint chi2", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.totalChi2 : NaN()); }},
    {"motherTopoNdf", "1", "PV refit", "8D total NDF", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.totalNdf : NaN()); }},
    {"motherTopoChi2Ndf", "1", "PV refit", "8D total chi2 / NDF", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? Ratio(c.motherTopo.totalChi2,c.motherTopo.totalNdf) : NaN()); }},
    {"flightS", "cm/(GeV/c)", "PV refit", "Signed flight parameter S = Larc / |p|", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.flightParameterS : NaN()); }},
    {"sigmaS", "cm/(GeV/c)", "PV refit", "sqrt(C_SS), includes fitted correlations", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? SafeSqrt(c.motherTopo.covariance.m[7][7]) : NaN()); }},
    {"SOverSigmaS", "1", "PV refit", "Signed S significance", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? Ratio(c.motherTopo.flightParameterS,SafeSqrt(c.motherTopo.covariance.m[7][7])) : NaN()); }},
    {"topoEnergyGeV", "GeV", "PV refit", "PV-constrained mother energy", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.energyGeV : NaN()); }},
    {"topoPxGeV", "GeV/c", "PV refit", "Constrained momentum at stored endpoint", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.momentumGeV[0] : NaN()); }},
    {"topoXCm", "cm", "PV refit", "Constrained position at stored endpoint (normally PV, not SV)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.positionCm[0] : NaN()); }},
    {"topoPyGeV", "GeV/c", "PV refit", "Constrained momentum at stored endpoint", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.momentumGeV[1] : NaN()); }},
    {"topoYCm", "cm", "PV refit", "Constrained position at stored endpoint (normally PV, not SV)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.positionCm[1] : NaN()); }},
    {"topoPzGeV", "GeV/c", "PV refit", "Constrained momentum at stored endpoint", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.momentumGeV[2] : NaN()); }},
    {"topoZCm", "cm", "PV refit", "Constrained position at stored endpoint (normally PV, not SV)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid ? c.motherTopo.positionCm[2] : NaN()); }},
    {"topoMassGeV", "GeV/c^2", "PV refit", "PV-constrained invariant mass", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid && c.topoFlight.valid ? c.topoFlight.massGeV : NaN()); }},
    {"topoSigmaMassGeV", "GeV/c^2", "PV refit", "PV-constrained mass uncertainty", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid && c.topoFlight.valid ? c.topoFlight.sigmaMassGeV : NaN()); }},
    {"topoLCm", "cm", "PV refit", "Signed arc length S |p|", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid && c.topoFlight.valid ? c.topoFlight.decayLengthCm : NaN()); }},
    {"topoSigmaLCm", "cm", "PV refit", "Full 8D propagated arc-length uncertainty, including p-S cross terms", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid && c.topoFlight.valid ? c.topoFlight.sigmaDecayLengthCm : NaN()); }},
    {"LOverSigmaL", "1", "PV refit", "Signed arc L / sigmaL after PV constraint", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid && c.topoFlight.valid ? Ratio(c.topoFlight.decayLengthCm,c.topoFlight.sigmaDecayLengthCm) : NaN()); }},
    {"ctauCm", "cm", "PV refit", "Proper decay length S m (not time in ns)", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid && c.topoFlight.valid ? c.topoFlight.ctauCm : NaN()); }},
    {"sigmaCtauCm", "cm", "PV refit", "Full covariance propagated proper-length uncertainty", [](const CeeV0Lite::Candidate& c,const AnalysisEvent& e)->double { (void)c; (void)e; return (c.motherTopo.valid && c.topoFlight.valid ? c.topoFlight.sigmaCtauCm : NaN()); }},
  };
  return values;
}
const CandidateObservable* FindCandidateObservable(const std::string& name) {
  const std::vector<CandidateObservable>& all=CandidateObservables();
  for (std::size_t i=0;i<all.size();++i) if (all[i].name==name) return &all[i];
  return 0;
}
double CandidateValue(const std::string& name,const CeeV0Lite::Candidate& c,const AnalysisEvent& e) {
  const CandidateObservable* o=FindCandidateObservable(name);
  if (!o) throw std::invalid_argument("unknown candidate variable: "+name);
  return o->get(c,e);
}
}
