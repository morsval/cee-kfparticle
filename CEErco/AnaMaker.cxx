// AnaMaker: TPC + MWDC + iTOF/eTOF Lambda reconstruction, full-reco input version
//
// Formal reconstruction inputs in this version:
//   event vertex       = reconstructed CeeTpcVertex (kVertexMode = 0)
//   TPC position       = reconstructed CeeTpcTrack::GetX/GetY/GetZ
//   TPC direction      = reconstructed CeeTpcTrack::GetPx/GetPy/GetPz
//   MWDC position      = reconstructed CeeMWDC_Track::r1 (mm -> cm)
//   MWDC direction     = reconstructed CeeMWDC_Track::p
//   iTOF/eTOF pos/time = reconstructed CeeiTOFHit/CeeEtofHit position and time
//   proton PID         = daughter DCA topology cut, with normal TOF matching
//   pion PID           = daughter DCA topology cut; if pion uses iTOF, a loose
//                        5 cm iTOF matching cut is applied. If pion uses eTOF,
//                        the source-dependent eTOF matching cut is still applied.
//
// MCTrack information is retained only for truth labels, mctruth/mctruth_hit,
// diagnostics, and optional comparison branches; it is not used as a formal
// reconstruction cut.
//
// Pair selection in this version:
//   keep pair_type == 0 : TPC proton  + TPC pion
//   keep pair_type == 1 : MWDC proton + TPC pion
//   reject pair_type == 2 : TPC proton  + MWDC pion
//   reject pair_type == 3 : MWDC proton + MWDC pion

// c++ class
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>

// root class Header
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TList.h"
#include "TLorentzVector.h"
#include "TProfile.h"
#include "TMath.h"

// FairRoot class
#include "FairEventHeader.h"
#include "FairLink.h"

// CeeRoot class headers
#include "CeeMCTrack.h"
#include "CeeTpcPoint.h"
#include "CeeMWDCPoint.h"
#include "CeeiTOFPoint.h"
#include "CeeeTOFPoint.h"
#include "CeeZdcPoint.h"
#include "CeeT0Hit.h"
#include "CeeTpcHit.h"
#include "CeeiTOFHit.h"
#include "CeeEtofHit.h"
#include "CeeZdcHit.h"
#include "CeeTpcTrack.h"
#include "CeeTpcVertex.h"
#include "CeeMWDC_Hit.h"
#include "CeeMWDC_Track.h"

#include "MyTrack.h"
#include "MyLambda.h"
#include "AnaMaker.h"
#include "TRandom.h"

using namespace std;

const double pmass  = 0.938272;
const double pimass = 0.13957;

// TPC-MWDC matching plane.
// Unit: cm. TPC is projected with MCTrack start position plus the
// reconstructed CeeTpcTrack direction; MWDC uses reconstructed r1 and p.
const Double_t kTpcMwdcMatchZ = 80.0;

// ============================================================
// Reconstruction mode switches
// ============================================================
// DCA-based daughter assignment; MC PDG/mother are retained only for truth tags.
const Int_t kPidMode = 5; // 5 = DCA PID + full-reco inputs + pion-iTOF 5 cm matching cut

// Full-reco default: all MC replacement switches are false.
// TPC position/momentum and MWDC position/momentum use reconstructed values.
const bool kUseTpcMCPos       = false;  // false: TPC position = CeeTpcTrack::GetX/GetY/GetZ
const bool kUseTpcMCMom       = false;  // false: TPC momentum = CeeTpcTrack::GetPx/GetPy/GetPz
const bool kUseMwdcMCPos      = false;  // false: MWDC position = CeeMWDC_Track::r1 * 0.1
const bool kUseMwdcMCMom      = false;  // false: MWDC momentum = CeeMWDC_Track::p

// Full-reco TOF default: reconstructed hit position and reconstructed hit time.
const bool kUseMCTofPointPos  = false;  // false: TOF position = CeeiTOFHit/CeeEtofHit position
const bool kUseMCTofPointTime = false;  // false: TOF time = CeeiTOFHit/CeeEtofHit time

// Common scan presets:
// V0 all MC baseline:
//   kUseTpcMCPos=true,  kUseTpcMCMom=true,  kUseMwdcMCPos=true,  kUseMwdcMCMom=true,
//   kUseMCTofPointPos=true,  kUseMCTofPointTime=true
// V1 only MWDC momentum reco:
//   kUseMwdcMCMom=false, all other switches true
// V2 only MWDC position reco:
//   kUseMwdcMCPos=false, all other switches true
// V3 MWDC position+momentum reco:
//   kUseMwdcMCPos=false, kUseMwdcMCMom=false, all other switches true
// V4 only TOF position reco:
//   kUseMCTofPointPos=false, all other switches true
// V5 only TOF time reco:
//   kUseMCTofPointTime=false, all other switches true
// V8 all reco but still pure-link:
//   all six switches false, kRequireTrackHitLink=true, kRequireHitPDG=true

// Truth purity requirements for this clean version.
const bool kRequireTrackHitLink = false; // require selected hit linkTrack == track mcId
const bool kRequireHitPDG       = false; // require selected hit linked MC PDG equals assumed daughter PDG

// Candidate-pair selection for this version.
// Keep only:
//   pair_type == 0 : TPC proton  + TPC pion
//   pair_type == 1 : MWDC proton + TPC pion
// Reject:
//   pair_type == 2 : TPC proton  + MWDC pion
//   pair_type == 3 : MWDC proton + MWDC pion
const bool kKeepOnlyTpcTpcAndMwdcPTpcPi = true;

// Event vertex mode. Full-reco default: use CeeTpcVertex.
const int kVertexMode = 0;
// 0 = reco vertex from CeeTpcVertex
// 1 = fixed target vertex
// 2 = fixed target vertex + Gaussian smear

// TOF projection cuts
// source: 0 = TPC, 1 = MWDC
// detType: 1 = iTOF, 2 = eTOF
// Control-variable version:
// keep TPC-iTOF, TPC-eTOF and MWDC-iTOF the same as before,
// only tighten MWDC-eTOF from the old 30 cm to 5 cm.
const double ToFMatchCut_TPC_iTOF  = 1.0e9;
const double ToFMatchCut_TPC_eTOF  = 1.0e9;
const double ToFMatchCut_MWDC_iTOF = 1.0e9;
const double ToFMatchCut_MWDC_eTOF = 1.0e9;

// Pion-specific loose iTOF matching cut.
// This is only applied when the selected pion TOF detector is iTOF.
// eTOF pion matching still uses the source-dependent eTOF cut above.
const double PionItofMatchCut = 5.0;

// Keep the old names for compatibility with old comments / quick tests.
const double ToFMatchCut_iTOF = ToFMatchCut_TPC_iTOF;
const double ToFMatchCut_eTOF = ToFMatchCut_TPC_eTOF;

// DCA-based PID cut
const double DaughterDcaCut_TPC  = 1.2;
const double DaughterDcaCut_MWDC = 5.0;

static Double_t GetDaughterDcaCutBySource(Int_t source)
{
    if (source == 1) return DaughterDcaCut_MWDC;
    return DaughterDcaCut_TPC;
}

// ============================================================
// Extra diagnostic branch arrays.
// They are file-scope arrays, so AnaMaker.h does not need to be changed.
// ============================================================
static const Int_t kMaxLambdaDiag = 500;

static Int_t   g_mctruth_hit[kMaxLambdaDiag];   // mctruth==1 and both selected TOF hits match their track MC ids
static Int_t   g_pid_mode[kMaxLambdaDiag];
static Int_t   g_pair_type[kMaxLambdaDiag];        // 0=TPC-TPC, 1=MWDC-TPC, 2=TPC-MWDC, 3=MWDC-MWDC
static Int_t   g_same_pid_mcid[kMaxLambdaDiag];
static Int_t   g_p_pid_mcid[kMaxLambdaDiag];
static Int_t   g_pi_pid_mcid[kMaxLambdaDiag];
static Int_t   g_p_mother_id[kMaxLambdaDiag];
static Int_t   g_pi_mother_id[kMaxLambdaDiag];
static Int_t   g_p_mother_pdg[kMaxLambdaDiag];
static Int_t   g_pi_mother_pdg[kMaxLambdaDiag];

static Int_t   g_p_tof_det[kMaxLambdaDiag];        // 1=iTOF, 2=eTOF
static Int_t   g_pi_tof_det[kMaxLambdaDiag];
static Int_t   g_p_tof_link[kMaxLambdaDiag];
static Int_t   g_pi_tof_link[kMaxLambdaDiag];
static Int_t   g_p_hit_pdg[kMaxLambdaDiag];
static Int_t   g_pi_hit_pdg[kMaxLambdaDiag];
static Int_t   g_p_hit_mother_id[kMaxLambdaDiag];
static Int_t   g_pi_hit_mother_id[kMaxLambdaDiag];
static Int_t   g_p_hit_mother_pdg[kMaxLambdaDiag];
static Int_t   g_pi_hit_mother_pdg[kMaxLambdaDiag];
static Int_t   g_p_track_hit_match[kMaxLambdaDiag];
static Int_t   g_pi_track_hit_match[kMaxLambdaDiag];
static Int_t   g_p_tof_correct[kMaxLambdaDiag];    // selected hit linkTrack == track mcId
static Int_t   g_pi_tof_correct[kMaxLambdaDiag];
static Float_t g_p_tof_dca[kMaxLambdaDiag];
static Float_t g_pi_tof_dca[kMaxLambdaDiag];
static Float_t g_p_tof_dt[kMaxLambdaDiag];         // reco hit time - matched point time
static Float_t g_pi_tof_dt[kMaxLambdaDiag];

// MWDC-proton reco-MC comparison
static Int_t   g_mwdc_p_mcid[kMaxLambdaDiag];
static Int_t   g_mwdc_p_pdg[kMaxLambdaDiag];
static Int_t   g_mwdc_p_mother_pdg[kMaxLambdaDiag];
static Float_t g_mwdc_p_purity[kMaxLambdaDiag];
static Float_t g_mwdc_p_dx[kMaxLambdaDiag];
static Float_t g_mwdc_p_dy[kMaxLambdaDiag];
static Float_t g_mwdc_p_dz[kMaxLambdaDiag];
static Float_t g_mwdc_p_dpx[kMaxLambdaDiag];
static Float_t g_mwdc_p_dpy[kMaxLambdaDiag];
static Float_t g_mwdc_p_dpz[kMaxLambdaDiag];
static Float_t g_mwdc_p_dangle[kMaxLambdaDiag];
static Float_t g_mwdc_p_line_dca_mcstart[kMaxLambdaDiag];
static Float_t g_mwdc_p_dca_reco[kMaxLambdaDiag];
static Float_t g_mwdc_p_dca_mc[kMaxLambdaDiag];

// MWDC-pion reco-MC comparison
static Int_t   g_mwdc_pi_mcid[kMaxLambdaDiag];
static Int_t   g_mwdc_pi_pdg[kMaxLambdaDiag];
static Int_t   g_mwdc_pi_mother_pdg[kMaxLambdaDiag];
static Float_t g_mwdc_pi_purity[kMaxLambdaDiag];
static Float_t g_mwdc_pi_dx[kMaxLambdaDiag];
static Float_t g_mwdc_pi_dy[kMaxLambdaDiag];
static Float_t g_mwdc_pi_dz[kMaxLambdaDiag];
static Float_t g_mwdc_pi_dpx[kMaxLambdaDiag];
static Float_t g_mwdc_pi_dpy[kMaxLambdaDiag];
static Float_t g_mwdc_pi_dpz[kMaxLambdaDiag];
static Float_t g_mwdc_pi_dangle[kMaxLambdaDiag];
static Float_t g_mwdc_pi_line_dca_mcstart[kMaxLambdaDiag];
static Float_t g_mwdc_pi_dca_reco[kMaxLambdaDiag];
static Float_t g_mwdc_pi_dca_mc[kMaxLambdaDiag];

//add mwdc res check
// ============================================================
// Independent MWDC track resolution tree
// One entry corresponds to one reconstructed MWDC track.
// These variables are file-scope, so AnaMaker.h does not need
// to be modified.
// ============================================================
static TTree* gMwdcTrackTree = nullptr;

// Event and track identity
static Int_t g_mwdc_evt_id       = -1;
static Int_t g_mwdc_track_index  = -1;
static Int_t g_mwdc_track_mcid   = -999;

// MC particle information
static Int_t g_mwdc_track_pdg        = -999999;
static Int_t g_mwdc_track_mother_id  = -999999;
static Int_t g_mwdc_track_mother_pdg = -999999;
static Int_t g_mwdc_track_is_primary = 0;

// MWDC reconstruction quality
static Int_t   g_mwdc_track_hit_num       = -1;
static Int_t   g_mwdc_track_num_true_hit  = -1;
static Float_t g_mwdc_track_purity        = -1.0;

// Reconstructed and MC momentum
static Float_t g_mwdc_track_p_reco     = -999.0;
static Float_t g_mwdc_track_p_mc       = -999.0;
static Float_t g_mwdc_track_theta_reco = -999.0;
static Float_t g_mwdc_track_theta_mc   = -999.0;

// Main resolution variables
static Float_t g_mwdc_track_angle_res = -999.0;
static Float_t g_mwdc_track_mom_res   = -999.0;

// ============================================================
// Independent TPC-MWDC matching tree.
// One entry corresponds to one TPC-track / MWDC-track pair.
// Only the minimum information needed for a z=80 cm cut scan is saved.
// ============================================================
static TTree* gTpcMwdcMatchTree = nullptr;

static Int_t   g_match_evt_id      = -1;
static Int_t   g_match_tpc_index   = -1;
static Int_t   g_match_mwdc_index  = -1;
static Int_t   g_match_tpc_mcid    = -999;
static Int_t   g_match_mwdc_mcid   = -999;
static Int_t   g_match_same_mcid   = 0;
static Float_t g_match_dx_z80      = -999.0;
static Float_t g_match_dy_z80      = -999.0;
static Float_t g_match_dr_z80      = -999.0;

//========res over


static Double_t NormalizeDeltaPhiDeg(Double_t dphi)
{
    dphi = std::fmod(dphi + 180.0, 360.0);
    if (dphi < 0) dphi += 360.0;
    return dphi - 180.0;
}

// ============================================================
// Straight-line projection to a fixed z plane.
// The returned x and y are in cm.
// ============================================================
static Bool_t ProjectTrackToZPlane(const TVector3& pos,
                                   const TVector3& mom,
                                   Double_t zPlane,
                                   TVector3& proj)
{
    proj.SetXYZ(-9999., -9999., zPlane);

    if (!std::isfinite(pos.X()) ||
        !std::isfinite(pos.Y()) ||
        !std::isfinite(pos.Z()) ||
        !std::isfinite(mom.X()) ||
        !std::isfinite(mom.Y()) ||
        !std::isfinite(mom.Z())) {
        return kFALSE;
    }

    if (fabs(mom.Z()) < 1.0e-12) return kFALSE;

    const Double_t scale = (zPlane - pos.Z()) / mom.Z();
    const Double_t x = pos.X() + scale * mom.X();
    const Double_t y = pos.Y() + scale * mom.Y();

    if (!std::isfinite(x) || !std::isfinite(y)) return kFALSE;

    proj.SetXYZ(x, y, zPlane);
    return kTRUE;
}

// ============================================================
// TOF projection helpers
// ============================================================
static Double_t GetTrackItofProjByModulePlane(const TVector3& pos,
                                              const TVector3& mom,
                                              const TVector3& hit,
                                              Int_t det,
                                              TVector3& proj)
{
    proj.SetXYZ(-9999., -9999., -9999.);
    const Double_t x0 = pos.X();
    const Double_t y0 = pos.Y();
    const Double_t z0 = pos.Z();
    const Double_t px = mom.X();
    const Double_t py = mom.Y();
    const Double_t pz = mom.Z();

    Double_t s = -9999.0;

    // side iTOF, use y-z plane
    if ((det >= 0 && det <= 8) || (det >= 15 && det <= 24)) {
        if (fabs(px) < 1.0e-12) return 1.0e9;
        s = (hit.X() - x0) / px;
        if (fabs(s) < 1.0e-12) return 1.0e9;
        proj.SetXYZ(x0 + s * px, y0 + s * py, z0 + s * pz);
        const Double_t dy = hit.Y() - proj.Y();
        const Double_t dz = hit.Z() - proj.Z();
        return sqrt(dy * dy + dz * dz);
    }

    // forward iTOF, use x-y plane
    if (det >= 9 && det <= 14) {
        if (fabs(pz) < 1.0e-12) return 1.0e9;
        s = (hit.Z() - z0) / pz;
        if (fabs(s) < 1.0e-12) return 1.0e9;
        proj.SetXYZ(x0 + s * px, y0 + s * py, z0 + s * pz);
        const Double_t dx = hit.X() - proj.X();
        const Double_t dy = hit.Y() - proj.Y();
        return sqrt(dx * dx + dy * dy);
    }

    // unknown det: fallback to DCA in xy
    const Double_t dx = pos.X() - hit.X();
    const Double_t dy = pos.Y() - hit.Y();
    const Double_t pMag2 = px * px + py * py;
    if (pMag2 < 1.0e-12) return 1.0e9;
    const Double_t cross = fabs(dx * py - dy * px);
    return cross / sqrt(pMag2);
}

static Double_t GetTrackEtofProjByPlane(const TVector3& pos,
                                        const TVector3& mom,
                                        const TVector3& hit,
                                        TVector3& proj)
{
    proj.SetXYZ(-9999., -9999., -9999.);
    const Double_t x0 = pos.X();
    const Double_t y0 = pos.Y();
    const Double_t z0 = pos.Z();
    const Double_t px = mom.X();
    const Double_t py = mom.Y();
    const Double_t pz = mom.Z();

    if (fabs(pz) < 1.0e-12) return 1.0e9;
    const Double_t s = (hit.Z() - z0) / pz;
    if (fabs(s) < 1.0e-12) return 1.0e9;

    proj.SetXYZ(x0 + s * px, y0 + s * py, z0 + s * pz);
    const Double_t dx = hit.X() - proj.X();
    const Double_t dy = hit.Y() - proj.Y();
    return sqrt(dx * dx + dy * dy);
}

struct TofMatchResult
{
    Bool_t found;
    Int_t detType;       // 0=none, 1=iTOF, 2=eTOF
    Double_t dca;
    Double_t mcDca;
    TVector3 hitPos;     // reconstructed hit position
    TVector3 mcHitPos;   // linked MC point position
    Double_t time;       // reconstructed hit time
    Double_t mcTime;     // linked MC point time
    Double_t mcLength;
    Bool_t hasPointInfo;
    Int_t linkTrack;
    Int_t linkPoint;

    TofMatchResult()
        : found(kFALSE), detType(0), dca(1.0e9), mcDca(1.0e9),
          hitPos(0., 0., 0.), mcHitPos(0., 0., 0.),
          time(0.), mcTime(0.), mcLength(0.), hasPointInfo(kFALSE),
          linkTrack(-999), linkPoint(-999) {}
};

static Bool_t IsGoodArrayIndex(TClonesArray* arr, Int_t idx)
{
    return (arr && idx >= 0 && idx < arr->GetEntries());
}

static TofMatchResult FindNearestItofHit(const TVector3& trkPos,
                                         const TVector3& trkMom,
                                         Int_t mcId,
                                         TClonesArray* itofHits,
                                         TClonesArray* itofPoints)
{
    TofMatchResult out;
    out.detType = 1;
    if (!itofHits) return out;

    for (Int_t ihit = 0; ihit < itofHits->GetEntries(); ++ihit) {
        CeeiTOFHit* hit = dynamic_cast<CeeiTOFHit*>(itofHits->At(ihit));
        if (!hit) continue;

        TVector3 hitPos(hit->GetX(), hit->GetY(), hit->GetZ());
        Int_t det_tmp = hit->GetDetectorID() / 100;

        TVector3 projPos;
        Double_t tmp_dca = GetTrackItofProjByModulePlane(trkPos, trkMom, hitPos, det_tmp, projPos);

        if (tmp_dca < out.dca) {
            FairLink link1 = hit->GetLink(0);
            FairLink link2 = hit->GetLink(1);

            out.found = kTRUE;
            out.dca = tmp_dca;
            out.linkTrack = link1.GetIndex();
            out.linkPoint = link2.GetIndex();
            out.time = hit->GetTime();
            out.hitPos.SetXYZ(hit->GetX(), hit->GetY(), hit->GetZ());

            if (out.linkTrack == mcId) out.mcDca = tmp_dca;

            // Fill point time/length whenever point link is valid.
            // This is useful for TPC tracks whose track mcid may be unavailable.
            out.mcTime = 0.;
            out.mcLength = 0.;
            if (IsGoodArrayIndex(itofPoints, out.linkPoint)) {
                CeeiTOFPoint* point = dynamic_cast<CeeiTOFPoint*>(itofPoints->At(out.linkPoint));
                if (point) {
                    out.mcTime = point->GetTime();
                    out.mcLength = point->GetLength();
                    out.mcHitPos.SetXYZ(point->GetX(), point->GetY(), point->GetZ());
                    out.hasPointInfo = kTRUE;
                }
            }
        }
    }
    return out;
}

static TofMatchResult FindNearestEtofHit(const TVector3& trkPos,
                                         const TVector3& trkMom,
                                         Int_t mcId,
                                         TClonesArray* etofHits,
                                         TClonesArray* etofPoints)
{
    TofMatchResult out;
    out.detType = 2;
    if (!etofHits) return out;

    for (Int_t ihit = 0; ihit < etofHits->GetEntries(); ++ihit) {
        CeeEtofHit* hit = dynamic_cast<CeeEtofHit*>(etofHits->At(ihit));
        if (!hit) continue;

        TVector3 hitPos(hit->GetX(), hit->GetY(), hit->GetZ());
        TVector3 projPos;
        Double_t tmp_dca = GetTrackEtofProjByPlane(trkPos, trkMom, hitPos, projPos);

        if (tmp_dca < out.dca) {
            FairLink link1 = hit->GetLink(0);
            FairLink link2 = hit->GetLink(1);

            out.found = kTRUE;
            out.dca = tmp_dca;
            out.linkTrack = link1.GetIndex();
            out.linkPoint = link2.GetIndex();
            out.time = hit->GetTime();
            out.hitPos.SetXYZ(hit->GetX(), hit->GetY(), hit->GetZ());

            if (out.linkTrack == mcId) out.mcDca = tmp_dca;

            out.mcTime = 0.;
            out.mcLength = 0.;
            if (IsGoodArrayIndex(etofPoints, out.linkPoint)) {
                CeeeTOFPoint* point = dynamic_cast<CeeeTOFPoint*>(etofPoints->At(out.linkPoint));
                if (point) {
                    out.mcTime = point->GetTime();
                    out.mcLength = point->GetLength();
                    out.mcHitPos.SetXYZ(point->GetX(), point->GetY(), point->GetZ());
                    out.hasPointInfo = kTRUE;
                }
            }
        }
    }
    return out;
}

static Double_t GetTofMatchCutBySourceDet(Int_t source, Int_t detType)
{
    // source: 0 = TPC, 1 = MWDC
    // detType: 1 = iTOF, 2 = eTOF
    if (source == 0 && detType == 1) return ToFMatchCut_TPC_iTOF;
    if (source == 0 && detType == 2) return ToFMatchCut_TPC_eTOF;
    if (source == 1 && detType == 1) return ToFMatchCut_MWDC_iTOF;
    if (source == 1 && detType == 2) return ToFMatchCut_MWDC_eTOF;
    return 1.0e9;
}

// Old det-only helper kept for compatibility.
static Double_t GetTofMatchCutByDetType(Int_t detType)
{
    if (detType == 1) return ToFMatchCut_iTOF;
    if (detType == 2) return ToFMatchCut_eTOF;
    return 1.0e9;
}

static TofMatchResult SelectTofPreferItof(const TofMatchResult& itof,
                                          const TofMatchResult& etof,
                                          Int_t source)
{
    Double_t cutItof = GetTofMatchCutBySourceDet(source, 1);
    Double_t cutEtof = GetTofMatchCutBySourceDet(source, 2);

    if (itof.found && fabs(itof.dca) <= cutItof) return itof;
    if (etof.found && fabs(etof.dca) <= cutEtof) return etof;
    if (itof.found && (!etof.found || itof.dca <= etof.dca)) return itof;
    return etof;
}

// Pion-specific TOF selection.
// Prefer the nearest reconstructed iTOF hit if it exists.
// The loose pion-iTOF matching cut is applied downstream.
// If there is no iTOF hit, use eTOF and keep the source-dependent eTOF cut.
static TofMatchResult SelectPionTofPreferItof(const TofMatchResult& itof,
                                              const TofMatchResult& etof)
{
    if (itof.found) return itof;
    return etof;
}

// ============================================================
// MC info and candidate list
// ============================================================
static Bool_t GetMCInfoFromMCID(TClonesArray* mcTracks,
                                Int_t mcid,
                                Int_t& pdg,
                                Int_t& motherId,
                                Int_t& motherPdg,
                                TVector3& mcPos,
                                TVector3& mcMom)
{
    pdg = -999999;
    motherId = -999999;
    motherPdg = -999999;
    mcPos.SetXYZ(0., 0., 0.);
    mcMom.SetXYZ(0., 0., 0.);

    if (!mcTracks) return kFALSE;
    if (mcid < 0 || mcid >= mcTracks->GetEntriesFast()) return kFALSE;

    CeeMCTrack* mc = dynamic_cast<CeeMCTrack*>(mcTracks->At(mcid));
    if (!mc) return kFALSE;

    pdg = mc->GetPdgCode();
    motherId = mc->GetMotherId();
    mcPos.SetXYZ(mc->GetStartX(), mc->GetStartY(), mc->GetStartZ());
    mc->GetMomentum(mcMom);

    if (motherId >= 0 && motherId < mcTracks->GetEntriesFast()) {
        CeeMCTrack* mother = dynamic_cast<CeeMCTrack*>(mcTracks->At(motherId));
        if (mother) motherPdg = mother->GetPdgCode();
    }
    return kTRUE;
}

struct RecoTrackCand
{
    Int_t source;   // 0 = TPC, 1 = MWDC
    Int_t index;
    Int_t mcId;     // MWDC: track_mcid. TPC: unavailable here, use TOF link as fallback.

    TVector3 pos;
    TVector3 mom;
    Double_t dedx;

    Int_t mcPdg;
    Int_t motherId;
    Int_t motherPdg;
    Double_t purity;
    Int_t hitNum;
    Int_t numMcHitTrue;

    TVector3 mcPos;
    TVector3 mcMom;

    RecoTrackCand()
        : source(-1), index(-1), mcId(-999),
          pos(0., 0., 0.), mom(0., 0., 0.), dedx(-999.),
          mcPdg(-999999), motherId(-999999), motherPdg(-999999),
          purity(-1.), hitNum(-1), numMcHitTrue(-1),
          mcPos(0., 0., 0.), mcMom(0., 0., 0.) {}
};

static void BuildRecoTrackCandList(std::vector<RecoTrackCand>& out,
                                   TClonesArray* tpcTracks,
                                   TClonesArray* mwdcTracks,
                                   TClonesArray* mcTracks)
{
    out.clear();

    // TPC tracks are kept to form TPC-TPC and MWDC-proton + TPC-pion pairs.
    if (tpcTracks) {
        for (Int_t i = 0; i < tpcTracks->GetEntries(); ++i) {
            CeeTpcTrack* trk = dynamic_cast<CeeTpcTrack*>(tpcTracks->At(i));
            if (!trk) continue;

            RecoTrackCand c;
            c.source = 0;
            c.index  = i;
            c.mcId   = trk->GetMCTrackID();
            c.pos.SetXYZ(trk->GetX(), trk->GetY(), trk->GetZ());
            c.mom.SetXYZ(trk->GetPx(), trk->GetPy(), trk->GetPz());
            c.dedx = -999.;

            TVector3 mcPos;
            TVector3 mcMom;
            Int_t pdg = -999999;
            Int_t motherId = -999999;
            Int_t motherPdg = -999999;
            if (GetMCInfoFromMCID(mcTracks, c.mcId, pdg, motherId, motherPdg, mcPos, mcMom)) {
                c.mcPdg = pdg;
                c.motherId = motherId;
                c.motherPdg = motherPdg;
                c.mcPos = mcPos;
                c.mcMom = mcMom;

                if (mcMom.Mag() > 1.0e-12) {
                    if (kUseTpcMCPos) c.pos = mcPos;
                    if (kUseTpcMCMom) c.mom = mcMom;
                }
            }

            if (c.mcId < 0) continue;
            if (c.mom.Mag() < 1.0e-12) continue;
            out.push_back(c);
        }
    }

    // MWDC tracks: use reco r,p by default.
    if (mwdcTracks) {
        for (Int_t i = 0; i < mwdcTracks->GetEntries(); ++i) {
            CeeMWDC_Track* trk = dynamic_cast<CeeMWDC_Track*>(mwdcTracks->At(i));
            if (!trk) continue;

            RecoTrackCand c;
            c.source = 1;
            c.index  = i;
            c.mcId   = trk->track_mcid;
            c.hitNum = trk->hit_num;
            c.numMcHitTrue = trk->num_mchit_true;

            if (trk->hit_num > 0) {
                c.purity = double(trk->num_mchit_true) / double(trk->hit_num);
            }

            c.pos = trk->r1 * 0.1; // mm -> cm
            c.mom = trk->p;
            c.dedx = -999.;

            TVector3 mcPos;
            TVector3 mcMom;
            Int_t pdg = -999999;
            Int_t motherId = -999999;
            Int_t motherPdg = -999999;

            if (GetMCInfoFromMCID(mcTracks, c.mcId, pdg, motherId, motherPdg, mcPos, mcMom)) {
                c.mcPdg = pdg;
                c.motherId = motherId;
                c.motherPdg = motherPdg;
                c.mcPos = mcPos;
                c.mcMom = mcMom;

                if (mcMom.Mag() > 1.0e-12) {
                    if (kUseMwdcMCPos) c.pos = mcPos;
                    if (kUseMwdcMCMom) c.mom = mcMom;
                }
            }

            if (c.mom.Mag() < 1.0e-12) continue;
            out.push_back(c);
        }
    }
}

static Int_t GetPairType(Int_t pSource, Int_t piSource)
{
    if (pSource == 1 && piSource == 0) return 1; // MWDC proton + TPC pion
    if (pSource == 0 && piSource == 1) return 2; // TPC proton + MWDC pion
    if (pSource == 1 && piSource == 1) return 3; // MWDC proton + MWDC pion
    return 0;                                    // TPC-TPC
}

static void FillMwdcDiagProton(Int_t n,
                               const RecoTrackCand& cand,
                               const TVector3& pos,
                               const TVector3& mom,
                               Double_t dcaReco,
                               Double_t dcaMC)
{
    if (n < 0 || n >= kMaxLambdaDiag) return;
    if (cand.source != 1 || cand.mcMom.Mag() < 1.0e-12) return;

    g_mwdc_p_mcid[n] = cand.mcId;
    g_mwdc_p_pdg[n] = cand.mcPdg;
    g_mwdc_p_mother_pdg[n] = cand.motherPdg;
    g_mwdc_p_purity[n] = cand.purity;

    g_mwdc_p_dx[n] = pos.X() - cand.mcPos.X();
    g_mwdc_p_dy[n] = pos.Y() - cand.mcPos.Y();
    g_mwdc_p_dz[n] = pos.Z() - cand.mcPos.Z();
    g_mwdc_p_dpx[n] = mom.X() - cand.mcMom.X();
    g_mwdc_p_dpy[n] = mom.Y() - cand.mcMom.Y();
    g_mwdc_p_dpz[n] = mom.Z() - cand.mcMom.Z();
    g_mwdc_p_dangle[n] = mom.Angle(cand.mcMom) * TMath::RadToDeg();
    g_mwdc_p_line_dca_mcstart[n] = (pos - cand.mcPos).Cross(mom.Unit()).Mag();
    g_mwdc_p_dca_reco[n] = dcaReco;
    g_mwdc_p_dca_mc[n] = dcaMC;
}

static void FillMwdcDiagPion(Int_t n,
                             const RecoTrackCand& cand,
                             const TVector3& pos,
                             const TVector3& mom,
                             Double_t dcaReco,
                             Double_t dcaMC)
{
    if (n < 0 || n >= kMaxLambdaDiag) return;
    if (cand.source != 1 || cand.mcMom.Mag() < 1.0e-12) return;

    g_mwdc_pi_mcid[n] = cand.mcId;
    g_mwdc_pi_pdg[n] = cand.mcPdg;
    g_mwdc_pi_mother_pdg[n] = cand.motherPdg;
    g_mwdc_pi_purity[n] = cand.purity;

    g_mwdc_pi_dx[n] = pos.X() - cand.mcPos.X();
    g_mwdc_pi_dy[n] = pos.Y() - cand.mcPos.Y();
    g_mwdc_pi_dz[n] = pos.Z() - cand.mcPos.Z();
    g_mwdc_pi_dpx[n] = mom.X() - cand.mcMom.X();
    g_mwdc_pi_dpy[n] = mom.Y() - cand.mcMom.Y();
    g_mwdc_pi_dpz[n] = mom.Z() - cand.mcMom.Z();
    g_mwdc_pi_dangle[n] = mom.Angle(cand.mcMom) * TMath::RadToDeg();
    g_mwdc_pi_line_dca_mcstart[n] = (pos - cand.mcPos).Cross(mom.Unit()).Mag();
    g_mwdc_pi_dca_reco[n] = dcaReco;
    g_mwdc_pi_dca_mc[n] = dcaMC;
}

// ============================================================
// AnaMaker
// ============================================================
AnaMaker::AnaMaker(const TString& inFileDir, int fileNum, int debug)
    : fOutFileName(""),
      fOutFile(nullptr),
      fChain(nullptr),
      fDebug(debug),
      fProtonTrack(nullptr),
      fPionTrack(nullptr),
      fProtonMCTrack(nullptr),
      fPionMCTrack(nullptr),
      fMCTrack(nullptr),
      fTpcPoint(nullptr),
      fMWDCPoint(nullptr),
      fiTOFPoint(nullptr),
      feTOFPoint(nullptr),
      fZdcPoint(nullptr),
      fT0Hit(nullptr),
      fTpcHit(nullptr),
      fiTOFHit(nullptr),
      fEtofHit(nullptr),
      fZdcHit(nullptr),
      fTpcTrack(nullptr),
      fTpcVertex(nullptr),
      fMWDC_Hit(nullptr),
      fMWDC_Track(nullptr),
      fEvtNo(0),
      fRotPi(false),
      fRotP(false),
      fEvtVtx(TVector3(0, 0, 0))
{
    fMCTrack    = new TClonesArray("CeeMCTrack");
    fTpcPoint   = new TClonesArray("CeeTpcPoint");
    fMWDCPoint  = new TClonesArray("CeeMWDCPoint");
    fiTOFPoint  = new TClonesArray("CeeiTOFPoint");
    feTOFPoint  = new TClonesArray("CeeeTOFPoint");
    fZdcPoint   = new TClonesArray("CeeZdcPoint");
    fT0Hit      = new TClonesArray("CeeT0Hit");
    fTpcHit     = new TClonesArray("CeeTpcHit");
    fiTOFHit    = new TClonesArray("CeeiTOFHit");
    fEtofHit    = new TClonesArray("CeeEtofHit");
    fZdcHit     = new TClonesArray("CeeZdcHit");
    fTpcTrack   = new TClonesArray("CeeTpcTrack");
    fTpcVertex  = new TClonesArray("CeeTpcVertex");
    fMWDC_Hit   = new TClonesArray("CeeMWDC_Hit");
    fMWDC_Track = new TClonesArray("CeeMWDC_Track");

    Init(inFileDir, fileNum);
}

AnaMaker::~AnaMaker()
{
    delete fOutFile;
    delete fChain;

    delete fProtonTrack;
    delete fPionTrack;

    delete mc_pty;
    delete mc_v1;
    delete rec_v1;
    delete tofdca;
    delete mc_tofdca;
    delete itofTime_diff;
    delete itofLength_diff;
}

bool AnaMaker::Init(TString inFileDir, int fileNum)
{
    fChain = new TChain("ceesim");

    if (inFileDir.IsNull()) {
        Error("Init", "Input dir is not defined");
        return false;
    }

    if (fileNum == 1 && inFileDir.EndsWith(".root")) {
        fChain->Add(inFileDir);
        if (fDebug > 0) Info("Init", "Add input file %s", inFileDir.Data());
    }
    else if (inFileDir.EndsWith(".list")) {
        ifstream fileList(inFileDir.Data());
        if (!fileList.is_open()) {
            Error("Init", "Could not open list %s", inFileDir.Data());
            return false;
        }
        string rootFilePath;
        while (std::getline(fileList, rootFilePath)) {
            TString tmpFile(rootFilePath);
            if (tmpFile.EndsWith(".root")) {
                if (fDebug > 0) Info("Init", "Add input file %s in the list %s", tmpFile.Data(), inFileDir.Data());
                fChain->Add(tmpFile);
            }
        }
    }
    else {
        TSystemDirectory* dir = new TSystemDirectory("dir", inFileDir);
        TList* fileList = dir->GetListOfFiles();

        if (!fileList) {
            Error("Init", "Failed to get the input list in %s", inFileDir.Data());
            return false;
        }

        TSystemFile* file;
        TString fileName;
        TIter next(fileList);
        int fileCount = 0;
        while ((file = (TSystemFile*)next())) {
            fileName = file->GetName();
            if (file->IsDirectory() || fileName == "." || fileName == ".." || !fileName.BeginsWith("Fsim_CC1p1"))
                continue;

            TString filePath = TString::Format("%s/%s", inFileDir.Data(), fileName.Data());
            fChain->Add(filePath);
            if (fDebug > 0) Info("Init", "Add input file %s", filePath.Data());

            fileCount++;
            if (fileCount >= fileNum) break;
        }
        Info("Init", "Total %d root files in %s are read :D", fileCount, inFileDir.Data());
        delete dir;
    }

    std::vector<TString> enabledBranches = {
        "MCTrack", "CeeiTOFHit", "CeeEtofHit", "CeeTpcTrack", "MWDC_Track",
        "CeeTpcVertex", "CeeT0Hit", "CeeiTOFPoint", "CeeeTOFPoint"
    };

    TObjArray* branches = fChain->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); i++) {
        TBranch* branch = (TBranch*)branches->At(i);
        TString branchName = branch->GetName();
        bool shouldEnable = (std::find(enabledBranches.begin(), enabledBranches.end(), branchName) != enabledBranches.end());
        fChain->SetBranchStatus(branchName, shouldEnable ? kTRUE : kFALSE);
    }

    if (fEvtNo == 0) fEvtNo = fChain->GetEntries();

    fChain->SetBranchAddress("MCTrack", &fMCTrack);
    fChain->SetBranchAddress("CeeiTOFPoint", &fiTOFPoint);
    fChain->SetBranchAddress("CeeeTOFPoint", &feTOFPoint);
    fChain->SetBranchAddress("CeeT0Hit", &fT0Hit);
    fChain->SetBranchAddress("CeeiTOFHit", &fiTOFHit);
    fChain->SetBranchAddress("CeeEtofHit", &fEtofHit);
    fChain->SetBranchAddress("CeeTpcTrack", &fTpcTrack);
    fChain->SetBranchAddress("CeeTpcVertex", &fTpcVertex);

    if (fChain->GetBranch("MWDC_Track")) {
        fChain->SetBranchAddress("MWDC_Track", &fMWDC_Track);
        Info("Init", "MWDC_Track branch is found and enabled.");
    }
    else {
        Warning("Init", "MWDC_Track branch is NOT found. This diagnostic run needs MWDC_Track.");
    }

    return true;
}

void AnaMaker::InitialHistgram()
{
    if (fOutFileName.IsNull()) fOutFileName = "ana_tpcTpc_and_mwdcPTpcPi_fullReco_pionItof5cm_mctruthHit.root";
    if (fDebug > 0) Info("InitialHistgram", "Output file: %s", fOutFileName.Data());

    gSystem->Exec("mkdir -p output");
    fOutFile = new TFile(TString("output/") + fOutFileName, "RECREATE");
    if (!fOutFile->IsOpen()) {
        std::cerr << "Error: Cannot create output file" << std::endl;
        return;
    }

    fOutFile->cd();
    fRecTree = new TTree("lambda", "reclambda");

    // ============================================================
    // Independent MWDC track resolution tree
    // ============================================================
    gMwdcTrackTree = new TTree("mwdcTrack", "MWDC reconstructed track resolution");

    // Event and track identity
    gMwdcTrackTree->Branch("evt_id", &g_mwdc_evt_id, "evt_id/I");
    gMwdcTrackTree->Branch("track_index", &g_mwdc_track_index, "track_index/I");
    gMwdcTrackTree->Branch("mcid", &g_mwdc_track_mcid, "mcid/I");
    // MC particle information
    gMwdcTrackTree->Branch("pdg", &g_mwdc_track_pdg, "pdg/I");
    gMwdcTrackTree->Branch("mother_id", &g_mwdc_track_mother_id, "mother_id/I");
    gMwdcTrackTree->Branch("mother_pdg", &g_mwdc_track_mother_pdg, "mother_pdg/I");
    gMwdcTrackTree->Branch("is_primary", &g_mwdc_track_is_primary, "is_primary/I");

    // MWDC reconstruction quality
    gMwdcTrackTree->Branch("hit_num", &g_mwdc_track_hit_num, "hit_num/I");

    gMwdcTrackTree->Branch("num_true_hit", &g_mwdc_track_num_true_hit, "num_true_hit/I");

    gMwdcTrackTree->Branch( "purity", &g_mwdc_track_purity, "purity/F");

    // Reconstructed and MC kinematics
    gMwdcTrackTree->Branch("p_reco", &g_mwdc_track_p_reco, "p_reco/F");

    gMwdcTrackTree->Branch("p_mc", &g_mwdc_track_p_mc, "p_mc/F");

    gMwdcTrackTree->Branch("theta_reco", &g_mwdc_track_theta_reco, "theta_reco/F");
    gMwdcTrackTree->Branch("angle_res", &g_mwdc_track_angle_res, "angle_res/F");
    gMwdcTrackTree->Branch("theta_mc", &g_mwdc_track_theta_mc, "theta_mc/F");

    // Main resolution branchesgMwdcTrackTree->Branch("angle_res", &g_mwdc_track_angle_res, "angle_res/F");
    
    gMwdcTrackTree->Branch("mom_res", &g_mwdc_track_mom_res, "mom_res/F");
    //====
    
    // ============================================================
    // Independent TPC-MWDC matching tree at z = 80 cm.
    // No TOF, PID, DCA, Lambda or event-vertex cut is applied here.
    // ============================================================
    gTpcMwdcMatchTree = new TTree(
        "tpcMwdcMatch",
        "TPC-MWDC straight-line matching at z=80 cm"
    );

    gTpcMwdcMatchTree->Branch("evt_id", &g_match_evt_id, "evt_id/I");
    gTpcMwdcMatchTree->Branch("tpc_index", &g_match_tpc_index, "tpc_index/I");
    gTpcMwdcMatchTree->Branch("mwdc_index", &g_match_mwdc_index, "mwdc_index/I");
    gTpcMwdcMatchTree->Branch("tpc_mcid", &g_match_tpc_mcid, "tpc_mcid/I");
    gTpcMwdcMatchTree->Branch("mwdc_mcid", &g_match_mwdc_mcid, "mwdc_mcid/I");
    gTpcMwdcMatchTree->Branch("same_mcid", &g_match_same_mcid, "same_mcid/I");
    gTpcMwdcMatchTree->Branch("dx_z80", &g_match_dx_z80, "dx_z80/F");
    gTpcMwdcMatchTree->Branch("dy_z80", &g_match_dy_z80, "dy_z80/F");
    gTpcMwdcMatchTree->Branch("dr_z80", &g_match_dr_z80, "dr_z80/F");

    fRecTree->Branch("evtId", &evtId, "evtId/I");
    fRecTree->Branch("psi", &psi, "psi/F");
    fRecTree->Branch("nLambda", &nLambda, "nLambda/I");
    fRecTree->Branch("mass", ld_mass, "ld_mass[nLambda]/F");
    fRecTree->Branch("mctruth", ld_mctruth, "ld_mctruth[nLambda]/I");
    fRecTree->Branch("mctruth_hit", g_mctruth_hit, "mctruth_hit[nLambda]/I");
    fRecTree->Branch("dca1to2", dca1to2, "dca1to2[nLambda]/F");
    fRecTree->Branch("decaylength", decaylength, "decaylength[nLambda]/F");
    fRecTree->Branch("dca", ld_dca, "ld_dca[nLambda]/F");
    fRecTree->Branch("px", ld_px, "ld_px[nLambda]/F");
    fRecTree->Branch("py", ld_py, "ld_py[nLambda]/F");
    fRecTree->Branch("pz", ld_pz, "ld_pz[nLambda]/F");
    fRecTree->Branch("mc_px", mc_ld_px, "mc_ld_px[nLambda]/F");
    fRecTree->Branch("mc_py", mc_ld_py, "mc_ld_py[nLambda]/F");
    fRecTree->Branch("mc_pz", mc_ld_pz, "mc_ld_pz[nLambda]/F");

    fRecTree->Branch("pi_pdg", pi_pdg, "pi_pdg[nLambda]/I");
    fRecTree->Branch("pi_dedx", pi_dedx, "pi_dedx[nLambda]/F");
    fRecTree->Branch("pi_dca", pi_dca, "pi_dca[nLambda]/F");
    fRecTree->Branch("pi_px", pi_px, "pi_px[nLambda]/F");
    fRecTree->Branch("pi_py", pi_py, "pi_py[nLambda]/F");
    fRecTree->Branch("pi_pz", pi_pz, "pi_pz[nLambda]/F");
    fRecTree->Branch("mc_pi_px", mc_pi_px, "mc_pi_px[nLambda]/F");
    fRecTree->Branch("mc_pi_py", mc_pi_py, "mc_pi_py[nLambda]/F");
    fRecTree->Branch("mc_pi_pz", mc_pi_pz, "mc_pi_pz[nLambda]/F");

    fRecTree->Branch("proton_pdg", proton_pdg, "proton_pdg[nLambda]/I");
    fRecTree->Branch("proton_dedx", proton_dedx, "proton_dedx[nLambda]/F");
    fRecTree->Branch("proton_dca", proton_dca, "proton_dca[nLambda]/F");
    fRecTree->Branch("proton_px", proton_px, "proton_px[nLambda]/F");
    fRecTree->Branch("proton_py", proton_py, "proton_py[nLambda]/F");
    fRecTree->Branch("proton_pz", proton_pz, "proton_pz[nLambda]/F");
    fRecTree->Branch("mc_proton_px", mc_proton_px, "mc_proton_px[nLambda]/F");
    fRecTree->Branch("mc_proton_py", mc_proton_py, "mc_proton_py[nLambda]/F");
    fRecTree->Branch("mc_proton_pz", mc_proton_pz, "mc_proton_pz[nLambda]/F");
    fRecTree->Branch("mc_proton_dca", mc_proton_dca, "mc_proton_dca[nLambda]/F");

    fRecTree->Branch("proton_source", proton_source, "proton_source[nLambda]/I");
    fRecTree->Branch("pi_source", pi_source, "pi_source[nLambda]/I");
    fRecTree->Branch("has_mwdc", has_mwdc, "has_mwdc[nLambda]/I");
    fRecTree->Branch("proton_theta", proton_theta, "proton_theta[nLambda]/F");
    fRecTree->Branch("pi_theta", pi_theta, "pi_theta[nLambda]/F");

    // Existing MWDC summary branches
    fRecTree->Branch("mwdc_true_p_mc_dca", mwdc_true_p_mc_dca, "mwdc_true_p_mc_dca[nLambda]/F");
    fRecTree->Branch("mwdc_true_pi_mc_dca", mwdc_true_pi_mc_dca, "mwdc_true_pi_mc_dca[nLambda]/F");
    fRecTree->Branch("mwdc_p_dp", mwdc_p_dp, "mwdc_p_dp[nLambda]/F");
    fRecTree->Branch("mwdc_p_dpt", mwdc_p_dpt, "mwdc_p_dpt[nLambda]/F");
    fRecTree->Branch("mwdc_p_dtheta", mwdc_p_dtheta, "mwdc_p_dtheta[nLambda]/F");
    fRecTree->Branch("mwdc_p_dphi", mwdc_p_dphi, "mwdc_p_dphi[nLambda]/F");
    fRecTree->Branch("mwdc_p_dr", mwdc_p_dr, "mwdc_p_dr[nLambda]/F");
    fRecTree->Branch("mwdc_p_ddca", mwdc_p_ddca, "mwdc_p_ddca[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dp", mwdc_pi_dp, "mwdc_pi_dp[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dpt", mwdc_pi_dpt, "mwdc_pi_dpt[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dtheta", mwdc_pi_dtheta, "mwdc_pi_dtheta[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dphi", mwdc_pi_dphi, "mwdc_pi_dphi[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dr", mwdc_pi_dr, "mwdc_pi_dr[nLambda]/F");
    fRecTree->Branch("mwdc_pi_ddca", mwdc_pi_ddca, "mwdc_pi_ddca[nLambda]/F");

    // New diagnostic branches
    fRecTree->Branch("pid_mode", g_pid_mode, "pid_mode[nLambda]/I");
    fRecTree->Branch("pair_type", g_pair_type, "pair_type[nLambda]/I");
    fRecTree->Branch("same_pid_mcid", g_same_pid_mcid, "same_pid_mcid[nLambda]/I");
    fRecTree->Branch("p_pid_mcid", g_p_pid_mcid, "p_pid_mcid[nLambda]/I");
    fRecTree->Branch("pi_pid_mcid", g_pi_pid_mcid, "pi_pid_mcid[nLambda]/I");
    fRecTree->Branch("p_mother_id", g_p_mother_id, "p_mother_id[nLambda]/I");
    fRecTree->Branch("pi_mother_id", g_pi_mother_id, "pi_mother_id[nLambda]/I");
    fRecTree->Branch("p_mother_pdg", g_p_mother_pdg, "p_mother_pdg[nLambda]/I");
    fRecTree->Branch("pi_mother_pdg", g_pi_mother_pdg, "pi_mother_pdg[nLambda]/I");

    fRecTree->Branch("p_tof_det", g_p_tof_det, "p_tof_det[nLambda]/I");
    fRecTree->Branch("pi_tof_det", g_pi_tof_det, "pi_tof_det[nLambda]/I");
    // User-friendly aliases for ROOT cut scans: 1=iTOF, 2=eTOF.
    fRecTree->Branch("proton_tofdet", g_p_tof_det, "proton_tofdet[nLambda]/I");
    fRecTree->Branch("pion_tofdet", g_pi_tof_det, "pion_tofdet[nLambda]/I");
    fRecTree->Branch("p_tof_link", g_p_tof_link, "p_tof_link[nLambda]/I");
    fRecTree->Branch("pi_tof_link", g_pi_tof_link, "pi_tof_link[nLambda]/I");
    fRecTree->Branch("p_hit_pdg", g_p_hit_pdg, "p_hit_pdg[nLambda]/I");
    fRecTree->Branch("pi_hit_pdg", g_pi_hit_pdg, "pi_hit_pdg[nLambda]/I");
    fRecTree->Branch("p_hit_mother_id", g_p_hit_mother_id, "p_hit_mother_id[nLambda]/I");
    fRecTree->Branch("pi_hit_mother_id", g_pi_hit_mother_id, "pi_hit_mother_id[nLambda]/I");
    fRecTree->Branch("p_hit_mother_pdg", g_p_hit_mother_pdg, "p_hit_mother_pdg[nLambda]/I");
    fRecTree->Branch("pi_hit_mother_pdg", g_pi_hit_mother_pdg, "pi_hit_mother_pdg[nLambda]/I");
    fRecTree->Branch("p_track_hit_match", g_p_track_hit_match, "p_track_hit_match[nLambda]/I");
    fRecTree->Branch("pi_track_hit_match", g_pi_track_hit_match, "pi_track_hit_match[nLambda]/I");
    fRecTree->Branch("p_tof_correct", g_p_tof_correct, "p_tof_correct[nLambda]/I");
    fRecTree->Branch("pi_tof_correct", g_pi_tof_correct, "pi_tof_correct[nLambda]/I");
    fRecTree->Branch("p_tof_dca", g_p_tof_dca, "p_tof_dca[nLambda]/F");
    fRecTree->Branch("pi_tof_dca", g_pi_tof_dca, "pi_tof_dca[nLambda]/F");
    // User-friendly aliases for ROOT cut scans.
    fRecTree->Branch("proton_tofdca", g_p_tof_dca, "proton_tofdca[nLambda]/F");
    fRecTree->Branch("pion_tofdca", g_pi_tof_dca, "pion_tofdca[nLambda]/F");
    fRecTree->Branch("p_tof_dt", g_p_tof_dt, "p_tof_dt[nLambda]/F");
    fRecTree->Branch("pi_tof_dt", g_pi_tof_dt, "pi_tof_dt[nLambda]/F");

    fRecTree->Branch("mwdc_p_mcid", g_mwdc_p_mcid, "mwdc_p_mcid[nLambda]/I");
    fRecTree->Branch("mwdc_p_pdg", g_mwdc_p_pdg, "mwdc_p_pdg[nLambda]/I");
    fRecTree->Branch("mwdc_p_mother_pdg", g_mwdc_p_mother_pdg, "mwdc_p_mother_pdg[nLambda]/I");
    fRecTree->Branch("mwdc_p_purity", g_mwdc_p_purity, "mwdc_p_purity[nLambda]/F");
    fRecTree->Branch("mwdc_p_dx", g_mwdc_p_dx, "mwdc_p_dx[nLambda]/F");
    fRecTree->Branch("mwdc_p_dy", g_mwdc_p_dy, "mwdc_p_dy[nLambda]/F");
    fRecTree->Branch("mwdc_p_dz", g_mwdc_p_dz, "mwdc_p_dz[nLambda]/F");
    fRecTree->Branch("mwdc_p_dpx", g_mwdc_p_dpx, "mwdc_p_dpx[nLambda]/F");
    fRecTree->Branch("mwdc_p_dpy", g_mwdc_p_dpy, "mwdc_p_dpy[nLambda]/F");
    fRecTree->Branch("mwdc_p_dpz", g_mwdc_p_dpz, "mwdc_p_dpz[nLambda]/F");
    fRecTree->Branch("mwdc_p_dangle", g_mwdc_p_dangle, "mwdc_p_dangle[nLambda]/F");
    fRecTree->Branch("mwdc_p_line_dca_mcstart", g_mwdc_p_line_dca_mcstart, "mwdc_p_line_dca_mcstart[nLambda]/F");
    fRecTree->Branch("mwdc_p_dca_reco", g_mwdc_p_dca_reco, "mwdc_p_dca_reco[nLambda]/F");
    fRecTree->Branch("mwdc_p_dca_mc", g_mwdc_p_dca_mc, "mwdc_p_dca_mc[nLambda]/F");

    fRecTree->Branch("mwdc_pi_mcid", g_mwdc_pi_mcid, "mwdc_pi_mcid[nLambda]/I");
    fRecTree->Branch("mwdc_pi_pdg", g_mwdc_pi_pdg, "mwdc_pi_pdg[nLambda]/I");
    fRecTree->Branch("mwdc_pi_mother_pdg", g_mwdc_pi_mother_pdg, "mwdc_pi_mother_pdg[nLambda]/I");
    fRecTree->Branch("mwdc_pi_purity", g_mwdc_pi_purity, "mwdc_pi_purity[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dx", g_mwdc_pi_dx, "mwdc_pi_dx[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dy", g_mwdc_pi_dy, "mwdc_pi_dy[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dz", g_mwdc_pi_dz, "mwdc_pi_dz[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dpx", g_mwdc_pi_dpx, "mwdc_pi_dpx[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dpy", g_mwdc_pi_dpy, "mwdc_pi_dpy[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dpz", g_mwdc_pi_dpz, "mwdc_pi_dpz[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dangle", g_mwdc_pi_dangle, "mwdc_pi_dangle[nLambda]/F");
    fRecTree->Branch("mwdc_pi_line_dca_mcstart", g_mwdc_pi_line_dca_mcstart, "mwdc_pi_line_dca_mcstart[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dca_reco", g_mwdc_pi_dca_reco, "mwdc_pi_dca_reco[nLambda]/F");
    fRecTree->Branch("mwdc_pi_dca_mc", g_mwdc_pi_dca_mc, "mwdc_pi_dca_mc[nLambda]/F");

    mc_pty = new TH2F("mc_pty", ";y;p_{T}(GeV/c)", 50, -1.4, 1.4, 50, 0, 1.6);
    mc_pt = new TH1F("mc_pt", ";p_{T}(GeV/c);counts", 20, 0, 1.4);
    mc_v1 = new TProfile("mc_v1", ";y;v_{1}", 6, -1, 1);
    rec_v1 = new TProfile("rec_v1", ";y;v_{1}", 6, -1, 1);

    itofTime_diff = new TH1F("itofTime_diff", ";#Delta t;counts", 100, -0.3, 0.3);
    itofLength_diff = new TH1F("itofLength_diff", ";#Delta L;counts", 100, -10, 20);
    tofdca = new TH1F("tofdca", ";track dca to tof;counts", 100, 0, 40);
    mc_tofdca = new TH1F("mc_tofdca", ";matched MC track dca to tof;counts", 100, 0, 40);

    if (fDebug > 0) Info("InitialHistgram", "MWDC diagnostic histograms and tree are defined.");
}

void AnaMaker::Finish()
{
    mc_pty->Write();
    mc_pt->Write();
    mc_v1->Write();
    rec_v1->Write();
    itofTime_diff->Write();
    itofLength_diff->Write();
    tofdca->Write();
    mc_tofdca->Write();
    fRecTree->Write();
    if (gMwdcTrackTree) {gMwdcTrackTree->Write();}
    if (gTpcMwdcMatchTree) {gTpcMwdcMatchTree->Write();}
    
    fOutFile->Close();

    if (fDebug > 0) Info("Finish", "Output saved.");
}

void AnaMaker::ProcessData()
{
    InitialHistgram();

    if (fDebug > 0)
        Info("ProcessData", "%d out of total entries %lld will be processed", fEvtNo, fChain->GetEntries());

    gRandom->SetSeed(12345);

    for (int ievent = 0; ievent < fEvtNo; ++ievent) {
        evtId = ievent;
        nLambda = 0;
        psi = 0;
        InitializeArrays();

        if (ievent % 1000 == 0) std::cout << "\r data process:" << ievent << std::flush;

        fChain->GetEntry(ievent);

        // ========================================================
        // TPC-MWDC matching diagnostic at z = 80 cm.
        //
        // TPC:
        //   position  = reconstructed CeeTpcTrack::GetX/GetY/GetZ
        //   direction = reconstructed CeeTpcTrack momentum
        //
        // MWDC:
        //   position  = reconstructed r1 (mm -> cm)
        //   direction = reconstructed p
        //
        // Every valid TPC-MWDC pair is stored. same_mcid is only
        // a truth label for evaluating the cut; it is not used to
        // select or reject a pair.
        // ========================================================
        if (gTpcMwdcMatchTree &&
            fTpcTrack &&
            fMWDC_Track &&
            fMCTrack) {

            for (Int_t itpc = 0;
                 itpc < fTpcTrack->GetEntriesFast();
                 ++itpc) {

                CeeTpcTrack* tpcTrack =
                    dynamic_cast<CeeTpcTrack*>(fTpcTrack->At(itpc));

                if (!tpcTrack) continue;

                const Int_t tpcMcid = tpcTrack->GetMCTrackID();

                if (tpcMcid < 0 ||
                    tpcMcid >= fMCTrack->GetEntriesFast()) {
                    continue;
                }

                CeeMCTrack* tpcMCTrack =
                    dynamic_cast<CeeMCTrack*>(fMCTrack->At(tpcMcid));

                if (!tpcMCTrack) continue;

                const TVector3 tpcPos(
                    tpcTrack->GetX(),
                    tpcTrack->GetY(),
                    tpcTrack->GetZ()
                );

                const TVector3 tpcMom(
                    tpcTrack->GetPx(),
                    tpcTrack->GetPy(),
                    tpcTrack->GetPz()
                );

                if (!std::isfinite(tpcPos.X()) ||
                    !std::isfinite(tpcPos.Y()) ||
                    !std::isfinite(tpcPos.Z())) {
                    continue;
                }

                TVector3 tpcProj;

                if (!ProjectTrackToZPlane(
                        tpcPos,
                        tpcMom,
                        kTpcMwdcMatchZ,
                        tpcProj)) {
                    continue;
                }

                for (Int_t imwdc = 0;
                     imwdc < fMWDC_Track->GetEntriesFast();
                     ++imwdc) {

                    CeeMWDC_Track* mwdcTrack =
                        dynamic_cast<CeeMWDC_Track*>(
                            fMWDC_Track->At(imwdc)
                        );

                    if (!mwdcTrack) continue;

                    const Int_t mwdcMcid =
                        mwdcTrack->track_mcid;

                    // A valid MC id is required only so that this
                    // pair can be labelled as true or false.
                    if (mwdcMcid < 0 ||
                        mwdcMcid >= fMCTrack->GetEntriesFast()) {
                        continue;
                    }

                    const TVector3 mwdcPos =
                        mwdcTrack->r1 * 0.1; // mm -> cm

                    const TVector3 mwdcMom =
                        mwdcTrack->p;

                    TVector3 mwdcProj;

                    if (!ProjectTrackToZPlane(
                            mwdcPos,
                            mwdcMom,
                            kTpcMwdcMatchZ,
                            mwdcProj)) {
                        continue;
                    }

                    const Double_t dx =
                        tpcProj.X() - mwdcProj.X();

                    const Double_t dy =
                        tpcProj.Y() - mwdcProj.Y();

                    g_match_evt_id = ievent;
                    g_match_tpc_index = itpc;
                    g_match_mwdc_index = imwdc;
                    g_match_tpc_mcid = tpcMcid;
                    g_match_mwdc_mcid = mwdcMcid;
                    g_match_same_mcid =
                        (tpcMcid == mwdcMcid) ? 1 : 0;
                    g_match_dx_z80 =
                        static_cast<Float_t>(dx);
                    g_match_dy_z80 =
                        static_cast<Float_t>(dy);
                    g_match_dr_z80 =
                        static_cast<Float_t>(
                            std::sqrt(dx * dx + dy * dy)
                        );

                    gTpcMwdcMatchTree->Fill();
                }
            }
        }

        // Check mwdc resolution
                if (gMwdcTrackTree && fMWDC_Track && fMCTrack) {

            for (Int_t imwdc = 0;
                 imwdc < fMWDC_Track->GetEntriesFast();
                 ++imwdc) {

                CeeMWDC_Track* mwdcTrack =
                    dynamic_cast<CeeMWDC_Track*>(
                        fMWDC_Track->At(imwdc)
                    );

                if (!mwdcTrack) continue;

                const Int_t mcid = mwdcTrack->track_mcid;

                // No valid MC truth link
                if (mcid < 0 ||
                    mcid >= fMCTrack->GetEntriesFast()) {
                    continue;
                }

                Int_t pdg       = -999999;
                Int_t motherId  = -999999;
                Int_t motherPdg = -999999;

                TVector3 mcPos(0., 0., 0.);
                TVector3 mcMom(0., 0., 0.);

                Bool_t hasMCInfo = GetMCInfoFromMCID(
                    fMCTrack,
                    mcid,
                    pdg,
                    motherId,
                    motherPdg,
                    mcPos,
                    mcMom
                );

                if (!hasMCInfo) continue;

                // Reconstructed momentum from MWDC track
                TVector3 recoMom = mwdcTrack->p;

                // Reject zero or invalid momentum
                if (recoMom.Mag() < 1.0e-12) continue;
                if (mcMom.Mag()   < 1.0e-12) continue;

                // ------------------------------------------------
                // Identity information
                // ------------------------------------------------
                g_mwdc_evt_id      = ievent;
                g_mwdc_track_index = imwdc;
                g_mwdc_track_mcid  = mcid;

                g_mwdc_track_pdg        = pdg;
                g_mwdc_track_mother_id  = motherId;
                g_mwdc_track_mother_pdg = motherPdg;
                g_mwdc_track_is_primary =
                    (motherId == -1) ? 1 : 0;

                // ------------------------------------------------
                // Track reconstruction quality
                // ------------------------------------------------
                g_mwdc_track_hit_num =
                    mwdcTrack->hit_num;

                g_mwdc_track_num_true_hit =
                    mwdcTrack->num_mchit_true;

                if (mwdcTrack->hit_num > 0) {
                    g_mwdc_track_purity =
                        static_cast<Float_t>(
                            mwdcTrack->num_mchit_true
                        ) /
                        static_cast<Float_t>(
                            mwdcTrack->hit_num
                        );
                }
                else {
                    g_mwdc_track_purity = -1.0;
                }

                // ------------------------------------------------
                // Reconstructed and MC momentum
                // ------------------------------------------------
                g_mwdc_track_p_reco =
                    recoMom.Mag();

                g_mwdc_track_p_mc =
                    mcMom.Mag();

                g_mwdc_track_theta_reco =
                    recoMom.Theta() * TMath::RadToDeg();

                g_mwdc_track_theta_mc =
                    mcMom.Theta() * TMath::RadToDeg();
                g_mwdc_track_angle_res =
                    recoMom.Angle(mcMom) * TMath::RadToDeg();
                // ------------------------------------------------
                // Angular resolution:
                // 3D opening angle between reco and MC momentum
                // Unit: degree
                // ------------------------------------------------
                g_mwdc_track_angle_res =
                    recoMom.Angle(mcMom) *
                    TMath::RadToDeg();

                // ------------------------------------------------
                // Relative momentum resolution:
                // (p_reco - p_MC) / p_MC
                // ------------------------------------------------
                g_mwdc_track_mom_res =
                    (recoMom.Mag() - mcMom.Mag()) /
                    mcMom.Mag();

                // One entry for one MWDC reconstructed track
                gMwdcTrackTree->Fill();
            }
        }

        
        
        // MC Lambda info
        for (Int_t i1 = 0; i1 < fMCTrack->GetEntries(); i1++) {
            iMCTrack = (CeeMCTrack*)fMCTrack->At(i1);
            int pdg = iMCTrack->GetPdgCode();
            int motherID = iMCTrack->GetMotherId();
            if (abs(pdg) != 3122 || motherID != -1) continue;
            double pt = iMCTrack->GetPt();
            double y = iMCTrack->GetRapidity();
            y = y - 0.7;
            mc_pty->Fill(y, pt);
            mc_pt->Fill(pt);
            TVector3 mom;
            iMCTrack->GetMomentum(mom);
            mc_v1->Fill(y, cos(mom.Phi() - psi));
        }

        // This version keeps TPC-TPC pairs as the baseline, so events are not
        // rejected just because there is no MWDC track. MWDC-related pairs are
        // naturally absent when MWDC_Track is empty.
        if (!fTpcTrack || fTpcTrack->GetEntries() < 1) continue;
        if ((fiTOFHit->GetEntries() + fEtofHit->GetEntries()) < 2) continue;

        std::vector<RecoTrackCand> recoTrackCands;
        BuildRecoTrackCandList(recoTrackCands, fTpcTrack, fMWDC_Track, fMCTrack);
        if (recoTrackCands.size() < 2) continue;

        // Event vertex
        bool hasRecoVtx = false;
        fEvtVtx.SetXYZ(0, 0, -35);

        if (kVertexMode == 0) {
            for (int itpcvtx = 0; itpcvtx < fTpcVertex->GetEntries(); itpcvtx++) {
                CeeTpcVertex* iTpcVertex = dynamic_cast<CeeTpcVertex*>(fTpcVertex->At(itpcvtx));
                if (!iTpcVertex) continue;
                fEvtVtx.SetXYZ(iTpcVertex->GetX()[0], iTpcVertex->GetY()[0], iTpcVertex->GetZ()[0]);
                hasRecoVtx = true;
                break;
            }
            if (!hasRecoVtx) continue;
        }
        else if (kVertexMode == 1) {
            fEvtVtx.SetXYZ(0, 0, -35);
        }
        else if (kVertexMode == 2) {
            const double sigmaX = 0.093;
            const double sigmaY = 0.081;
            const double sigmaZ = 0.112;
            fEvtVtx.SetXYZ(gRandom->Gaus(0.0, sigmaX),
                           gRandom->Gaus(0.0, sigmaY),
                           gRandom->Gaus(-35., sigmaZ));
        }

        if (fDebug > 0) {
            Info("ProcessData",
                 "Event vertex: x=%1.2f y=%1.2f z=%1.2f, TPC tracks=%d, MWDC tracks=%d, selected=%d, iTOF hits=%d, eTOF hits=%d",
                 fEvtVtx.X(), fEvtVtx.Y(), fEvtVtx.Z(),
                 fTpcTrack ? fTpcTrack->GetEntries() : 0,
                 fMWDC_Track ? fMWDC_Track->GetEntries() : 0,
                 (int)recoTrackCands.size(),
                 fiTOFHit ? fiTOFHit->GetEntries() : 0,
                 fEtofHit ? fEtofHit->GetEntries() : 0);
        }

        // global iTOF hit time/length diagnostic
        for (int ihit = 0; ihit < fiTOFHit->GetEntries(); ihit++) {
            CeeiTOFHit* iitofHit = dynamic_cast<CeeiTOFHit*>(fiTOFHit->At(ihit));
            if (!iitofHit) continue;
            TVector3 length = TVector3(iitofHit->GetX(), iitofHit->GetY(), iitofHit->GetZ() + 35);
            Double_t time = iitofHit->GetTime();
            FairLink tmplink2 = iitofHit->GetLink(1);
            Int_t MC2_id = tmplink2.GetIndex();
            CeeiTOFPoint* iitofPoint = dynamic_cast<CeeiTOFPoint*>(fiTOFPoint->At(MC2_id));
            if (!iitofPoint) continue;
            itofTime_diff->Fill(time - iitofPoint->GetTime());
            itofLength_diff->Fill(length.Mag() - iitofPoint->GetLength());
        }

        // proton loop
        for (size_t itrack = 0; itrack < recoTrackCands.size(); itrack++) {
            const RecoTrackCand& pCand = recoTrackCands[itrack];
            if (pCand.mcId < 0) continue;

            TVector3 pos_p = pCand.pos;
            TVector3 mom_p = pCand.mom;
            Int_t MC_id1 = pCand.mcId;

            // Track PDG PID: the reconstructed TPC/MWDC track itself must be proton.
            Int_t mcPdg1 = pCand.mcPdg;
            Int_t motherId1 = pCand.motherId;
            Int_t motherPdg1 = pCand.motherPdg;
            TVector3 mc_pos_p = pCand.mcPos;
            TVector3 mc_mom_p = pCand.mcMom;
            
            //if (mcPdg1 != 2212) continue;
            // DCA-based proton assignment.
            // Smaller daughter DCA is treated as proton.
            // TPC uses 1.2 cm, MWDC uses DaughterDcaCut_MWDC.
            Double_t dca1 = GetDca(pos_p, mom_p, fEvtVtx);
            Double_t dcaCut1 = GetDaughterDcaCutBySource(pCand.source);
            if (dca1 > dcaCut1) continue;
            

            TofMatchResult tofMatch1 = SelectTofPreferItof(
                FindNearestItofHit(pos_p, mom_p, MC_id1, fiTOFHit, fiTOFPoint),
                FindNearestEtofHit(pos_p, mom_p, MC_id1, fEtofHit, feTOFPoint),
                pCand.source);

            if (!tofMatch1.found) continue;
            Double_t tofMatchCut1 = GetTofMatchCutBySourceDet(pCand.source, tofMatch1.detType);
            if (fabs(tofMatch1.dca) > tofMatchCut1) continue;

            // Hit link PDG: the selected hit must be linked to a proton MC track.
            Int_t pHitPdg = -999999;
            Int_t pHitMotherId = -999999;
            Int_t pHitMotherPdg = -999999;
            TVector3 pHitMcPos(0., 0., 0.);
            TVector3 pHitMcMom(0., 0., 0.);
            Bool_t pHitInfoOK = GetMCInfoFromMCID(fMCTrack, tofMatch1.linkTrack,
                                                pHitPdg, pHitMotherId, pHitMotherPdg,
                                                pHitMcPos, pHitMcMom);
            if (kRequireHitPDG) {
                if (!pHitInfoOK) continue;
                if (pHitPdg != 2212) continue;
            }

            // Track-hit truth match: selected hit must belong to this track's MC id.
            Bool_t pTrackHitMatch = (tofMatch1.linkTrack == MC_id1);
            if (kRequireTrackHitLink && !pTrackHitMatch) continue;
            if ((kUseMCTofPointPos || kUseMCTofPointTime) && !tofMatch1.hasPointInfo) continue;

            tofdca->Fill(tofMatch1.dca);
            mc_tofdca->Fill(tofMatch1.mcDca);

           // Double_t dca1 = GetDca(pos_p, mom_p, fEvtVtx);

            TVector3 tofPos1 = (kUseMCTofPointPos ? tofMatch1.mcHitPos : tofMatch1.hitPos);
            Double_t tofTime1 = (kUseMCTofPointTime ? tofMatch1.mcTime : tofMatch1.time);

            MyTrack* myproton = new MyTrack(pos_p, tofPos1, fEvtVtx, tofTime1, pmass, mom_p);
            if (!myproton || myproton->mom().Mag() < 0.005) {
                delete myproton;
                continue;
            }

            if (fRotP) {
                gRandom->SetSeed(0);
                Double_t angle = gRandom->Uniform(150, 210) / 180.0 * TMath::Pi();
                myproton->RotateXY(angle);
            }

            // pion loop
            for (size_t jtrack = 0; jtrack < recoTrackCands.size(); jtrack++) {
                if (itrack == jtrack) continue;
                const RecoTrackCand& piCand = recoTrackCands[jtrack];
                if (piCand.mcId < 0) continue;

                const Int_t pairType = GetPairType(pCand.source, piCand.source);

                // Keep only two combinations:
                //   pair_type == 0 : TPC proton  + TPC pion
                //   pair_type == 1 : MWDC proton + TPC pion
                // Reject all combinations with MWDC pion:
                //   pair_type == 2 : TPC proton  + MWDC pion
                //   pair_type == 3 : MWDC proton + MWDC pion
                if (kKeepOnlyTpcTpcAndMwdcPTpcPi &&
                    pairType != 0 &&
                    pairType != 1) continue;

                TVector3 pos_pi = piCand.pos;
                TVector3 mom_pi = piCand.mom;
                Int_t MC_id2 = piCand.mcId;

                // Track PDG PID: the reconstructed TPC/MWDC track itself must be pi-.
                Int_t mcPdg2 = piCand.mcPdg;
                Int_t motherId2 = piCand.motherId;
                Int_t motherPdg2 = piCand.motherPdg;
                TVector3 mc_pos_pi = piCand.mcPos;
                TVector3 mc_mom_pi = piCand.mcMom;
                //if (mcPdg2 != -211) continue;
                // DCA-based pion assignment.
                // Larger daughter DCA is treated as pion.
                Double_t dca2 = GetDca(pos_pi, mom_pi, fEvtVtx);
                Double_t dcaCut2 = GetDaughterDcaCutBySource(piCand.source);
                if (dca2 <= dcaCut2) continue; 
                
                // Pion iTOF/eTOF selection:
                // if the pion selects iTOF, apply a loose 5 cm matching cut;
                // if it falls back to eTOF, keep the source-dependent eTOF matching cut.
                TofMatchResult itofMatch2 =
                    FindNearestItofHit(pos_pi, mom_pi, MC_id2, fiTOFHit, fiTOFPoint);
                TofMatchResult etofMatch2 =
                    FindNearestEtofHit(pos_pi, mom_pi, MC_id2, fEtofHit, feTOFPoint);
                TofMatchResult tofMatch2 =
                    SelectPionTofPreferItof(itofMatch2, etofMatch2);

                if (!tofMatch2.found) continue;
                if (tofMatch2.hitPos.Mag() < 10) continue;

                Double_t tofMatchCut2 = 1.0e9;
                if (tofMatch2.detType == 1) {
                    tofMatchCut2 = PionItofMatchCut;
                }
                else if (tofMatch2.detType == 2) {
                    tofMatchCut2 = GetTofMatchCutBySourceDet(piCand.source, 2);
                }
                if (fabs(tofMatch2.dca) > tofMatchCut2) continue;

                // Hit link PDG: the selected hit must be linked to a pi- MC track.
                Int_t piHitPdg = -999999;
                Int_t piHitMotherId = -999999;
                Int_t piHitMotherPdg = -999999;
                TVector3 piHitMcPos(0., 0., 0.);
                TVector3 piHitMcMom(0., 0., 0.);
                Bool_t piHitInfoOK = GetMCInfoFromMCID(fMCTrack, tofMatch2.linkTrack,
                                                  piHitPdg, piHitMotherId, piHitMotherPdg,
                                                  piHitMcPos, piHitMcMom);
                if (kRequireHitPDG) {
                    if (!piHitInfoOK) continue;
                    if (piHitPdg != -211) continue;
                }

                // Track-hit truth match: selected hit must belong to this track's MC id.
                Bool_t piTrackHitMatch = (tofMatch2.linkTrack == MC_id2);
                if (kRequireTrackHitLink && !piTrackHitMatch) continue;
                if ((kUseMCTofPointPos || kUseMCTofPointTime) && !tofMatch2.hasPointInfo) continue;

               // Double_t dca2 = GetDca(pos_pi, mom_pi, fEvtVtx);

                TVector3 tofPos2 = (kUseMCTofPointPos ? tofMatch2.mcHitPos : tofMatch2.hitPos);
                Double_t tofTime2 = (kUseMCTofPointTime ? tofMatch2.mcTime : tofMatch2.time);

                MyTrack* mypion = new MyTrack(pos_pi, tofPos2, fEvtVtx, tofTime2, pimass, mom_pi);
                if (!mypion || mypion->mom().Mag() < 0.001) {
                    delete mypion;
                    continue;
                }

                MyLambda* myld = new MyLambda(myproton, mypion, fRotPi);
                if (!myld) {
                    delete mypion;
                    continue;
                }

                TVector3 LdMom = myld->mom();
                TVector3 LdDecayVtx = myld->decayVertex();
                Double_t LdDca = GetDca(LdDecayVtx, LdMom, fEvtVtx);
                Double_t L = (LdDecayVtx - fEvtVtx).Mag();
                Double_t Dca1to2 = myld->dca1to2();

                TLorentzVector Ldp4(myproton->px() + mypion->px(),
                                    myproton->py() + mypion->py(),
                                    myproton->pz() + mypion->pz(),
                                    myproton->e()  + mypion->e());
                double tmp_mass = Ldp4.M();
                if (tmp_mass < 0.7 || tmp_mass > 2.0) {
                    delete myld;
                    delete mypion;
                    continue;
                }

                if (nLambda >= kMaxLambdaDiag) {
                    delete myld;
                    delete mypion;
                    continue;
                }

                ld_mass[nLambda] = tmp_mass;
                proton_source[nLambda] = pCand.source;
                pi_source[nLambda] = piCand.source;
                has_mwdc[nLambda] = (pairType > 0) ? 1 : 0;
                proton_theta[nLambda] = mom_p.Theta() * TMath::RadToDeg();
                pi_theta[nLambda] = mom_pi.Theta() * TMath::RadToDeg();

                Double_t pDcaMC = -999.;
                if (pCand.mcMom.Mag() > 1.0e-12) pDcaMC = GetDca(pCand.mcPos, pCand.mcMom, fEvtVtx);
                Double_t piDcaMC = -999.;
                if (piCand.mcMom.Mag() > 1.0e-12) piDcaMC = GetDca(piCand.mcPos, piCand.mcMom, fEvtVtx);

                mwdc_true_p_mc_dca[nLambda] = (pCand.source == 1) ? pDcaMC : -999.;
                mwdc_true_pi_mc_dca[nLambda] = (piCand.source == 1) ? piDcaMC : -999.;

                if (pCand.source == 1 && pCand.mcMom.Mag() > 1e-12) {
                    mwdc_p_dp[nLambda] = mom_p.Mag() - pCand.mcMom.Mag();
                    mwdc_p_dpt[nLambda] = mom_p.Pt() - pCand.mcMom.Pt();
                    mwdc_p_dtheta[nLambda] = mom_p.Theta() * TMath::RadToDeg()
                                           - pCand.mcMom.Theta() * TMath::RadToDeg();
                    mwdc_p_dphi[nLambda] = NormalizeDeltaPhiDeg(mom_p.Phi() * TMath::RadToDeg()
                                                               - pCand.mcMom.Phi() * TMath::RadToDeg());
                    mwdc_p_dr[nLambda] = (pos_p - pCand.mcPos).Mag();
                    mwdc_p_ddca[nLambda] = dca1 - pDcaMC;
                    FillMwdcDiagProton(nLambda, pCand, pos_p, mom_p, dca1, pDcaMC);
                }

                if (piCand.source == 1 && piCand.mcMom.Mag() > 1e-12) {
                    mwdc_pi_dp[nLambda] = mom_pi.Mag() - piCand.mcMom.Mag();
                    mwdc_pi_dpt[nLambda] = mom_pi.Pt() - piCand.mcMom.Pt();
                    mwdc_pi_dtheta[nLambda] = mom_pi.Theta() * TMath::RadToDeg()
                                            - piCand.mcMom.Theta() * TMath::RadToDeg();
                    mwdc_pi_dphi[nLambda] = NormalizeDeltaPhiDeg(mom_pi.Phi() * TMath::RadToDeg()
                                                                - piCand.mcMom.Phi() * TMath::RadToDeg());
                    mwdc_pi_dr[nLambda] = (pos_pi - piCand.mcPos).Mag();
                    mwdc_pi_ddca[nLambda] = dca2 - piDcaMC;
                    FillMwdcDiagPion(nLambda, piCand, pos_pi, mom_pi, dca2, piDcaMC);
                }

                ld_px[nLambda] = LdMom.X();
                ld_py[nLambda] = LdMom.Y();
                ld_pz[nLambda] = LdMom.Z();
                ld_dca[nLambda] = LdDca;
                decaylength[nLambda] = L;
                dca1to2[nLambda] = Dca1to2;

                int tmp_truth = 0;
                TVector3 mc_ldmom(0., 0., 0.);
                if (mcPdg1 == 2212 && mcPdg2 == -211 &&
                    motherPdg1 == 3122 && motherPdg2 == 3122 &&
                    motherId1 >= 0 && motherId1 == motherId2 &&
                    motherId1 < fMCTrack->GetEntriesFast()) {
                    CeeMCTrack* mother = dynamic_cast<CeeMCTrack*>(fMCTrack->At(motherId1));
                    if (mother && mother->GetPdgCode() == 3122 && mother->GetMotherId() == -1) {
                        tmp_truth = 1;
                        mother->GetMomentum(mc_ldmom);
                    }
                }
                Int_t tmp_truth_hit = 0;
                if (tmp_truth == 1 && pTrackHitMatch && piTrackHitMatch) {
                    tmp_truth_hit = 1;
                }
                ld_mctruth[nLambda] = tmp_truth;
                g_mctruth_hit[nLambda] = tmp_truth_hit;

                pi_pdg[nLambda] = mcPdg2;
                pi_dca[nLambda] = dca2;
                pi_dedx[nLambda] = dca2;
                pi_px[nLambda] = mypion->mom().X();
                pi_py[nLambda] = mypion->mom().Y();
                pi_pz[nLambda] = mypion->mom().Z();
                mc_pi_px[nLambda] = mc_mom_pi.X();
                mc_pi_py[nLambda] = mc_mom_pi.Y();
                mc_pi_pz[nLambda] = mc_mom_pi.Z();

                proton_pdg[nLambda] = mcPdg1;
                proton_dca[nLambda] = dca1;
                proton_dedx[nLambda] = dca1;
                proton_px[nLambda] = myproton->mom().X();
                proton_py[nLambda] = myproton->mom().Y();
                proton_pz[nLambda] = myproton->mom().Z();
                mc_proton_px[nLambda] = mc_mom_p.X();
                mc_proton_py[nLambda] = mc_mom_p.Y();
                mc_proton_pz[nLambda] = mc_mom_p.Z();
                mc_proton_dca[nLambda] = pDcaMC;

                mc_ld_px[nLambda] = mc_ldmom.X();
                mc_ld_py[nLambda] = mc_ldmom.Y();
                mc_ld_pz[nLambda] = mc_ldmom.Z();

                g_pid_mode[nLambda] = kPidMode;
                g_pair_type[nLambda] = pairType;
                g_same_pid_mcid[nLambda] = (MC_id1 >= 0 && MC_id1 == MC_id2) ? 1 : 0;
                g_p_pid_mcid[nLambda] = MC_id1;
                g_pi_pid_mcid[nLambda] = MC_id2;
                g_p_mother_id[nLambda] = motherId1;
                g_pi_mother_id[nLambda] = motherId2;
                g_p_mother_pdg[nLambda] = motherPdg1;
                g_pi_mother_pdg[nLambda] = motherPdg2;

                g_p_tof_det[nLambda] = tofMatch1.detType;
                g_pi_tof_det[nLambda] = tofMatch2.detType;
                g_p_tof_link[nLambda] = tofMatch1.linkTrack;
                g_pi_tof_link[nLambda] = tofMatch2.linkTrack;
                g_p_hit_pdg[nLambda] = pHitPdg;
                g_pi_hit_pdg[nLambda] = piHitPdg;
                g_p_hit_mother_id[nLambda] = pHitMotherId;
                g_pi_hit_mother_id[nLambda] = piHitMotherId;
                g_p_hit_mother_pdg[nLambda] = pHitMotherPdg;
                g_pi_hit_mother_pdg[nLambda] = piHitMotherPdg;
                g_p_track_hit_match[nLambda] = pTrackHitMatch ? 1 : 0;
                g_pi_track_hit_match[nLambda] = piTrackHitMatch ? 1 : 0;
                g_p_tof_correct[nLambda] = pTrackHitMatch ? 1 : 0;
                g_pi_tof_correct[nLambda] = piTrackHitMatch ? 1 : 0;
                g_p_tof_dca[nLambda] = tofMatch1.dca;
                g_pi_tof_dca[nLambda] = tofMatch2.dca;
                g_p_tof_dt[nLambda] = tofMatch1.time - tofMatch1.mcTime;
                g_pi_tof_dt[nLambda] = tofMatch2.time - tofMatch2.mcTime;

                double y = Ldp4.Rapidity() - 0.7;
                if (tmp_truth == 1) rec_v1->Fill(y, cos(LdMom.Phi() - psi));

                nLambda++;

                delete myld;
                delete mypion;
            }

            delete myproton;
        }

        fRecTree->Fill();
    }

    Finish();
}

Double_t AnaMaker::GetDca(const TVector3& pos, const TVector3& p, const TVector3& vertex)
{
    if (p.Mag() < 1.0e-12) return 1.0e9;
    TVector3 diff = pos - vertex;
    TVector3 unitMom = p.Unit();
    Double_t length = diff.Mag();
    Double_t projLength = diff.Dot(unitMom);
    Double_t d2 = length * length - projLength * projLength;
    if (d2 <= 0) return 0.;
    return sqrt(d2);
}

Double_t AnaMaker::GetDcaXY(const TVector3& pos, const TVector3& p, const TVector3& vertex)
{
    Double_t dx = pos.X() - vertex.X();
    Double_t dy = pos.Y() - vertex.Y();
    Double_t px = p.X();
    Double_t py = p.Y();
    Double_t pMag2 = px * px + py * py;
    if (pMag2 < 1.0e-12) return 1.0e9;
    Double_t cross = fabs(dx * py - dy * px);
    return cross / sqrt(pMag2);
}

void AnaMaker::InitializeArrays()
{
    int size = 500;
    for (int i = 0; i < size; ++i) {
        pi_pdg[i] = 0;
        pi_dca[i] = 0.0;
        pi_px[i] = 0.0;
        pi_py[i] = 0.0;
        pi_pz[i] = 0.0;
        mc_pi_px[i] = 0.0;
        mc_pi_py[i] = 0.0;
        mc_pi_pz[i] = 0.0;

        proton_pdg[i] = 0;
        proton_dca[i] = 0.0;
        proton_px[i] = 0.0;
        proton_py[i] = 0.0;
        proton_pz[i] = 0.0;
        mc_proton_px[i] = 0.0;
        mc_proton_py[i] = 0.0;
        mc_proton_pz[i] = 0.0;
        mc_proton_dca[i] = 0.0;

        ld_mass[i] = 0.0;
        ld_mctruth[i] = 0;
        g_mctruth_hit[i] = 0;
        ld_dca[i] = 0.0;
        dca1to2[i] = 0.0;
        decaylength[i] = 0.0;
        ld_px[i] = 0.0;
        ld_py[i] = 0.0;
        ld_pz[i] = 0.0;
        mc_ld_px[i] = 0.0;
        mc_ld_py[i] = 0.0;
        mc_ld_pz[i] = 0.0;

        proton_source[i] = -1;
        pi_source[i] = -1;
        has_mwdc[i] = 0;
        proton_theta[i] = -999.;
        pi_theta[i] = -999.;

        mwdc_true_p_mc_dca[i] = -999.;
        mwdc_true_pi_mc_dca[i] = -999.;
        mwdc_p_dp[i] = -999.;
        mwdc_p_dpt[i] = -999.;
        mwdc_p_dtheta[i] = -999.;
        mwdc_p_dphi[i] = -999.;
        mwdc_p_dr[i] = -999.;
        mwdc_p_ddca[i] = -999.;
        mwdc_pi_dp[i] = -999.;
        mwdc_pi_dpt[i] = -999.;
        mwdc_pi_dtheta[i] = -999.;
        mwdc_pi_dphi[i] = -999.;
        mwdc_pi_dr[i] = -999.;
        mwdc_pi_ddca[i] = -999.;

        g_pid_mode[i] = -1;
        g_pair_type[i] = -1;
        g_same_pid_mcid[i] = -1;
        g_p_pid_mcid[i] = -999;
        g_pi_pid_mcid[i] = -999;
        g_p_mother_id[i] = -999;
        g_pi_mother_id[i] = -999;
        g_p_mother_pdg[i] = -999;
        g_pi_mother_pdg[i] = -999;

        g_p_tof_det[i] = -1;
        g_pi_tof_det[i] = -1;
        g_p_tof_link[i] = -999;
        g_pi_tof_link[i] = -999;
        g_p_hit_pdg[i] = -999999;
        g_pi_hit_pdg[i] = -999999;
        g_p_hit_mother_id[i] = -999;
        g_pi_hit_mother_id[i] = -999;
        g_p_hit_mother_pdg[i] = -999999;
        g_pi_hit_mother_pdg[i] = -999999;
        g_p_track_hit_match[i] = -1;
        g_pi_track_hit_match[i] = -1;
        g_p_tof_correct[i] = -1;
        g_pi_tof_correct[i] = -1;
        g_p_tof_dca[i] = -999.;
        g_pi_tof_dca[i] = -999.;
        g_p_tof_dt[i] = -999.;
        g_pi_tof_dt[i] = -999.;

        g_mwdc_p_mcid[i] = -999;
        g_mwdc_p_pdg[i] = -999999;
        g_mwdc_p_mother_pdg[i] = -999999;
        g_mwdc_p_purity[i] = -999.;
        g_mwdc_p_dx[i] = -999.;
        g_mwdc_p_dy[i] = -999.;
        g_mwdc_p_dz[i] = -999.;
        g_mwdc_p_dpx[i] = -999.;
        g_mwdc_p_dpy[i] = -999.;
        g_mwdc_p_dpz[i] = -999.;
        g_mwdc_p_dangle[i] = -999.;
        g_mwdc_p_line_dca_mcstart[i] = -999.;
        g_mwdc_p_dca_reco[i] = -999.;
        g_mwdc_p_dca_mc[i] = -999.;

        g_mwdc_pi_mcid[i] = -999;
        g_mwdc_pi_pdg[i] = -999999;
        g_mwdc_pi_mother_pdg[i] = -999999;
        g_mwdc_pi_purity[i] = -999.;
        g_mwdc_pi_dx[i] = -999.;
        g_mwdc_pi_dy[i] = -999.;
        g_mwdc_pi_dz[i] = -999.;
        g_mwdc_pi_dpx[i] = -999.;
        g_mwdc_pi_dpy[i] = -999.;
        g_mwdc_pi_dpz[i] = -999.;
        g_mwdc_pi_dangle[i] = -999.;
        g_mwdc_pi_line_dca_mcstart[i] = -999.;
        g_mwdc_pi_dca_reco[i] = -999.;
        g_mwdc_pi_dca_mc[i] = -999.;
    }
}
