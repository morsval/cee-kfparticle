#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TTree.h"
#include "TString.h"

namespace
{
void Style(TH1* histogram, int color, int width = 2)
{
  histogram->SetLineColor(color);
  histogram->SetLineWidth(width);
  histogram->SetStats(0);
}

void DrawTriplet(TTree* tree, const char* expression, const char* baseName,
                 const char* title, int bins, double minimum, double maximum,
                 bool hasTruth)
{
  TString allName = TString(baseName) + "All";
  TString truthName = TString(baseName) + "Truth";
  TString backgroundName = TString(baseName) + "Background";
  TH1D* all = new TH1D(allName, title, bins, minimum, maximum);
  TH1D* truth = new TH1D(truthName, title, bins, minimum, maximum);
  TH1D* background = new TH1D(backgroundName, title, bins, minimum, maximum);
  tree->Draw(TString::Format("%s>>%s", expression, allName.Data()), "", "goff");
  if (hasTruth) {
    tree->Draw(TString::Format("%s>>%s", expression, truthName.Data()),
               "truthAvailable==1&&truthMatched==1", "goff");
    tree->Draw(TString::Format("%s>>%s", expression, backgroundName.Data()),
               "truthAvailable==1&&truthMatched==0", "goff");
  }
  Style(all, kBlack);
  Style(truth, kRed + 1);
  Style(background, kBlue + 1);
  all->Draw("hist");
  if (hasTruth) {
    truth->Draw("hist same");
    background->Draw("hist same");
    TLegend* legend = new TLegend(0.57, 0.70, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->AddEntry(all, "all", "l");
    legend->AddEntry(truth, "truth matched", "l");
    legend->AddEntry(background, "truth background", "l");
    legend->Draw();
  }
}
}

void plotCeeV0(const char* input = "CeeV0Lambda.root",
               const char* outputPdf = "CeeV0QA.pdf")
{
  gStyle->SetOptStat(0);
  TFile* file = TFile::Open(input, "READ");
  if (!file || file->IsZombie()) {
    Error("plotCeeV0", "Cannot open %s", input);
    return;
  }
  TTree* tree = dynamic_cast<TTree*>(file->Get("LambdaCandidates"));
  if (!tree) {
    Error("plotCeeV0", "LambdaCandidates is absent");
    file->Close();
    return;
  }
  const bool hasTruth = tree->GetEntries("truthAvailable==1") > 0;

  TCanvas first("cCeeV0First", "CEE Lambda QA", 1200, 900);
  first.Divide(2, 2);
  first.cd(1);
  DrawTriplet(tree, "massGeV", "hPlotMass",
              "#Lambda candidates;inv. mass (GeV/c^{2});candidates",
              260, 0.7, 2.0, hasTruth);
  first.cd(2);
  DrawTriplet(tree, "daughterDcaCm", "hPlotDca",
              "Daughter-pair DCA;DCA (cm);candidates",
              200, 0., 10., hasTruth);
  first.cd(3);
  DrawTriplet(tree, "decayLengthCm", "hPlotLength",
              "Decay length;L (cm);candidates", 240, 0., 60., hasTruth);
  first.cd(4);
  DrawTriplet(tree, "pointingCosine", "hPlotPointing",
              "Pointing cosine;cos(#theta);candidates", 220, -1.1, 1.1,
              hasTruth);
  first.Print(TString(outputPdf) + "(");

  TCanvas second("cCeeV0Second", "CEE Lambda chi2 QA", 1200, 900);
  second.Divide(2, 2);
  second.cd(1);
  DrawTriplet(tree, "geometricChi2Ndf", "hPlotGeo",
              "Geometric #chi^{2}/ndf;#chi^{2}/ndf;candidates",
              250, 0., 50., hasTruth);
  second.cd(2);
  DrawTriplet(tree, "topologicalChi2Ndf", "hPlotTopo",
              "STAR-like total topological #chi^{2}/ndf;#chi^{2}/ndf;candidates",
              250, 0., 50., hasTruth);
  second.cd(3);
  DrawTriplet(tree, "primaryVertexConstraintChi2Ndf", "hPlotPvConstraint",
              "PV-only constraint #chi^{2}/ndf;#chi^{2}/ndf;candidates",
              250, 0., 50., hasTruth);
  second.cd(4);
  DrawTriplet(tree, "min(daughterPvChi20,daughterPvChi21)",
              "hPlotMinimumDaughterPvChi2",
              "Minimum daughter PV #chi^{2};min #chi^{2};candidates",
              250, 0., 200., hasTruth);
  second.Print(TString(outputPdf) + ")");
  file->Close();
}
