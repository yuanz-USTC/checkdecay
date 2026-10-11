#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TStyle.h"

#include "./run_checkdecay.C"

void run_jpsi_eff_tree_h_reference()
{
  const char *outputRoot = "./output/jpsi_eff_tree_h_reference.root";
  const char *outputPdf = "./output/jpsi_eff_tree_h_reference_pt.pdf";
  const char *outputPng = "./output/jpsi_eff_tree_h_reference_pt.png";

  run_checkdecay("tree",
                 1000000000000LL,
                 "./input/singlegeff.root",
                 outputRoot,
                 80,
                 0.0,
                 20.0,
                 80,
                 -4.0,
                 4.0,
                 "./input/LHC25b16_mcrecgen.root",
                 "DF_2363808317891232/O2rtdiemtreerec",
                 12345,
                 0.0,
                 0.0,
                 0.0,
                 "h_reference",
                 "h_reference",
                 "h_reference",
                 "rec",
                 0,
                 0,
                 true);

  TFile input(outputRoot, "READ");
  if (input.IsZombie()) {
    return;
  }
  TH1D *hEfficiencyPt = static_cast<TH1D *>(input.Get("hJpsiPidEfficiencyPt"));
  if (!hEfficiencyPt) {
    return;
  }

  gStyle->SetOptStat(0);
  TCanvas canvas("canvas", "J/psi efficiency from h_reference", 850, 650);
  canvas.SetLeftMargin(0.12);
  canvas.SetBottomMargin(0.12);
  hEfficiencyPt->SetTitle("Decay-tree J/#psi PID efficiency from h_reference;p_{T} (GeV/#it{c});efficiency");
  hEfficiencyPt->SetLineColor(kBlue + 1);
  hEfficiencyPt->SetMarkerColor(kBlue + 1);
  hEfficiencyPt->SetMarkerStyle(20);
  hEfficiencyPt->SetLineWidth(2);
  hEfficiencyPt->GetYaxis()->SetRangeUser(0.0, 1.0);
  hEfficiencyPt->Draw("E");
  canvas.SaveAs(outputPdf);
  canvas.SaveAs(outputPng);
}
