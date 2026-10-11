#include <algorithm>
#include <cmath>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

namespace {

bool PassBit(unsigned int value, int bit)
{
  if (bit < 0) {
    return true;
  }
  if (bit > 31) {
    return false;
  }
  return (value & (1u << bit)) != 0u;
}

TH1D *MakeWidthWeightedReference(const TH1 *source, const TH1D *targetBinning, const char *newName)
{
  TH1D *rebinned = static_cast<TH1D *>(targetBinning->Clone(newName));
  rebinned->Reset("ICES");
  rebinned->SetTitle("h1_eff_PID rebinned;p_{T} (GeV/#it{c});efficiency");

  constexpr double kTolerance = 1e-6;
  for (int ibin = 1; ibin <= rebinned->GetNbinsX(); ++ibin) {
    const double low = rebinned->GetXaxis()->GetBinLowEdge(ibin);
    const double high = rebinned->GetXaxis()->GetBinUpEdge(ibin);
    double weightedSum = 0.0;
    double weightSum = 0.0;
    double err2 = 0.0;

    for (int jbin = 1; jbin <= source->GetNbinsX(); ++jbin) {
      const double srcLow = source->GetXaxis()->GetBinLowEdge(jbin);
      const double srcHigh = source->GetXaxis()->GetBinUpEdge(jbin);
      if (srcLow + kTolerance < low || srcHigh > high + kTolerance) {
        continue;
      }
      const double width = srcHigh - srcLow;
      weightedSum += width * source->GetBinContent(jbin);
      err2 += width * width * source->GetBinError(jbin) * source->GetBinError(jbin);
      weightSum += width;
    }

    if (weightSum <= 0.0) {
      continue;
    }
    rebinned->SetBinContent(ibin, weightedSum / weightSum);
    rebinned->SetBinError(ibin, std::sqrt(err2) / weightSum);
  }

  return rebinned;
}

void StyleHist(TH1 *hist, int color, int marker, int lineStyle)
{
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.0);
  hist->SetLineStyle(lineStyle);
  hist->SetLineWidth(2);
}

} // namespace

void compare_filterbit_pid_efficiency(const char *treeFile = "input/LHC25b16_mcrecgen.root",
                                      const char *treeName = "DF_2363808317891232/O2rtdiemtreerec",
                                      const char *referenceFile = "input/LHC25b16_wacc_efficiency.root",
                                      const char *referenceHistName = "h1_eff_PID",
                                      const char *outputRoot = "output/filterbit1_over_bit0_pid_efficiency.root",
                                      const char *outputPdf = "output/filterbit1_over_bit0_pid_efficiency.pdf",
                                      bool requireCorrectAssoc = true,
                                      bool requireNotAmbiguous = true,
                                      const char *ambiguousBranchName = "fIsAmbiguous")
{
  const double ptBins[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0};
  const int nPtBins = sizeof(ptBins) / sizeof(ptBins[0]) - 1;

  TFile input(treeFile, "READ");
  if (input.IsZombie()) {
    std::cerr << "Failed to open tree file: " << treeFile << std::endl;
    return;
  }
  TTree *tree = static_cast<TTree *>(input.Get(treeName));
  if (!tree) {
    std::cerr << "Failed to find tree: " << treeName << std::endl;
    return;
  }

  float fPt = 0.0f;
  unsigned int fFilterMap = 0;
  bool fIsCorrectAssoc = false;
  bool fIsAmbiguous = false;

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("fPt", 1);
  tree->SetBranchStatus("fFilterMap", 1);
  tree->SetBranchAddress("fPt", &fPt);
  tree->SetBranchAddress("fFilterMap", &fFilterMap);

  const bool hasCorrectAssoc = tree->GetBranch("fIsCorrectAssoc") != nullptr;
  if (requireCorrectAssoc && hasCorrectAssoc) {
    tree->SetBranchStatus("fIsCorrectAssoc", 1);
    tree->SetBranchAddress("fIsCorrectAssoc", &fIsCorrectAssoc);
  } else if (requireCorrectAssoc) {
    std::cout << "WARNING: fIsCorrectAssoc branch is missing; correct-association cut is not applied." << std::endl;
  }

  const bool hasAmbiguous = tree->GetBranch(ambiguousBranchName) != nullptr;
  if (requireNotAmbiguous && hasAmbiguous) {
    tree->SetBranchStatus(ambiguousBranchName, 1);
    tree->SetBranchAddress(ambiguousBranchName, &fIsAmbiguous);
  } else if (requireNotAmbiguous) {
    std::cout << "WARNING: ambiguous branch " << ambiguousBranchName
              << " is missing; IsAmbiguous==false cut is not applied." << std::endl;
  }

  TH1D hDen("hFilterBit0Denominator", "Filter bit0 denominator;p_{T} (GeV/#it{c});counts", nPtBins, ptBins);
  TH1D hNum("hFilterBit1Numerator", "Filter bit1 numerator;p_{T} (GeV/#it{c});counts", nPtBins, ptBins);
  TH1D hEff("hEffFilterBit1OverBit0", "Tree filter bit1 / bit0 efficiency;p_{T} (GeV/#it{c});efficiency", nPtBins, ptBins);

  Long64_t selectedBase = 0;
  Long64_t denominator = 0;
  Long64_t numerator = 0;
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if (requireCorrectAssoc && hasCorrectAssoc && !fIsCorrectAssoc) {
      continue;
    }
    if (requireNotAmbiguous && hasAmbiguous && fIsAmbiguous) {
      continue;
    }
    ++selectedBase;

    const bool passBit0 = PassBit(fFilterMap, 0);
    const bool passBit1 = PassBit(fFilterMap, 1);
    if (passBit0) {
      hDen.Fill(fPt);
      ++denominator;
    }
    if (passBit0 && passBit1) {
      hNum.Fill(fPt);
      ++numerator;
    }
  }

  for (int i = 1; i <= hEff.GetNbinsX(); ++i) {
    const double den = hDen.GetBinContent(i);
    if (den <= 0.0) {
      continue;
    }
    const double num = hNum.GetBinContent(i);
    const double eff = num / den;
    hEff.SetBinContent(i, eff);
    hEff.SetBinError(i, std::sqrt(std::max(0.0, eff * (1.0 - eff) / den)));
  }

  TFile refInput(referenceFile, "READ");
  if (refInput.IsZombie()) {
    std::cerr << "Failed to open reference file: " << referenceFile << std::endl;
    return;
  }
  TH1 *hRefNative = static_cast<TH1 *>(refInput.Get(referenceHistName));
  if (!hRefNative) {
    std::cerr << "Failed to find reference histogram: " << referenceHistName << std::endl;
    return;
  }
  TH1D *hRef = MakeWidthWeightedReference(hRefNative, &hEff, "hReferencePidEfficiencyPtRebinned");
  TH1D *hRatio = static_cast<TH1D *>(hEff.Clone("hFilterBitEfficiencyOverReference"));
  hRatio->SetDirectory(nullptr);
  hRatio->Divide(hRef);
  hRatio->SetTitle(";p_{T} (GeV/#it{c});filter bit ratio / h1_eff_PID");

  TFile output(outputRoot, "RECREATE");
  hDen.Write();
  hNum.Write();
  hEff.Write();
  hRefNative->Write("hReferencePidEfficiencyNative");
  hRef->Write();
  hRatio->Write();
  output.Close();

  gStyle->SetOptStat(0);
  StyleHist(&hEff, kRed + 1, 20, 1);
  StyleHist(hRef, kBlack, 20, 2);
  StyleHist(hRatio, kRed + 1, 20, 1);

  TCanvas canvas("canvas", "filter bit PID efficiency", 900, 850);
  TPad top("top", "top", 0.0, 0.32, 1.0, 1.0);
  TPad bottom("bottom", "bottom", 0.0, 0.0, 1.0, 0.32);
  top.SetBottomMargin(0.02);
  top.SetLeftMargin(0.12);
  top.SetRightMargin(0.04);
  bottom.SetTopMargin(0.03);
  bottom.SetBottomMargin(0.30);
  bottom.SetLeftMargin(0.12);
  bottom.SetRightMargin(0.04);
  top.Draw();
  bottom.Draw();

  top.cd();
  hEff.GetYaxis()->SetRangeUser(0.25, 0.95);
  hEff.GetXaxis()->SetLabelSize(0.0);
  hEff.Draw("E");
  hRef->Draw("E SAME");
  TLegend legend(0.18, 0.18, 0.58, 0.34);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.AddEntry(&hEff, "Tree fFilterMap bit1 / bit0", "lp");
  legend.AddEntry(hRef, "LHC25b16_wacc h1_eff_PID", "lp");
  legend.Draw();

  bottom.cd();
  hRatio->GetYaxis()->SetRangeUser(0.82, 1.16);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTitleSize(0.095);
  hRatio->GetYaxis()->SetTitleOffset(0.55);
  hRatio->GetYaxis()->SetLabelSize(0.085);
  hRatio->GetXaxis()->SetTitleSize(0.105);
  hRatio->GetXaxis()->SetLabelSize(0.085);
  hRatio->Draw("E");
  TLine line(hRatio->GetXaxis()->GetXmin(), 1.0, hRatio->GetXaxis()->GetXmax(), 1.0);
  line.SetLineStyle(2);
  line.Draw("SAME");

  canvas.SaveAs(outputPdf);

  std::cout << "Selected base entries: " << selectedBase << std::endl;
  std::cout << "Denominator bit0 entries: " << denominator << std::endl;
  std::cout << "Numerator bit1 entries: " << numerator << std::endl;
  std::cout << "pT bin, bit1/bit0, h1_eff_PID, ratio" << std::endl;
  for (int i = 1; i <= hEff.GetNbinsX(); ++i) {
    std::cout << hEff.GetXaxis()->GetBinLowEdge(i) << "-" << hEff.GetXaxis()->GetBinUpEdge(i)
              << "  " << hEff.GetBinContent(i)
              << "  " << hRef->GetBinContent(i)
              << "  " << hRatio->GetBinContent(i) << std::endl;
  }
  std::cout << "Wrote " << outputRoot << std::endl;
  std::cout << "Wrote " << outputPdf << std::endl;

  delete hRef;
  delete hRatio;
}
