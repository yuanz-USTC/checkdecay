#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"

namespace {

TH1D *ReadEfficiency(const char *fileName, const char *histName, const char *newName)
{
  TFile input(fileName, "READ");
  if (input.IsZombie()) {
    std::cerr << "Failed to open " << fileName << std::endl;
    return nullptr;
  }
  TH1D *hist = static_cast<TH1D *>(input.Get(histName));
  if (!hist) {
    std::cerr << "Missing histogram " << histName << " in " << fileName << std::endl;
    return nullptr;
  }
  TH1D *copy = static_cast<TH1D *>(hist->Clone(newName));
  copy->SetDirectory(nullptr);
  return copy;
}

TH1 *ReadNativeReference(const char *fileName, const char *histName, const char *newName)
{
  TFile input(fileName, "READ");
  if (input.IsZombie()) {
    std::cerr << "Failed to open " << fileName << std::endl;
    return nullptr;
  }
  TH1 *hist = static_cast<TH1 *>(input.Get(histName));
  if (!hist) {
    std::cerr << "Missing histogram " << histName << " in " << fileName << std::endl;
    return nullptr;
  }
  TH1 *copy = static_cast<TH1 *>(hist->Clone(newName));
  copy->SetDirectory(nullptr);
  return copy;
}

TH1D *MakeWidthWeightedReference(const TH1 *source, const TH1D *targetBinning, const char *newName)
{
  if (!source || !targetBinning) {
    return nullptr;
  }

  TH1D *rebinned = static_cast<TH1D *>(targetBinning->Clone(newName));
  rebinned->Reset("ICES");
  rebinned->SetTitle("LHC25b16 wacc h1_eff_PID;p_{T} (GeV/#it{c});PID efficiency");

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

TH1D *MakeRatio(const TH1D *hist, const TH1D *reference, const char *name)
{
  TH1D *ratio = static_cast<TH1D *>(hist->Clone(name));
  ratio->SetDirectory(nullptr);
  ratio->Divide(reference);
  ratio->SetTitle(";p_{T} (GeV/#it{c});ratio to h1_eff_PID");
  return ratio;
}

} // namespace

void plot_pid_efficiency_four_way_comparison(const char *outputPdf = "output/jpsi_pid_efficiency_four_way_comparison.pdf",
                                             const char *outputRoot = "output/jpsi_pid_efficiency_four_way_comparison.root")
{
  gStyle->SetOptStat(0);

  TH1D *hTreeCorr = ReadEfficiency("output/jpsi_eff_tree_corr.root", "hJpsiPidEfficiencyPt", "hTreeCorr");
  TH1D *hTreeAvg = ReadEfficiency("output/jpsi_eff_tree_avg.root", "hJpsiPidEfficiencyPt", "hTreeAvg");
  TH1D *hTreeHReference = ReadEfficiency("output/jpsi_eff_tree_hreference.root", "hJpsiPidEfficiencyPt", "hTreeHReference");
  TH1D *hSimCorr = ReadEfficiency("output/jpsi_eff_sim_weighted_corr.root", "hJpsiPidEfficiencyPt", "hSimCorr");
  TH1D *hSimAvg = ReadEfficiency("output/jpsi_eff_sim_weighted_avg.root", "hJpsiPidEfficiencyPt", "hSimAvg");
  TH1 *hRefNative = ReadNativeReference("input/LHC25b16_wacc_efficiency.root", "h1_eff_PID", "hRefNative");

  if (!hTreeCorr || !hTreeAvg || !hTreeHReference || !hSimCorr || !hSimAvg || !hRefNative) {
    return;
  }

  TH1D *hRef = MakeWidthWeightedReference(hRefNative, hTreeCorr, "hReferencePidEfficiencyPtRebinned");
  if (!hRef) {
    return;
  }

  StyleHist(hTreeCorr, kRed + 1, 20, 1);
  StyleHist(hTreeAvg, kOrange + 7, 21, 1);
  StyleHist(hTreeHReference, kMagenta + 2, 22, 1);
  StyleHist(hSimCorr, kBlue + 1, 24, 1);
  StyleHist(hSimAvg, kAzure + 7, 25, 1);
  StyleHist(hRef, kBlack, 20, 2);

  TH1D *rTreeCorr = MakeRatio(hTreeCorr, hRef, "hRatioTreeCorrToReference");
  TH1D *rTreeAvg = MakeRatio(hTreeAvg, hRef, "hRatioTreeAvgToReference");
  TH1D *rTreeHReference = MakeRatio(hTreeHReference, hRef, "hRatioTreeHReferenceToReference");
  TH1D *rSimCorr = MakeRatio(hSimCorr, hRef, "hRatioSimCorrToReference");
  TH1D *rSimAvg = MakeRatio(hSimAvg, hRef, "hRatioSimAvgToReference");

  StyleHist(rTreeCorr, kRed + 1, 20, 1);
  StyleHist(rTreeAvg, kOrange + 7, 21, 1);
  StyleHist(rTreeHReference, kMagenta + 2, 22, 1);
  StyleHist(rSimCorr, kBlue + 1, 24, 1);
  StyleHist(rSimAvg, kAzure + 7, 25, 1);

  TCanvas canvas("canvas", "J/psi PID efficiency comparison", 900, 850);
  TPad padTop("padTop", "padTop", 0.0, 0.32, 1.0, 1.0);
  TPad padBottom("padBottom", "padBottom", 0.0, 0.0, 1.0, 0.32);
  padTop.SetBottomMargin(0.02);
  padTop.SetLeftMargin(0.12);
  padTop.SetRightMargin(0.04);
  padBottom.SetTopMargin(0.03);
  padBottom.SetBottomMargin(0.30);
  padBottom.SetLeftMargin(0.12);
  padBottom.SetRightMargin(0.04);
  padTop.Draw();
  padBottom.Draw();

  padTop.cd();
  hTreeCorr->SetTitle(";p_{T} (GeV/#it{c});J/#psi PID efficiency");
  hTreeCorr->GetYaxis()->SetRangeUser(0.25, 0.95);
  hTreeCorr->GetYaxis()->SetTitleSize(0.050);
  hTreeCorr->GetYaxis()->SetLabelSize(0.043);
  hTreeCorr->GetXaxis()->SetLabelSize(0.0);
  hTreeCorr->Draw("E");
  hTreeAvg->Draw("E SAME");
  hTreeHReference->Draw("E SAME");
  hSimCorr->Draw("E SAME");
  hSimAvg->Draw("E SAME");
  hRef->Draw("E SAME");

  TLegend legend(0.16, 0.13, 0.70, 0.42);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.030);
  legend.AddEntry(hTreeCorr, "Input tree + correlated maps", "lp");
  legend.AddEntry(hTreeAvg, "Input tree + average maps", "lp");
  legend.AddEntry(hTreeHReference, "Input tree + h_reference map", "lp");
  legend.AddEntry(hSimCorr, "Sim decay weighted by input p_{T}-y + correlated maps", "lp");
  legend.AddEntry(hSimAvg, "Sim decay weighted by input p_{T}-y + average maps", "lp");
  legend.AddEntry(hRef, "LHC25b16_wacc h1_eff_PID", "lp");
  legend.Draw();

  padBottom.cd();
  rTreeCorr->GetYaxis()->SetRangeUser(0.82, 1.16);
  rTreeCorr->GetYaxis()->SetNdivisions(505);
  rTreeCorr->GetYaxis()->SetTitleSize(0.095);
  rTreeCorr->GetYaxis()->SetTitleOffset(0.55);
  rTreeCorr->GetYaxis()->SetLabelSize(0.085);
  rTreeCorr->GetXaxis()->SetTitleSize(0.105);
  rTreeCorr->GetXaxis()->SetLabelSize(0.085);
  rTreeCorr->Draw("E");
  rTreeAvg->Draw("E SAME");
  rTreeHReference->Draw("E SAME");
  rSimCorr->Draw("E SAME");
  rSimAvg->Draw("E SAME");

  TLine line(rTreeCorr->GetXaxis()->GetXmin(), 1.0, rTreeCorr->GetXaxis()->GetXmax(), 1.0);
  line.SetLineColor(kBlack);
  line.SetLineStyle(2);
  line.Draw("SAME");

  canvas.SaveAs(outputPdf);

  TFile output(outputRoot, "RECREATE");
  hTreeCorr->Write("hJpsiPidEfficiencyPtTreeCorrelatedMaps");
  hTreeAvg->Write("hJpsiPidEfficiencyPtTreeAverageMaps");
  hTreeHReference->Write("hJpsiPidEfficiencyPtTreeHReferenceMap");
  hSimCorr->Write("hJpsiPidEfficiencyPtSimWeightedCorrelatedMaps");
  hSimAvg->Write("hJpsiPidEfficiencyPtSimWeightedAverageMaps");
  hRefNative->Write("hReferencePidEfficiencyNative");
  hRef->Write("hReferencePidEfficiencyPtRebinned");
  rTreeCorr->Write();
  rTreeAvg->Write();
  rTreeHReference->Write();
  rSimCorr->Write();
  rSimAvg->Write();
  output.Close();

  std::cout << "Wrote " << outputPdf << std::endl;
  std::cout << "Wrote " << outputRoot << std::endl;

  delete hTreeCorr;
  delete hTreeAvg;
  delete hTreeHReference;
  delete hSimCorr;
  delete hSimAvg;
  delete hRefNative;
  delete hRef;
  delete rTreeCorr;
  delete rTreeAvg;
  delete rTreeHReference;
  delete rSimCorr;
  delete rSimAvg;
}
