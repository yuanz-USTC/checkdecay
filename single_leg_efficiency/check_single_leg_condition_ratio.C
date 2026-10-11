#include <algorithm>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"

void check_single_leg_condition_ratio(const char *inputFile = "../input/pid_efficiency_maps.root",
                                      const char *outputRoot = "output/single_leg_if_other_pass_over_fail.root",
                                      const char *outputPdf = "output/single_leg_if_other_pass_over_fail.pdf")
{
  TFile input(inputFile, "READ");
  if (input.IsZombie()) {
    std::cerr << "Failed to open input file: " << inputFile << std::endl;
    return;
  }

  TH2D *hPass = static_cast<TH2D *>(input.Get("hPidEffIfOtherPass"));
  TH2D *hFail = static_cast<TH2D *>(input.Get("hPidEffIfOtherFail"));
  if (!hPass || !hFail) {
    std::cerr << "Missing hPidEffIfOtherPass or hPidEffIfOtherFail" << std::endl;
    return;
  }

  TH2D *hRatio = static_cast<TH2D *>(hPass->Clone("hPidEffIfOtherPassOverFail"));
  hRatio->SetDirectory(nullptr);
  hRatio->Reset("ICES");
  hRatio->SetTitle("Single-leg PID efficiency ratio: other pass / other fail;p_{in} (GeV/#it{c});#eta;ratio");

  double minRatio = 1e30;
  double maxRatio = -1e30;
  int validBins = 0;

  for (int ix = 1; ix <= hRatio->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hRatio->GetNbinsY(); ++iy) {
      const double pass = hPass->GetBinContent(ix, iy);
      const double fail = hFail->GetBinContent(ix, iy);
      if (fail <= 0.0) {
        continue;
      }
      const double ratio = pass / fail;
      hRatio->SetBinContent(ix, iy, ratio);
      minRatio = std::min(minRatio, ratio);
      maxRatio = std::max(maxRatio, ratio);
      ++validBins;
    }
  }

  TFile output(outputRoot, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create output file: " << outputRoot << std::endl;
    delete hRatio;
    return;
  }
  hRatio->Write();
  output.Close();

  gStyle->SetOptStat(0);
  TCanvas canvas("canvas", "single leg condition ratio", 900, 750);
  canvas.SetRightMargin(0.15);
  hRatio->GetZaxis()->SetRangeUser(0.85, 1.15);
  hRatio->Draw("COLZ");
  canvas.SaveAs(outputPdf);

  std::cout << "Valid ratio bins: " << validBins << std::endl;
  if (validBins > 0) {
    std::cout << "Ratio min/max: " << minRatio << " / " << maxRatio << std::endl;
  }
  std::cout << "Wrote " << outputRoot << std::endl;
  std::cout << "Wrote " << outputPdf << std::endl;

  delete hRatio;
}
