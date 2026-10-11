#include <algorithm>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"

void compare_leg_eff_to_external(const char *newMapFile = "../input/pid_efficiency_maps.root",
                                 const char *newMapHist = "hPidEffAvg",
                                 const char *externalFile = "../input/singlegeff.root",
                                 const char *externalHist = "hEff",
                                 const char *outputRoot = "output/tree_leg_eff_over_singlegeff.root",
                                 const char *outputPdf = "output/tree_leg_eff_over_singlegeff.pdf")
{
  TFile inputNew(newMapFile, "READ");
  if (inputNew.IsZombie()) {
    std::cerr << "Failed to open new map file: " << newMapFile << std::endl;
    return;
  }
  TH2D *hNew = static_cast<TH2D *>(inputNew.Get(newMapHist));
  if (!hNew) {
    std::cerr << "Missing new map histogram: " << newMapHist << std::endl;
    return;
  }

  TFile inputExternal(externalFile, "READ");
  if (inputExternal.IsZombie()) {
    std::cerr << "Failed to open external file: " << externalFile << std::endl;
    return;
  }
  TH2D *hExternal = static_cast<TH2D *>(inputExternal.Get(externalHist));
  if (!hExternal) {
    std::cerr << "Missing external histogram: " << externalHist << std::endl;
    return;
  }

  if (hNew->GetNbinsX() != hExternal->GetNbinsX() || hNew->GetNbinsY() != hExternal->GetNbinsY()) {
    std::cerr << "Histogram bin counts do not match." << std::endl;
    return;
  }

  TH2D *hRatio = static_cast<TH2D *>(hNew->Clone("hTreeLegEffOverExternalEff"));
  hRatio->SetDirectory(nullptr);
  hRatio->Reset("ICES");
  hRatio->SetTitle("Tree-derived leg efficiency / external hEff;p_{in} (GeV/#it{c});#eta;ratio");

  double minRatio = 1e30;
  double maxRatio = -1e30;
  int validBins = 0;
  for (int ix = 1; ix <= hRatio->GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hRatio->GetNbinsY(); ++iy) {
      const double den = hExternal->GetBinContent(ix, iy);
      if (den <= 0.0) {
        continue;
      }
      const double ratio = hNew->GetBinContent(ix, iy) / den;
      hRatio->SetBinContent(ix, iy, ratio);
      minRatio = std::min(minRatio, ratio);
      maxRatio = std::max(maxRatio, ratio);
      ++validBins;
    }
  }

  TFile output(outputRoot, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create output: " << outputRoot << std::endl;
    delete hRatio;
    return;
  }
  hRatio->Write();
  output.Close();

  gStyle->SetOptStat(0);
  TCanvas canvas("canvas", "tree leg eff over external hEff", 900, 750);
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
