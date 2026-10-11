#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"

namespace {

std::vector<double> AllEdges(const TAxis *axis)
{
  std::vector<double> edges;
  edges.reserve(axis->GetNbins() + 1);
  for (int ibin = 1; ibin <= axis->GetNbins(); ++ibin) {
    edges.push_back(axis->GetBinLowEdge(ibin));
  }
  edges.push_back(axis->GetBinUpEdge(axis->GetNbins()));
  return edges;
}

TH2D *ReadHist2D(const char *fileName, const char *histName, const char *cloneName)
{
  TFile input(fileName, "READ");
  if (input.IsZombie()) {
    std::cerr << "Failed to open " << fileName << std::endl;
    return nullptr;
  }

  TH2 *hist = dynamic_cast<TH2 *>(input.Get(histName));
  if (!hist) {
    std::cerr << "Missing 2D histogram " << histName << " in " << fileName << std::endl;
    return nullptr;
  }

  TH2D *copy = dynamic_cast<TH2D *>(hist->Clone(cloneName));
  if (!copy) {
    const std::vector<double> xEdges = AllEdges(hist->GetXaxis());
    const std::vector<double> yEdges = AllEdges(hist->GetYaxis());
    copy = new TH2D(cloneName, hist->GetTitle(), hist->GetNbinsX(), xEdges.data(), hist->GetNbinsY(),
                    yEdges.data());
    for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
      for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
        copy->SetBinContent(ix, iy, hist->GetBinContent(ix, iy));
        copy->SetBinError(ix, iy, hist->GetBinError(ix, iy));
      }
    }
  }
  copy->SetDirectory(nullptr);
  return copy;
}

std::vector<double> MakeEdges(const TAxis *axis, int firstBin, int lastBin)
{
  std::vector<double> edges;
  edges.reserve(lastBin - firstBin + 2);
  for (int ibin = firstBin; ibin <= lastBin; ++ibin) {
    edges.push_back(axis->GetBinLowEdge(ibin));
  }
  edges.push_back(axis->GetBinUpEdge(lastBin));
  return edges;
}

bool FindOverlapBins(const TAxis *axis, const TAxis *referenceAxis, int &firstBin, int &lastBin)
{
  firstBin = 0;
  lastBin = -1;
  constexpr double eps = 1e-9;

  for (int ibin = 1; ibin <= axis->GetNbins(); ++ibin) {
    const double center = axis->GetBinCenter(ibin);
    const int refBin = referenceAxis->FindFixBin(center);
    if (refBin < 1 || refBin > referenceAxis->GetNbins()) {
      continue;
    }
    if (center < referenceAxis->GetBinLowEdge(refBin) - eps ||
        center >= referenceAxis->GetBinUpEdge(refBin) + eps) {
      continue;
    }
    if (firstBin == 0) {
      firstBin = ibin;
    }
    lastBin = ibin;
  }

  return firstBin > 0 && lastBin >= firstBin;
}

double RatioError(double num, double numErr, double den, double denErr)
{
  if (den == 0.0) {
    return 0.0;
  }
  const double ratio = num / den;
  double rel2 = 0.0;
  if (num != 0.0) {
    rel2 += (numErr / num) * (numErr / num);
  }
  if (den != 0.0) {
    rel2 += (denErr / den) * (denErr / den);
  }
  return std::abs(ratio) * std::sqrt(rel2);
}

void PrintAxis(const char *label, const TH2D *hist)
{
  std::cout << label << ": "
            << "nx=" << hist->GetNbinsX() << " x=[" << hist->GetXaxis()->GetXmin() << ", "
            << hist->GetXaxis()->GetXmax() << "] dx~=" << hist->GetXaxis()->GetBinWidth(1)
            << ", ny=" << hist->GetNbinsY() << " y=[" << hist->GetYaxis()->GetXmin() << ", "
            << hist->GetYaxis()->GetXmax() << "] dy~=" << hist->GetYaxis()->GetBinWidth(1)
            << std::endl;
}

} // namespace

void compare_pid_effavg_to_reference(const char *pidMapFile = "./input/pid_efficiency_maps.root",
                                     const char *pidMapHist = "hPidEffAvg",
                                     const char *referenceFile = "./input/singlegeff.root",
                                     const char *referenceHist = "h_reference",
                                     const char *outputRoot = "./output/hPidEffAvg_over_h_reference.root",
                                     const char *outputPdf = "./output/hPidEffAvg_over_h_reference.pdf",
                                     const char *outputPng = "./output/hPidEffAvg_over_h_reference.png")
{
  gStyle->SetOptStat(0);

  TH2D *hPid = ReadHist2D(pidMapFile, pidMapHist, "hPidEffAvg_copy");
  TH2D *hReference = ReadHist2D(referenceFile, referenceHist, "hReference_copy");
  if (!hPid || !hReference) {
    delete hPid;
    delete hReference;
    return;
  }

  PrintAxis("hPidEffAvg", hPid);
  PrintAxis("h_reference", hReference);

  int firstX = 0;
  int lastX = -1;
  int firstY = 0;
  int lastY = -1;
  if (!FindOverlapBins(hPid->GetXaxis(), hReference->GetXaxis(), firstX, lastX) ||
      !FindOverlapBins(hPid->GetYaxis(), hReference->GetYaxis(), firstY, lastY)) {
    std::cerr << "No overlapping bin centers between the two histograms." << std::endl;
    delete hPid;
    delete hReference;
    return;
  }

  const std::vector<double> xEdges = MakeEdges(hPid->GetXaxis(), firstX, lastX);
  const std::vector<double> yEdges = MakeEdges(hPid->GetYaxis(), firstY, lastY);
  TH2D *hRatio = new TH2D("hPidEffAvgOverReference", "hPidEffAvg / h_reference;p_{in} (GeV/#it{c});#eta;ratio",
                          static_cast<int>(xEdges.size()) - 1, xEdges.data(),
                          static_cast<int>(yEdges.size()) - 1, yEdges.data());
  hRatio->SetDirectory(nullptr);

  double minRatio = std::numeric_limits<double>::max();
  double maxRatio = -std::numeric_limits<double>::max();
  int validBins = 0;
  int skippedBins = 0;
  for (int ix = firstX; ix <= lastX; ++ix) {
    const double x = hPid->GetXaxis()->GetBinCenter(ix);
    const int refX = hReference->GetXaxis()->FindFixBin(x);
    for (int iy = firstY; iy <= lastY; ++iy) {
      const double y = hPid->GetYaxis()->GetBinCenter(iy);
      const int refY = hReference->GetYaxis()->FindFixBin(y);
      const double den = hReference->GetBinContent(refX, refY);
      if (den <= 0.0) {
        ++skippedBins;
        continue;
      }

      const double num = hPid->GetBinContent(ix, iy);
      const double ratio = num / den;
      const int outX = ix - firstX + 1;
      const int outY = iy - firstY + 1;
      hRatio->SetBinContent(outX, outY, ratio);
      hRatio->SetBinError(outX, outY,
                          RatioError(num, hPid->GetBinError(ix, iy), den, hReference->GetBinError(refX, refY)));
      minRatio = std::min(minRatio, ratio);
      maxRatio = std::max(maxRatio, ratio);
      ++validBins;
    }
  }

  TFile output(outputRoot, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create " << outputRoot << std::endl;
    delete hPid;
    delete hReference;
    delete hRatio;
    return;
  }
  hPid->Write("hPidEffAvg");
  hReference->Write("h_reference");
  hRatio->Write();
  output.Close();

  TCanvas canvas("canvas", "hPidEffAvg over h_reference", 900, 750);
  canvas.SetRightMargin(0.15);
  if (validBins > 0) {
    const double margin = std::max(0.02, 0.10 * (maxRatio - minRatio));
    hRatio->GetZaxis()->SetRangeUser(minRatio - margin, maxRatio + margin);
  }
  hRatio->Draw("COLZ");
  canvas.SaveAs(outputPdf);
  canvas.SaveAs(outputPng);

  std::cout << "Overlap bins from hPidEffAvg: x " << firstX << "-" << lastX << ", y " << firstY << "-"
            << lastY << std::endl;
  std::cout << "Valid ratio bins: " << validBins << ", skipped denominator<=0 bins: " << skippedBins << std::endl;
  if (validBins > 0) {
    std::cout << "Ratio min/max: " << minRatio << " / " << maxRatio << std::endl;
  }
  std::cout << "Wrote " << outputRoot << std::endl;
  std::cout << "Wrote " << outputPdf << std::endl;
  std::cout << "Wrote " << outputPng << std::endl;

  delete hPid;
  delete hReference;
  delete hRatio;
}
