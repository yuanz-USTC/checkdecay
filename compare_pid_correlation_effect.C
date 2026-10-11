#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TH1.h"

namespace {

TH1D *ReadPtEfficiency(TFile &file, const char *histName, const char *newName)
{
  TH1D *hist = static_cast<TH1D *>(file.Get(histName));
  if (!hist) {
    std::cerr << "Missing histogram: " << histName << std::endl;
    return nullptr;
  }
  TH1D *copy = static_cast<TH1D *>(hist->Clone(newName));
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
  rebinned->SetTitle("Reference PID efficiency rebinned to comparison p_{T};p_{T} (GeV/#it{c});efficiency");

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

} // namespace

void compare_pid_correlation_effect(const char *correlatedFile = "output/jpsi_pid_efficiency_rec.root",
                                    const char *averageMapFile = "output/jpsi_pid_efficiency_rec_avgmaps.root",
                                    const char *referenceFile = "input/LHC25b16_wacc_efficiency.root",
                                    const char *referenceHistName = "h1_eff_PID",
                                    const char *outputFile = "output/jpsi_pid_efficiency_pt_correlation_comparison.root")
{
  TFile inputCorrelated(correlatedFile, "READ");
  if (inputCorrelated.IsZombie()) {
    std::cerr << "Failed to open correlated input: " << correlatedFile << std::endl;
    return;
  }

  TFile inputAverage(averageMapFile, "READ");
  if (inputAverage.IsZombie()) {
    std::cerr << "Failed to open average-map input: " << averageMapFile << std::endl;
    return;
  }

  TH1D *hCorrelated = ReadPtEfficiency(inputCorrelated, "hJpsiPidEfficiencyPt", "hJpsiPidEfficiencyPtCorrelatedMaps");
  TH1D *hAverage = ReadPtEfficiency(inputAverage, "hJpsiPidEfficiencyPt", "hJpsiPidEfficiencyPtAverageMaps");
  if (!hCorrelated || !hAverage) {
    delete hCorrelated;
    delete hAverage;
    return;
  }

  TFile inputReference(referenceFile, "READ");
  TH1 *hReferenceNative = nullptr;
  TH1D *hReferenceRebinned = nullptr;
  if (inputReference.IsZombie()) {
    std::cerr << "Failed to open reference input: " << referenceFile << std::endl;
  } else {
    TH1 *reference = static_cast<TH1 *>(inputReference.Get(referenceHistName));
    if (!reference) {
      std::cerr << "Missing reference histogram: " << referenceHistName << std::endl;
    } else {
      hReferenceNative = static_cast<TH1 *>(reference->Clone("hReferencePidEfficiencyNative"));
      hReferenceNative->SetDirectory(nullptr);
      hReferenceRebinned = MakeWidthWeightedReference(reference, hCorrelated, "hReferencePidEfficiencyPtRebinned");
    }
  }

  TH1D *hDiff = static_cast<TH1D *>(hAverage->Clone("hJpsiPidEfficiencyPtAverageMinusCorrelated"));
  hDiff->SetDirectory(nullptr);
  hDiff->Add(hCorrelated, -1.0);
  hDiff->SetTitle("J/#psi PID efficiency: average maps - correlated maps;p_{T} (GeV/#it{c});#Delta efficiency");

  TH1D *hRatio = static_cast<TH1D *>(hAverage->Clone("hJpsiPidEfficiencyPtAverageOverCorrelated"));
  hRatio->SetDirectory(nullptr);
  hRatio->Divide(hCorrelated);
  hRatio->SetTitle("J/#psi PID efficiency: average maps / correlated maps;p_{T} (GeV/#it{c});ratio");

  TH1D *hCorrMinusRef = nullptr;
  TH1D *hAvgMinusRef = nullptr;
  TH1D *hCorrOverRef = nullptr;
  TH1D *hAvgOverRef = nullptr;
  if (hReferenceRebinned) {
    hCorrMinusRef = static_cast<TH1D *>(hCorrelated->Clone("hJpsiPidEfficiencyPtCorrelatedMinusReference"));
    hCorrMinusRef->SetDirectory(nullptr);
    hCorrMinusRef->Add(hReferenceRebinned, -1.0);
    hCorrMinusRef->SetTitle("Correlated-map PID efficiency - reference;p_{T} (GeV/#it{c});#Delta efficiency");

    hAvgMinusRef = static_cast<TH1D *>(hAverage->Clone("hJpsiPidEfficiencyPtAverageMinusReference"));
    hAvgMinusRef->SetDirectory(nullptr);
    hAvgMinusRef->Add(hReferenceRebinned, -1.0);
    hAvgMinusRef->SetTitle("Average-map PID efficiency - reference;p_{T} (GeV/#it{c});#Delta efficiency");

    hCorrOverRef = static_cast<TH1D *>(hCorrelated->Clone("hJpsiPidEfficiencyPtCorrelatedOverReference"));
    hCorrOverRef->SetDirectory(nullptr);
    hCorrOverRef->Divide(hReferenceRebinned);
    hCorrOverRef->SetTitle("Correlated-map PID efficiency / reference;p_{T} (GeV/#it{c});ratio");

    hAvgOverRef = static_cast<TH1D *>(hAverage->Clone("hJpsiPidEfficiencyPtAverageOverReference"));
    hAvgOverRef->SetDirectory(nullptr);
    hAvgOverRef->Divide(hReferenceRebinned);
    hAvgOverRef->SetTitle("Average-map PID efficiency / reference;p_{T} (GeV/#it{c});ratio");
  }

  TFile output(outputFile, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create output file: " << outputFile << std::endl;
    delete hCorrelated;
    delete hAverage;
    delete hDiff;
    delete hRatio;
    delete hReferenceNative;
    delete hReferenceRebinned;
    delete hCorrMinusRef;
    delete hAvgMinusRef;
    delete hCorrOverRef;
    delete hAvgOverRef;
    return;
  }

  hCorrelated->Write();
  hAverage->Write();
  if (hReferenceNative) {
    hReferenceNative->Write();
  }
  if (hReferenceRebinned) {
    hReferenceRebinned->Write();
  }
  hDiff->Write();
  hRatio->Write();
  if (hCorrMinusRef) {
    hCorrMinusRef->Write();
    hAvgMinusRef->Write();
    hCorrOverRef->Write();
    hAvgOverRef->Write();
  }
  output.Close();

  std::cout << "pT bin, correlated, average maps, reference, avg-corr, corr-ref, avg-ref, avg/corr, corr/ref, avg/ref" << std::endl;
  for (int i = 1; i <= hCorrelated->GetNbinsX(); ++i) {
    const double corr = hCorrelated->GetBinContent(i);
    const double avg = hAverage->GetBinContent(i);
    const double ref = hReferenceRebinned ? hReferenceRebinned->GetBinContent(i) : 0.0;
    const double diff = hDiff->GetBinContent(i);
    const double ratio = hRatio->GetBinContent(i);
    std::cout << hCorrelated->GetXaxis()->GetBinLowEdge(i) << "-"
              << hCorrelated->GetXaxis()->GetBinUpEdge(i) << "  "
              << corr << "  " << avg << "  " << ref << "  "
              << diff << "  "
              << (hCorrMinusRef ? hCorrMinusRef->GetBinContent(i) : 0.0) << "  "
              << (hAvgMinusRef ? hAvgMinusRef->GetBinContent(i) : 0.0) << "  "
              << ratio << "  "
              << (hCorrOverRef ? hCorrOverRef->GetBinContent(i) : 0.0) << "  "
              << (hAvgOverRef ? hAvgOverRef->GetBinContent(i) : 0.0) << std::endl;
  }

  std::cout << "Wrote comparison file: " << outputFile << std::endl;

  delete hCorrelated;
  delete hAverage;
  delete hDiff;
  delete hRatio;
  delete hReferenceNative;
  delete hReferenceRebinned;
  delete hCorrMinusRef;
  delete hAvgMinusRef;
  delete hCorrOverRef;
  delete hAvgOverRef;
}
