#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "TAxis.h"
#include "TFile.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"

namespace {

struct LegInfo {
  float p = 0.0f;
  float eta = 0.0f;
  float nSigmaEl = 999.0f;
  float nSigmaPi = -999.0f;
  float nSigmaPr = -999.0f;
};

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

bool PassPid(const LegInfo &leg,
             double nSigmaElMin,
             double nSigmaElMax,
             double nSigmaPiMin,
             double nSigmaPrMin)
{
  return leg.nSigmaEl >= nSigmaElMin &&
         leg.nSigmaEl <= nSigmaElMax &&
         leg.nSigmaPi > nSigmaPiMin &&
         leg.nSigmaPr > nSigmaPrMin;
}

void FillLeg(const LegInfo &leg,
             bool legPassPid,
             bool otherLegPassPid,
             TH2D &hDenAvg,
             TH2D &hNumAvg,
             TH2D &hDenIfOtherPass,
             TH2D &hNumIfOtherPass,
             TH2D &hDenIfOtherFail,
             TH2D &hNumIfOtherFail)
{
  hDenAvg.Fill(leg.p, leg.eta);
  if (legPassPid) {
    hNumAvg.Fill(leg.p, leg.eta);
  }

  if (otherLegPassPid) {
    hDenIfOtherPass.Fill(leg.p, leg.eta);
    if (legPassPid) {
      hNumIfOtherPass.Fill(leg.p, leg.eta);
    }
  } else {
    hDenIfOtherFail.Fill(leg.p, leg.eta);
    if (legPassPid) {
      hNumIfOtherFail.Fill(leg.p, leg.eta);
    }
  }
}

void MakeEfficiency(TH2D &hEff, const TH2D &hNum, const TH2D &hDen)
{
  hEff.Reset("ICES");
  for (int ix = 1; ix <= hEff.GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hEff.GetNbinsY(); ++iy) {
      const double den = hDen.GetBinContent(ix, iy);
      if (den <= 0.0) {
        continue;
      }
      const double num = hNum.GetBinContent(ix, iy);
      const double eff = num / den;
      hEff.SetBinContent(ix, iy, eff);
      hEff.SetBinError(ix, iy, std::sqrt(std::max(0.0, eff * (1.0 - eff) / den)));
    }
  }
}

std::vector<double> BuildDefaultPinBins()
{
  std::vector<double> edges;
  for (int i = 0; i <= 20; ++i) {
    edges.push_back(0.2 * i);
  }
  edges.push_back(6.0);
  edges.push_back(10.0);
  edges.push_back(15.0);
  edges.push_back(20.0);
  return edges;
}

std::vector<double> BuildUniformBins(double min, double max, int nBins)
{
  std::vector<double> edges;
  edges.reserve(nBins + 1);
  const double width = (max - min) / nBins;
  for (int i = 0; i <= nBins; ++i) {
    edges.push_back(min + width * i);
  }
  return edges;
}

std::vector<double> BuildAxisBins(const TAxis *axis)
{
  std::vector<double> edges;
  if (!axis) {
    return edges;
  }
  edges.reserve(axis->GetNbins() + 1);
  for (int i = 1; i <= axis->GetNbins(); ++i) {
    edges.push_back(axis->GetBinLowEdge(i));
  }
  edges.push_back(axis->GetBinUpEdge(axis->GetNbins()));
  return edges;
}

bool LoadReferenceBinning(const char *fileName,
                          const char *histName,
                          std::vector<double> &pBins,
                          std::vector<double> &etaBins)
{
  if (!fileName || !histName || TString(fileName).IsNull()) {
    return false;
  }
  TFile input(fileName, "READ");
  if (input.IsZombie()) {
    std::cerr << "Failed to open reference binning file: " << fileName << std::endl;
    return false;
  }
  TH2D *hist = static_cast<TH2D *>(input.Get(histName));
  if (!hist) {
    std::cerr << "Failed to find reference TH2D: " << histName << " in " << fileName << std::endl;
    return false;
  }
  pBins = BuildAxisBins(hist->GetXaxis());
  etaBins = BuildAxisBins(hist->GetYaxis());
  return !pBins.empty() && !etaBins.empty();
}

} // namespace

void make_single_leg_efficiency_maps(const char *inputFile = "../input/LHC25b16_mcrecgen.root",
                                     const char *outputFile = "../input/pid_efficiency_maps.root",
                                     const char *directoryName = "DF_2363808317891232",
                                     const char *treeName = "O2rtdiemtreerec",
                                     Long64_t maxEntries = -1,
                                     int filterMapBit = 0,
                                     bool requireCorrectAssoc = true,
                                     int nPBins = -1,
                                     double pMin = 0.0,
                                     double pMax = 20.0,
                                     int nEtaBins = -1,
                                     double etaMin = -0.8,
                                     double etaMax = 0.8,
                                     double nSigmaElMin = -2.5,
                                     double nSigmaElMax = 3.0,
                                     double nSigmaPiMin = 3.0,
                                     double nSigmaPrMin = 3.0,
                                     const char *referenceBinningFile = "",
                                     const char *referenceBinningHist = "hEff")
{
  TFile input(inputFile, "READ");
  if (input.IsZombie()) {
    std::cerr << "Failed to open input file: " << inputFile << std::endl;
    return;
  }

  TString treePath = TString::Format("%s/%s", directoryName, treeName);
  TTree *tree = static_cast<TTree *>(input.Get(treePath));
  if (!tree) {
    std::cerr << "Failed to find tree: " << treePath << std::endl;
    return;
  }

  unsigned int fFilterMap = 0;
  bool fIsCorrectAssoc = false;
  LegInfo leg1;
  LegInfo leg2;

  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("fFilterMap", 1);
  tree->SetBranchStatus("fIsCorrectAssoc", 1);
  tree->SetBranchStatus("fPin1", 1);
  tree->SetBranchStatus("fEta1", 1);
  tree->SetBranchStatus("fTPCNSigmaEl1", 1);
  tree->SetBranchStatus("fTPCNSigmaPi1", 1);
  tree->SetBranchStatus("fTPCNSigmaPr1", 1);
  tree->SetBranchStatus("fPin2", 1);
  tree->SetBranchStatus("fEta2", 1);
  tree->SetBranchStatus("fTPCNSigmaEl2", 1);
  tree->SetBranchStatus("fTPCNSigmaPi2", 1);
  tree->SetBranchStatus("fTPCNSigmaPr2", 1);

  tree->SetBranchAddress("fFilterMap", &fFilterMap);
  tree->SetBranchAddress("fIsCorrectAssoc", &fIsCorrectAssoc);
  tree->SetBranchAddress("fPin1", &leg1.p);
  tree->SetBranchAddress("fEta1", &leg1.eta);
  tree->SetBranchAddress("fTPCNSigmaEl1", &leg1.nSigmaEl);
  tree->SetBranchAddress("fTPCNSigmaPi1", &leg1.nSigmaPi);
  tree->SetBranchAddress("fTPCNSigmaPr1", &leg1.nSigmaPr);
  tree->SetBranchAddress("fPin2", &leg2.p);
  tree->SetBranchAddress("fEta2", &leg2.eta);
  tree->SetBranchAddress("fTPCNSigmaEl2", &leg2.nSigmaEl);
  tree->SetBranchAddress("fTPCNSigmaPi2", &leg2.nSigmaPi);
  tree->SetBranchAddress("fTPCNSigmaPr2", &leg2.nSigmaPr);

  std::vector<double> pBins;
  std::vector<double> etaBins;
  if (!LoadReferenceBinning(referenceBinningFile, referenceBinningHist, pBins, etaBins)) {
    pBins = nPBins > 0 ? BuildUniformBins(pMin, pMax, nPBins) : BuildDefaultPinBins();
    etaBins = nEtaBins > 0 ? BuildUniformBins(etaMin, etaMax, nEtaBins) : BuildUniformBins(etaMin, etaMax, 16);
  }
  const int nPinBinsToUse = static_cast<int>(pBins.size()) - 1;
  const int nEtaBinsToUse = static_cast<int>(etaBins.size()) - 1;

  TH2D hDenAvg("hDenPidEffAvg", "Average PID denominator;p_{in} (GeV/#it{c});#eta", nPinBinsToUse, pBins.data(), nEtaBinsToUse, etaBins.data());
  TH2D hNumAvg("hNumPidEffAvg", "Average PID numerator;p_{in} (GeV/#it{c});#eta", nPinBinsToUse, pBins.data(), nEtaBinsToUse, etaBins.data());
  TH2D hPidEffAvg("hPidEffAvg", "Average leg PID efficiency;p_{in} (GeV/#it{c});#eta", nPinBinsToUse, pBins.data(), nEtaBinsToUse, etaBins.data());

  TH2D hDenIfOtherPass("hDenPidEffIfOtherPass", "PID denominator if other leg passes;p_{in} (GeV/#it{c});#eta", nPinBinsToUse, pBins.data(), nEtaBinsToUse, etaBins.data());
  TH2D hNumIfOtherPass("hNumPidEffIfOtherPass", "PID numerator if other leg passes;p_{in} (GeV/#it{c});#eta", nPinBinsToUse, pBins.data(), nEtaBinsToUse, etaBins.data());
  TH2D hPidEffIfOtherPass("hPidEffIfOtherPass", "Leg PID efficiency if other leg passes;p_{in} (GeV/#it{c});#eta", nPinBinsToUse, pBins.data(), nEtaBinsToUse, etaBins.data());

  TH2D hDenIfOtherFail("hDenPidEffIfOtherFail", "PID denominator if other leg fails;p_{in} (GeV/#it{c});#eta", nPinBinsToUse, pBins.data(), nEtaBinsToUse, etaBins.data());
  TH2D hNumIfOtherFail("hNumPidEffIfOtherFail", "PID numerator if other leg fails;p_{in} (GeV/#it{c});#eta", nPinBinsToUse, pBins.data(), nEtaBinsToUse, etaBins.data());
  TH2D hPidEffIfOtherFail("hPidEffIfOtherFail", "Leg PID efficiency if other leg fails;p_{in} (GeV/#it{c});#eta", nPinBinsToUse, pBins.data(), nEtaBinsToUse, etaBins.data());

  const Long64_t entries = tree->GetEntries();
  const Long64_t entriesToRead = maxEntries < 0 ? entries : std::min(maxEntries, entries);
  Long64_t selectedEvents = 0;
  Long64_t selectedLegs = 0;
  Long64_t pidPassedLegs = 0;

  for (Long64_t i = 0; i < entriesToRead; ++i) {
    tree->GetEntry(i);

    if (!PassBit(fFilterMap, filterMapBit)) {
      continue;
    }
    if (requireCorrectAssoc && !fIsCorrectAssoc) {
      continue;
    }

    const bool leg1PassPid = PassPid(leg1, nSigmaElMin, nSigmaElMax, nSigmaPiMin, nSigmaPrMin);
    const bool leg2PassPid = PassPid(leg2, nSigmaElMin, nSigmaElMax, nSigmaPiMin, nSigmaPrMin);

    FillLeg(leg1, leg1PassPid, leg2PassPid,
            hDenAvg, hNumAvg,
            hDenIfOtherPass, hNumIfOtherPass,
            hDenIfOtherFail, hNumIfOtherFail);
    FillLeg(leg2, leg2PassPid, leg1PassPid,
            hDenAvg, hNumAvg,
            hDenIfOtherPass, hNumIfOtherPass,
            hDenIfOtherFail, hNumIfOtherFail);

    ++selectedEvents;
    selectedLegs += 2;
    pidPassedLegs += static_cast<Long64_t>(leg1PassPid) + static_cast<Long64_t>(leg2PassPid);
  }

  MakeEfficiency(hPidEffAvg, hNumAvg, hDenAvg);
  MakeEfficiency(hPidEffIfOtherPass, hNumIfOtherPass, hDenIfOtherPass);
  MakeEfficiency(hPidEffIfOtherFail, hNumIfOtherFail, hDenIfOtherFail);

  TFile output(outputFile, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create output file: " << outputFile << std::endl;
    return;
  }

  hDenAvg.Write();
  hNumAvg.Write();
  hPidEffAvg.Write();
  hDenIfOtherPass.Write();
  hNumIfOtherPass.Write();
  hPidEffIfOtherPass.Write();
  hDenIfOtherFail.Write();
  hNumIfOtherFail.Write();
  hPidEffIfOtherFail.Write();
  output.Close();

  std::cout << "Read entries: " << entriesToRead << std::endl;
  std::cout << "Selected events: " << selectedEvents << std::endl;
  std::cout << "Selected legs: " << selectedLegs << std::endl;
  std::cout << "PID-passed legs: " << pidPassedLegs << std::endl;
  std::cout << "Wrote PID maps to: " << outputFile << std::endl;
}
