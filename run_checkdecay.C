#include <algorithm>
#include <cmath>
#include <iostream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLeaf.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"

#include "JpsiDecay.C"
#include "PidEfficiencyMaps.C"

namespace {

double LeafValue(TTree *tree, const char *leafName)
{
  if (!tree || !leafName) {
    return 0.0;
  }
  TLeaf *leaf = tree->GetLeaf(leafName);
  if (!leaf) {
    return 0.0;
  }
  return leaf->GetValue();
}

bool HasLeaf(TTree *tree, const char *leafName)
{
  return tree && leafName && tree->GetLeaf(leafName);
}

bool PassBitSelection(TTree *tree, const char *leafName, int bit)
{
  if (bit < 0) {
    return true;
  }
  if (bit > 31) {
    std::cerr << "Bit index out of range for " << leafName << ": " << bit << std::endl;
    return false;
  }
  const unsigned int value = static_cast<unsigned int>(LeafValue(tree, leafName));
  return (value & (1u << bit)) != 0u;
}

bool PassBoolSelection(TTree *tree, const char *leafName, bool requiredValue)
{
  return (LeafValue(tree, leafName) != 0.0) == requiredValue;
}

bool CheckRequiredLeaves(TTree *tree, const char *const *requiredLeaves, int nLeaves)
{
  for (int i = 0; i < nLeaves; ++i) {
    if (!HasLeaf(tree, requiredLeaves[i])) {
      std::cerr << "Missing required branch/leaf: " << requiredLeaves[i] << std::endl;
      return false;
    }
  }
  return true;
}

void FillEfficiencyHistograms(const CheckDecay::JpsiDielectron &event,
                              const CheckDecay::PidEfficiencyMaps &pidMaps,
                              TH2D &hGenerated,
                              TH2D &hPassedWeighted,
                              TH1D &hGeneratedPt,
                              TH1D &hPassedWeightedPt,
                              TH2D &hClosureDiff,
                              TH2D &hClosureCount,
                              double eventWeight = 1.0)
{
  const double pt = event.jpsi.Pt();
  const double y = event.jpsi.Rapidity();
  const double passProb = pidMaps.BothPassProbability(event.leg1, event.leg2);
  const double closure12 = pidMaps.LegAverageClosureDifference(event.leg1, event.leg2);
  const double closure21 = pidMaps.LegAverageClosureDifference(event.leg2, event.leg1);

  hGenerated.Fill(pt, y, eventWeight);
  hPassedWeighted.Fill(pt, y, eventWeight * passProb);
  hGeneratedPt.Fill(pt, eventWeight);
  hPassedWeightedPt.Fill(pt, eventWeight * passProb);
  hClosureDiff.Fill(pt, y, eventWeight * 0.5 * (closure12 + closure21));
  hClosureCount.Fill(pt, y, eventWeight);
}

void FinalizeEfficiency(TH2D &hEfficiency, const TH2D &hPassedWeighted, const TH2D &hGenerated)
{
  hEfficiency.Reset("ICES");
  for (int ix = 1; ix <= hEfficiency.GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hEfficiency.GetNbinsY(); ++iy) {
      const double generated = hGenerated.GetBinContent(ix, iy);
      if (generated <= 0.0) {
        continue;
      }
      const double passed = hPassedWeighted.GetBinContent(ix, iy);
      hEfficiency.SetBinContent(ix, iy, passed / generated);
      hEfficiency.SetBinError(ix, iy, 0.0);
    }
  }
}

void FinalizeEfficiency(TH1D &hEfficiency, const TH1D &hPassedWeighted, const TH1D &hGenerated)
{
  hEfficiency.Reset("ICES");
  for (int ix = 1; ix <= hEfficiency.GetNbinsX(); ++ix) {
    const double generated = hGenerated.GetBinContent(ix);
    if (generated <= 0.0) {
      continue;
    }
    const double passed = hPassedWeighted.GetBinContent(ix);
    hEfficiency.SetBinContent(ix, passed / generated);
    hEfficiency.SetBinError(ix, 0.0);
  }
}

void FinalizeAverage(TH2D &hAverage, const TH2D &hSum, const TH2D &hCount)
{
  hAverage.Reset("ICES");
  for (int ix = 1; ix <= hAverage.GetNbinsX(); ++ix) {
    for (int iy = 1; iy <= hAverage.GetNbinsY(); ++iy) {
      const double count = hCount.GetBinContent(ix, iy);
      if (count <= 0.0) {
        continue;
      }
      hAverage.SetBinContent(ix, iy, hSum.GetBinContent(ix, iy) / count);
      hAverage.SetBinError(ix, iy, 0.0);
    }
  }
}

double HistogramWeight(const TH2D *hist, double x, double y)
{
  if (!hist) {
    return 1.0;
  }
  const int xbin = hist->GetXaxis()->FindFixBin(x);
  const int ybin = hist->GetYaxis()->FindFixBin(y);
  if (xbin < 1 || xbin > hist->GetNbinsX() || ybin < 1 || ybin > hist->GetNbinsY()) {
    return 0.0;
  }
  return hist->GetBinContent(xbin, ybin);
}

bool PassSimAcceptance(const CheckDecay::JpsiDielectron &event,
                       double jpsiPtMin,
                       double jpsiPtMax,
                       double jpsiYMin,
                       double jpsiYMax,
                       double legPtMin,
                       double legAbsEtaMax)
{
  const double jpsiPt = event.jpsi.Pt();
  const double jpsiY = event.jpsi.Rapidity();
  if (jpsiPt < jpsiPtMin || jpsiPt >= jpsiPtMax || jpsiY < jpsiYMin || jpsiY > jpsiYMax) {
    return false;
  }
  if (event.leg1.Pt() <= legPtMin || std::abs(event.leg1.Eta()) >= legAbsEtaMax) {
    return false;
  }
  if (event.leg2.Pt() <= legPtMin || std::abs(event.leg2.Eta()) >= legAbsEtaMax) {
    return false;
  }
  return true;
}

TH2D *BuildRecJpsiPtYWeightHistogram(const char *treeFile,
                                     const char *treeName,
                                     int nPtBins,
                                     double ptMin,
                                     double ptMax,
                                     int nYBins,
                                     double yMin,
                                     double yMax,
                                     int mcDecisionBit,
                                     int filterMapBit,
                                     bool requireCorrectAssoc)
{
  TFile input(treeFile, "READ");
  if (input.IsZombie()) {
    std::cerr << "Failed to open weight tree input file: " << treeFile << std::endl;
    return nullptr;
  }
  TTree *tree = static_cast<TTree *>(input.Get(treeName));
  if (!tree) {
    std::cerr << "Failed to find weight tree " << treeName << " in " << treeFile << std::endl;
    return nullptr;
  }

  const char *requiredLeaves[] = {
    "fMass", "fPt", "fEta", "fPhi", "fMcDecision", "fFilterMap", "fIsCorrectAssoc"
  };
  if (!CheckRequiredLeaves(tree, requiredLeaves, sizeof(requiredLeaves) / sizeof(requiredLeaves[0]))) {
    return nullptr;
  }

  TH2D *hist = new TH2D("hInputJpsiPtYWeight", "Input rec J/#psi p_{T}-y weight;p_{T} (GeV/#it{c});y",
                        nPtBins, ptMin, ptMax, nYBins, yMin, yMax);
  hist->SetDirectory(nullptr);

  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    if (!PassBitSelection(tree, "fMcDecision", mcDecisionBit)) {
      continue;
    }
    if (!PassBitSelection(tree, "fFilterMap", filterMapBit)) {
      continue;
    }
    if (requireCorrectAssoc && !PassBoolSelection(tree, "fIsCorrectAssoc", true)) {
      continue;
    }
    const TLorentzVector jpsi =
      CheckDecay::MakeParticleFromPtEtaPhiM(LeafValue(tree, "fPt"),
                                            LeafValue(tree, "fEta"),
                                            LeafValue(tree, "fPhi"),
                                            LeafValue(tree, "fMass"));
    hist->Fill(jpsi.Pt(), jpsi.Rapidity());
  }

  return hist;
}

} // namespace

void run_checkdecay(const char *mode = "simulate",
                    Long64_t nJpsi = 1000000,
                    const char *pidMapFile = "input/pid_efficiency_maps.root",
                    const char *outputFile = "output/jpsi_pid_efficiency.root",
                    int nPtBins = 80,
                    double ptMin = 0.0,
                    double ptMax = 20.0,
                    int nYBins = 80,
                    double yMin = -4.0,
                    double yMax = 4.0,
                    const char *treeFile = "input/jpsi_to_ee_tree.root",
                    const char *treeName = "jpsiTree",
                    UInt_t seed = 12345,
                    double fallbackAvgEff = 0.85,
                    double fallbackIfOtherPassEff = 0.90,
                    double fallbackIfOtherFailEff = 0.70,
                    const char *avgMapName = "hPidEffAvg",
                    const char *passMapName = "hPidEffIfOtherPass",
                    const char *failMapName = "hPidEffIfOtherFail",
                    const char *treeLevel = "gen",
                    int mcDecisionBit = 0,
                    int filterMapBit = 0,
                    bool requireCorrectAssoc = true,
                    bool useInputRecPtYWeightForSimulation = false,
                    const char *weightTreeFile = "input/LHC25b16_mcrecgen.root",
                    const char *weightTreeName = "DF_2363808317891232/O2rtdiemtreerec",
                    int weightMcDecisionBit = 0,
                    int weightFilterMapBit = 0,
                    bool weightRequireCorrectAssoc = true,
                    bool requireSimAcceptance = true,
                    double simAccJpsiPtMin = 0.0,
                    double simAccJpsiPtMax = 20.0,
                    double simAccJpsiYMin = -0.8,
                    double simAccJpsiYMax = 0.8,
                    double simAccLegPtMin = 1.0,
                    double simAccLegAbsEtaMax = 0.8)
{
  if (nJpsi <= 0) {
    std::cerr << "nJpsi must be positive." << std::endl;
    return;
  }

  CheckDecay::PidEfficiencyMaps pidMaps;
  pidMaps.SetFallbackEfficiencies(fallbackAvgEff, fallbackIfOtherPassEff, fallbackIfOtherFailEff);
  pidMaps.Load(pidMapFile, avgMapName, passMapName, failMapName);

  TH2D hGenerated("hJpsiGenerated", "Generated J/#psi;p_{T} (GeV/#it{c});y", nPtBins, ptMin, ptMax, nYBins, yMin, yMax);
  TH2D hPassedWeighted("hJpsiPidPassedWeighted", "Weighted PID-passed J/#psi;p_{T} (GeV/#it{c});y", nPtBins, ptMin, ptMax, nYBins, yMin, yMax);
  TH2D hEfficiency("hJpsiPidEfficiency", "J/#psi PID efficiency;p_{T} (GeV/#it{c});y", nPtBins, ptMin, ptMax, nYBins, yMin, yMax);
  const double ptBinsFor1D[] = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0};
  const int nPtBinsFor1D = sizeof(ptBinsFor1D) / sizeof(ptBinsFor1D[0]) - 1;
  TH1D hGeneratedPt("hJpsiGeneratedPt", "Generated J/#psi;p_{T} (GeV/#it{c});counts", nPtBinsFor1D, ptBinsFor1D);
  TH1D hPassedWeightedPt("hJpsiPidPassedWeightedPt", "Weighted PID-passed J/#psi;p_{T} (GeV/#it{c});weighted counts", nPtBinsFor1D, ptBinsFor1D);
  TH1D hEfficiencyPt("hJpsiPidEfficiencyPt", "J/#psi PID efficiency;p_{T} (GeV/#it{c});efficiency", nPtBinsFor1D, ptBinsFor1D);
  TH2D hClosureDiff("hPidMapClosureDiffSum", "Sum of conditional-map closure differences;p_{T} (GeV/#it{c});y", nPtBins, ptMin, ptMax, nYBins, yMin, yMax);
  TH2D hClosureCount("hPidMapClosureCount", "Events used for conditional-map closure;p_{T} (GeV/#it{c});y", nPtBins, ptMin, ptMax, nYBins, yMin, yMax);
  TH2D hClosureMean("hPidMapClosureDiffMean", "Mean conditional-map closure difference;p_{T} (GeV/#it{c});y", nPtBins, ptMin, ptMax, nYBins, yMin, yMax);

  TString sourceMode(mode);
  sourceMode.ToLower();

  Long64_t processed = 0;
  Long64_t selected = 0;
  if (sourceMode == "simulate") {
    TH2D *hInputJpsiPtYWeight = nullptr;
    if (useInputRecPtYWeightForSimulation) {
      hInputJpsiPtYWeight =
        BuildRecJpsiPtYWeightHistogram(weightTreeFile, weightTreeName,
                                       nPtBins, ptMin, ptMax, nYBins, yMin, yMax,
                                       weightMcDecisionBit, weightFilterMapBit, weightRequireCorrectAssoc);
      if (!hInputJpsiPtYWeight) {
        return;
      }
    }

    TRandom3 rng(seed);
    for (Long64_t i = 0; i < nJpsi; ++i) {
      const CheckDecay::JpsiDielectron event =
        CheckDecay::GenerateUniformJpsiDecay(rng, ptMin, ptMax, yMin, yMax);
      if (requireSimAcceptance &&
          !PassSimAcceptance(event,
                             simAccJpsiPtMin,
                             simAccJpsiPtMax,
                             simAccJpsiYMin,
                             simAccJpsiYMax,
                             simAccLegPtMin,
                             simAccLegAbsEtaMax)) {
        continue;
      }
      const double eventWeight = HistogramWeight(hInputJpsiPtYWeight, event.jpsi.Pt(), event.jpsi.Rapidity());
      FillEfficiencyHistograms(event, pidMaps, hGenerated, hPassedWeighted, hGeneratedPt, hPassedWeightedPt, hClosureDiff, hClosureCount, eventWeight);
      ++selected;
    }
    processed = nJpsi;
    if (hInputJpsiPtYWeight) {
      delete hInputJpsiPtYWeight;
      hInputJpsiPtYWeight = nullptr;
    }
  } else if (sourceMode == "tree") {
    TFile input(treeFile, "READ");
    if (input.IsZombie()) {
      std::cerr << "Failed to open tree input file: " << treeFile << std::endl;
      return;
    }
    TTree *tree = static_cast<TTree *>(input.Get(treeName));
    if (!tree) {
      std::cerr << "Failed to find tree " << treeName << " in " << treeFile << std::endl;
      return;
    }

    TString level(treeLevel);
    level.ToLower();

    const char *requiredGenLeaves[] = {
      "fE", "fPt", "fEta", "fPhi",
      "fPtMC1", "fEtaMC1", "fPhiMC1",
      "fPtMC2", "fEtaMC2", "fPhiMC2",
      "fMcDecision"
    };
    const char *requiredRecLeaves[] = {
      "fMass", "fPt", "fEta", "fPhi",
      "fPt1", "fEta1", "fPhi1",
      "fPt2", "fEta2", "fPhi2",
      "fMcDecision", "fFilterMap", "fIsCorrectAssoc"
    };

    if (level == "gen") {
      if (!CheckRequiredLeaves(tree, requiredGenLeaves, sizeof(requiredGenLeaves) / sizeof(requiredGenLeaves[0]))) {
        return;
      }
    } else if (level == "rec") {
      if (!CheckRequiredLeaves(tree, requiredRecLeaves, sizeof(requiredRecLeaves) / sizeof(requiredRecLeaves[0]))) {
        return;
      }
    } else {
      std::cerr << "Unknown treeLevel: " << treeLevel << ". Use gen or rec." << std::endl;
      return;
    }

    const Long64_t entriesToRead = std::min(nJpsi, tree->GetEntries());
    for (Long64_t i = 0; i < entriesToRead; ++i) {
      tree->GetEntry(i);
      CheckDecay::JpsiDielectron event;

      if (!PassBitSelection(tree, "fMcDecision", mcDecisionBit)) {
        continue;
      }

      if (level == "gen") {
        event.jpsi = CheckDecay::MakeParticleFromPtEtaPhiE(LeafValue(tree, "fPt"),
                                                            LeafValue(tree, "fEta"),
                                                            LeafValue(tree, "fPhi"),
                                                            LeafValue(tree, "fE"));
        event.leg1 = CheckDecay::MakeParticleFromPtEtaPhiM(LeafValue(tree, "fPtMC1"),
                                                            LeafValue(tree, "fEtaMC1"),
                                                            LeafValue(tree, "fPhiMC1"),
                                                            CheckDecay::kMassElectron);
        event.leg2 = CheckDecay::MakeParticleFromPtEtaPhiM(LeafValue(tree, "fPtMC2"),
                                                            LeafValue(tree, "fEtaMC2"),
                                                            LeafValue(tree, "fPhiMC2"),
                                                            CheckDecay::kMassElectron);
      } else {
        if (!PassBitSelection(tree, "fFilterMap", filterMapBit)) {
          continue;
        }
        if (requireCorrectAssoc && !PassBoolSelection(tree, "fIsCorrectAssoc", true)) {
          continue;
        }

        event.jpsi = CheckDecay::MakeParticleFromPtEtaPhiM(LeafValue(tree, "fPt"),
                                                            LeafValue(tree, "fEta"),
                                                            LeafValue(tree, "fPhi"),
                                                            LeafValue(tree, "fMass"));
        event.leg1 = CheckDecay::MakeParticleFromPtEtaPhiM(LeafValue(tree, "fPt1"),
                                                            LeafValue(tree, "fEta1"),
                                                            LeafValue(tree, "fPhi1"),
                                                            CheckDecay::kMassElectron);
        event.leg2 = CheckDecay::MakeParticleFromPtEtaPhiM(LeafValue(tree, "fPt2"),
                                                            LeafValue(tree, "fEta2"),
                                                            LeafValue(tree, "fPhi2"),
                                                            CheckDecay::kMassElectron);
      }

      FillEfficiencyHistograms(event, pidMaps, hGenerated, hPassedWeighted, hGeneratedPt, hPassedWeightedPt, hClosureDiff, hClosureCount);
      ++selected;
    }
    processed = entriesToRead;
  } else {
    std::cerr << "Unknown mode: " << mode << ". Use simulate or tree." << std::endl;
    return;
  }

  FinalizeEfficiency(hEfficiency, hPassedWeighted, hGenerated);
  FinalizeEfficiency(hEfficiencyPt, hPassedWeightedPt, hGeneratedPt);
  FinalizeAverage(hClosureMean, hClosureDiff, hClosureCount);

  TFile output(outputFile, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create output file: " << outputFile << std::endl;
    return;
  }

  hGenerated.Write();
  hPassedWeighted.Write();
  hEfficiency.Write();
  hGeneratedPt.Write();
  hPassedWeightedPt.Write();
  hEfficiencyPt.Write();
  hClosureDiff.Write();
  hClosureCount.Write();
  hClosureMean.Write();
  output.Close();

  std::cout << "Processed " << processed << " J/psi candidates in mode " << sourceMode << std::endl;
  std::cout << "Selected " << selected << " J/psi candidates for efficiency filling" << std::endl;
  std::cout << "Wrote output file: " << outputFile << std::endl;
}
