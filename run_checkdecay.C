#include <cmath>
#include <iostream>

#include "TFile.h"
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

void FillEfficiencyHistograms(const CheckDecay::JpsiDielectron &event,
                              const CheckDecay::PidEfficiencyMaps &pidMaps,
                              TH2D &hGenerated,
                              TH2D &hPassedWeighted,
                              TH2D &hClosureDiff,
                              TH2D &hClosureCount)
{
  const double pt = event.jpsi.Pt();
  const double y = event.jpsi.Rapidity();
  const double passProb = pidMaps.BothPassProbability(event.leg1, event.leg2);
  const double closure12 = pidMaps.LegAverageClosureDifference(event.leg1, event.leg2);
  const double closure21 = pidMaps.LegAverageClosureDifference(event.leg2, event.leg1);

  hGenerated.Fill(pt, y);
  hPassedWeighted.Fill(pt, y, passProb);
  hClosureDiff.Fill(pt, y, 0.5 * (closure12 + closure21));
  hClosureCount.Fill(pt, y);
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
                    const char *failMapName = "hPidEffIfOtherFail")
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
  TH2D hClosureDiff("hPidMapClosureDiffSum", "Sum of conditional-map closure differences;p_{T} (GeV/#it{c});y", nPtBins, ptMin, ptMax, nYBins, yMin, yMax);
  TH2D hClosureCount("hPidMapClosureCount", "Events used for conditional-map closure;p_{T} (GeV/#it{c});y", nPtBins, ptMin, ptMax, nYBins, yMin, yMax);
  TH2D hClosureMean("hPidMapClosureDiffMean", "Mean conditional-map closure difference;p_{T} (GeV/#it{c});y", nPtBins, ptMin, ptMax, nYBins, yMin, yMax);

  TString sourceMode(mode);
  sourceMode.ToLower();

  Long64_t processed = 0;
  if (sourceMode == "simulate") {
    TRandom3 rng(seed);
    for (Long64_t i = 0; i < nJpsi; ++i) {
      const CheckDecay::JpsiDielectron event =
        CheckDecay::GenerateUniformJpsiDecay(rng, ptMin, ptMax, yMin, yMax);
      FillEfficiencyHistograms(event, pidMaps, hGenerated, hPassedWeighted, hClosureDiff, hClosureCount);
    }
    processed = nJpsi;
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

    const char *requiredLeaves[] = {
      "jpsi_pt", "jpsi_eta", "jpsi_phi", "jpsi_energy",
      "leg1_pt", "leg1_eta", "leg1_phi",
      "leg2_pt", "leg2_eta", "leg2_phi"
    };
    for (const char *leafName : requiredLeaves) {
      if (!HasLeaf(tree, leafName)) {
        std::cerr << "Missing required branch/leaf: " << leafName << std::endl;
        return;
      }
    }

    const Long64_t entriesToRead = std::min(nJpsi, tree->GetEntries());
    for (Long64_t i = 0; i < entriesToRead; ++i) {
      tree->GetEntry(i);
      CheckDecay::JpsiDielectron event;
      event.jpsi = CheckDecay::MakeParticleFromPtEtaPhiE(LeafValue(tree, "jpsi_pt"),
                                                          LeafValue(tree, "jpsi_eta"),
                                                          LeafValue(tree, "jpsi_phi"),
                                                          LeafValue(tree, "jpsi_energy"));
      event.leg1 = CheckDecay::MakeParticleFromPtEtaPhiM(LeafValue(tree, "leg1_pt"),
                                                          LeafValue(tree, "leg1_eta"),
                                                          LeafValue(tree, "leg1_phi"),
                                                          CheckDecay::kMassElectron);
      event.leg2 = CheckDecay::MakeParticleFromPtEtaPhiM(LeafValue(tree, "leg2_pt"),
                                                          LeafValue(tree, "leg2_eta"),
                                                          LeafValue(tree, "leg2_phi"),
                                                          CheckDecay::kMassElectron);
      FillEfficiencyHistograms(event, pidMaps, hGenerated, hPassedWeighted, hClosureDiff, hClosureCount);
    }
    processed = entriesToRead;
  } else {
    std::cerr << "Unknown mode: " << mode << ". Use simulate or tree." << std::endl;
    return;
  }

  FinalizeEfficiency(hEfficiency, hPassedWeighted, hGenerated);
  FinalizeAverage(hClosureMean, hClosureDiff, hClosureCount);

  TFile output(outputFile, "RECREATE");
  if (output.IsZombie()) {
    std::cerr << "Failed to create output file: " << outputFile << std::endl;
    return;
  }

  hGenerated.Write();
  hPassedWeighted.Write();
  hEfficiency.Write();
  hClosureDiff.Write();
  hClosureCount.Write();
  hClosureMean.Write();
  output.Close();

  std::cout << "Processed " << processed << " J/psi candidates in mode " << sourceMode << std::endl;
  std::cout << "Wrote output file: " << outputFile << std::endl;
}
