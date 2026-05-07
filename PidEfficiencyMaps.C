#ifndef CHECKDECAY_PIDEFFICIENCYMAPS_C
#define CHECKDECAY_PIDEFFICIENCYMAPS_C

#include <algorithm>
#include <iostream>

#include "TAxis.h"
#include "TFile.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TObject.h"
#include "TString.h"

namespace CheckDecay {

class PidEfficiencyMaps {
public:
  PidEfficiencyMaps() = default;

  ~PidEfficiencyMaps()
  {
    if (m_file) {
      m_file->Close();
      delete m_file;
      m_file = nullptr;
    }
  }

  bool Load(const char *fileName,
            const char *avgMapName = "hPidEffAvg",
            const char *passMapName = "hPidEffIfOtherPass",
            const char *failMapName = "hPidEffIfOtherFail")
  {
    if (!fileName || TString(fileName).IsNull()) {
      std::cout << "No PID map file configured. Using fallback efficiencies." << std::endl;
      return false;
    }

    m_file = TFile::Open(fileName, "READ");
    if (!m_file || m_file->IsZombie()) {
      std::cout << "Cannot open PID map file " << fileName
                << ". Using fallback efficiencies." << std::endl;
      return false;
    }

    m_avg = ReadTH2(avgMapName);
    m_ifOtherPass = ReadTH2(passMapName);
    m_ifOtherFail = ReadTH2(failMapName);

    if (!m_avg || !m_ifOtherPass || !m_ifOtherFail) {
      std::cout << "At least one PID map is missing. Missing maps use fallback efficiencies." << std::endl;
    }
    return m_avg || m_ifOtherPass || m_ifOtherFail;
  }

  void SetFallbackEfficiencies(double avg, double ifOtherPass, double ifOtherFail)
  {
    m_fallbackAvg = ClampEfficiency(avg);
    m_fallbackIfOtherPass = ClampEfficiency(ifOtherPass);
    m_fallbackIfOtherFail = ClampEfficiency(ifOtherFail);
  }

  double Average(const TLorentzVector &leg) const
  {
    return Value(m_avg, leg.P(), leg.Eta(), m_fallbackAvg);
  }

  double IfOtherPass(const TLorentzVector &leg) const
  {
    return Value(m_ifOtherPass, leg.P(), leg.Eta(), m_fallbackIfOtherPass);
  }

  double IfOtherFail(const TLorentzVector &leg) const
  {
    return Value(m_ifOtherFail, leg.P(), leg.Eta(), m_fallbackIfOtherFail);
  }

  double BothPassProbability(const TLorentzVector &leg1, const TLorentzVector &leg2) const
  {
    const double order12 = Average(leg1) * IfOtherPass(leg2);
    const double order21 = Average(leg2) * IfOtherPass(leg1);
    return ClampEfficiency(0.5 * (order12 + order21));
  }

  double LegAverageClosureDifference(const TLorentzVector &firstLeg, const TLorentzVector &secondLeg) const
  {
    const double firstAvg = Average(firstLeg);
    const double predictedSecondAvg =
      firstAvg * IfOtherPass(secondLeg) + (1.0 - firstAvg) * IfOtherFail(secondLeg);
    return predictedSecondAvg - Average(secondLeg);
  }

private:
  TH2 *ReadTH2(const char *name)
  {
    if (!m_file || !name || TString(name).IsNull()) {
      return nullptr;
    }
    TObject *obj = m_file->Get(name);
    if (!obj || !obj->InheritsFrom(TH2::Class())) {
      std::cout << "PID map " << name << " not found as TH2." << std::endl;
      return nullptr;
    }
    return static_cast<TH2 *>(obj);
  }

  static double ClampEfficiency(double eff)
  {
    if (!std::isfinite(eff)) {
      return 0.0;
    }
    return std::max(0.0, std::min(1.0, eff));
  }

  static int FindClampedBin(const TAxis *axis, double value)
  {
    if (!axis) {
      return 1;
    }
    const int nbins = axis->GetNbins();
    if (value < axis->GetXmin()) {
      return 1;
    }
    if (value >= axis->GetXmax()) {
      return nbins;
    }
    return axis->FindBin(value);
  }

  static double Value(const TH2 *hist, double p, double eta, double fallback)
  {
    if (!hist) {
      return fallback;
    }
    const int xbin = FindClampedBin(hist->GetXaxis(), p);
    const int ybin = FindClampedBin(hist->GetYaxis(), eta);
    const double eff = hist->GetBinContent(xbin, ybin);
    return ClampEfficiency(eff);
  }

  TFile *m_file = nullptr;
  TH2 *m_avg = nullptr;
  TH2 *m_ifOtherPass = nullptr;
  TH2 *m_ifOtherFail = nullptr;
  double m_fallbackAvg = 0.85;
  double m_fallbackIfOtherPass = 0.90;
  double m_fallbackIfOtherFail = 0.70;
};

} // namespace CheckDecay

#endif
