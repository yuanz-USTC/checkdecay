#ifndef CHECKDECAY_JPSIDECAY_C
#define CHECKDECAY_JPSIDECAY_C

#include <cmath>

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TVector3.h"

namespace CheckDecay {

constexpr double kMassJpsi = 3.0969;
constexpr double kMassElectron = 0.00051099895;

struct JpsiDielectron {
  TLorentzVector jpsi;
  TLorentzVector leg1;
  TLorentzVector leg2;
};

TLorentzVector MakeJpsiFromPtYPhi(double pt, double y, double phi)
{
  const double mt = std::sqrt(pt * pt + kMassJpsi * kMassJpsi);
  const double px = pt * std::cos(phi);
  const double py = pt * std::sin(phi);
  const double pz = mt * std::sinh(y);
  const double e = mt * std::cosh(y);
  return TLorentzVector(px, py, pz, e);
}

TLorentzVector MakeParticleFromPtEtaPhiM(double pt, double eta, double phi, double mass)
{
  const double px = pt * std::cos(phi);
  const double py = pt * std::sin(phi);
  const double pz = pt * std::sinh(eta);
  const double p2 = px * px + py * py + pz * pz;
  const double e = std::sqrt(p2 + mass * mass);
  return TLorentzVector(px, py, pz, e);
}

TLorentzVector MakeParticleFromPtEtaPhiE(double pt, double eta, double phi, double energy)
{
  const double px = pt * std::cos(phi);
  const double py = pt * std::sin(phi);
  const double pz = pt * std::sinh(eta);
  return TLorentzVector(px, py, pz, energy);
}

JpsiDielectron DecayJpsiToDielectron(const TLorentzVector &jpsi, TRandom3 &rng)
{
  const double eRest = kMassJpsi / 2.0;
  const double pRest = std::sqrt(std::max(0.0, eRest * eRest - kMassElectron * kMassElectron));
  const double cosTheta = rng.Uniform(-1.0, 1.0);
  const double sinTheta = std::sqrt(std::max(0.0, 1.0 - cosTheta * cosTheta));
  const double phi = rng.Uniform(0.0, 2.0 * M_PI);

  TLorentzVector leg1(pRest * sinTheta * std::cos(phi),
                      pRest * sinTheta * std::sin(phi),
                      pRest * cosTheta,
                      eRest);
  TLorentzVector leg2(-leg1.Px(), -leg1.Py(), -leg1.Pz(), eRest);

  const TVector3 boost = jpsi.BoostVector();
  leg1.Boost(boost);
  leg2.Boost(boost);

  return {jpsi, leg1, leg2};
}

JpsiDielectron GenerateUniformJpsiDecay(TRandom3 &rng,
                                        double ptMin,
                                        double ptMax,
                                        double yMin,
                                        double yMax)
{
  const double pt = rng.Uniform(ptMin, ptMax);
  const double y = rng.Uniform(yMin, yMax);
  const double phi = rng.Uniform(0.0, 2.0 * M_PI);
  return DecayJpsiToDielectron(MakeJpsiFromPtYPhi(pt, y, phi), rng);
}

} // namespace CheckDecay

#endif
