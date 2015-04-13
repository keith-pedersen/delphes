#include "KDPClasses.h"

CompBase TaggingEfficiencyJet::fgCompare = CompPT<TaggingEfficiencyJet>::Instance();
CompBase TaggingEfficiencyMuon::fgCompare = CompPT<TaggingEfficiencyMuon>::Instance();

const Float_t TaggingEfficiencyMuon::Mass = 0.1056583715;

TLorentzVector TaggingEfficiencyJet::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

TLorentzVector TaggingEfficiencyMuon::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}
