/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


/** \class MuXboostedBTagger
 *
 *  MuXboostedBTagger tags high-pT heavy flavor jets using a muonic tag
 *  (see arXiv:1511.xxxxx for a full discussion of the underlying physics).
 *  Check for updates @ <https://github.com/keith-pedersen/delphes/tree/MuXboostedBTag>
 *
 *
 *  One of the crucial features of MuXboostedBTagger is that it
 *  ALTERS the pT of tagged jets, which come from semi-leptonic decay
 *  by construction, and thus are always missing neutrino energy.
 *
 *  MuXboostedBTagger accounts for this missing energy by ESTIMATING
 *  the neutrino's energy (currently using the simplest choice pNu=pMu,
 *  from the  shared boost). Thus, "taggable" muons are added a second
 *  time to jets that are tagged.
 *
 *  Therefore, in addition to altering the BTag bit of the original jet,
 *  (for which no neutrino estimation is performed), MuXboostedBTagger
 *  also clones the list of jets (adding neutrino energy when they
 *  pass the tag).
 *
 *
 *
 *  NOTE: This module gives more accurate light jet fake rates when
 *  used in conjunction with AllParticlePropagator (which simulates the
 *  initial bending of in-flight Pion/Kaon decays to muons).
 *
 *
 *
 *  The algorithm
 *
 *  1. Clone each jet into the output array, regardless of whether it's tagged.
 *
 *  2. Look for jets passing a preliminary pT cut (no less than half
 *     the final cut, since neutrino estimation could potentially
 *     double a jet's pT)
 *       2a. (ptJet > 0.5*fMinJetPt)
 *
 *  3. Look inside the jet for at least one "taggable" muon (muon pT >= fMinMuonPt).
 *     WARNING: This requires muons to be clustered into jets during
 *     the initial jet clustering; MuXboostedBTagger does not take a
 *     muon input array.
 *       3a. (ptMu > fMinMuonPt)
 *       3b. After all "taggable" neutrinos are found, ensure that the
 *           final jet pT (with neutrino estimation) will pass the pT cut
 *           (ptJet [with neutrinos] > fMinJetPt).
 *
 *  4. Reconstruct the subjet of the semi-leptonic decay
 *       4a. Recluster the jet (discarding towers whose pT ratio, to the
 *           original jet pT, is less than fMinTowerPtRatio), using
 *           anti-kt with radius fCoreAntiktR.
 *       4b. The mass of the core is poorly measured, constrain it
 *           to fCoreMassHypothesis.
 *       4c. Search through the list of candidate cores and, if their
 *           boost is greater than fCoreMinBoost, calculate the mass
 *           of the resulting subjet (using only the hardest muon).
 *       4d. The "correct" core is the one which produces mSubjet closest
 *           to fSubjetMassHypothesis.
 *
 *  5. Calculate x
 *       5a. x = tan(thetaLab)*Esubjet/mSubjet. However, restrict mSubjet
 *           to be no larger than fMaxSubjetMass (in case of a
 *           poorly reconstructed subjet).
 *
 *  6. Calculate fSubjet
 *       6a. fSubjet = pTsubjet/pTjet (account for neutrino estimation).
 *
 *  7. Tag the jet
 *       7a. If (x <= fMaxX) and (f >= fSubjet), tag both the original jet and
 *           the clone (using fBitNumber in the BTag field).
 *       7b. If the jet is tagged, add the neutrino pT to the clone.
 *
 *  \author K. Pedersen - Illinois Institute of Technology
 *
 */

#include "modules/MuXboostedBTagger.h"
#include "classes/DelphesClasses.h"

#include "TObjArray.h"
#include "TLorentzVector.h"
//~ #include "TDirectory.h"
//~ #include "TFile.h"
//~ #include "TH1I.h"
//~ #include "TH1F.h"
//~ #include "TH2F.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>
//#include <iostream>
//#include <iomanip>
//#include <sstream>

using namespace std;

const Double_t PI = acos(-1.);

//------------------------------------------------------------------------------

MuXboostedBTagger::MuXboostedBTagger():
   fAllParticles(0), fCoreDefinition(0), fJetInputArray(0), fItJetInputArray(0), fJetOutputArray(0)
{}

//------------------------------------------------------------------------------

MuXboostedBTagger::~MuXboostedBTagger()
{ /* Memory clean-up handled in Finish() */}

//------------------------------------------------------------------------------

void MuXboostedBTagger::Init()
{
   // read parameters

   fBitNumber = GetInt("BitNumber", 3);
   fTaggedBit = 1<<fBitNumber;

   fMaxX = GetDouble("MaxX", 3.);
   fMinCorePtRatio = GetDouble("fMinCorePtRatio", .5);
   
   fMinJetPt = GetDouble("MinJetPt", 300.);
   fMinMuonPt = GetDouble("MinMuonPt", 10.);

   fMinTowerPtRatio = GetDouble("MinTowerPtRatio", .05);
   fCoreAntiktR = GetDouble("CoreAntiktR", .04);
   
   fCoreMinBoost = GetDouble("CoreMinBoost", 1.);
   fCoreMinBoostSquared = Squared(fCoreMinBoost);

   fCoreMassHypothesis = GetDouble("CoreMassHypothesis", 2.0);
   fCoreMassHypothesisSquared = Squared(fCoreMassHypothesis);

   fSubjetMassHypothesis = GetDouble("SubjetMassHypothesis", 5.3);

   fMaxSubjetMass = GetDouble("MaxSubjetMass", 12.);
   
   fAllParticles = ImportArray("Delphes/allParticles");
   
   fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
   fItJetInputArray = fJetInputArray->MakeIterator();
   
   fJetOutputArray = ExportArray(GetString("JetOutputArray", "jets"));
   
   // Create the core definition
   fCoreDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm, fCoreAntiktR);
}

//------------------------------------------------------------------------------

void MuXboostedBTagger::Finish()
{
   delete fCoreDefinition;
   delete fItJetInputArray;
}

//------------------------------------------------------------------------------

void MuXboostedBTagger::Process()
{
   Candidate* originalJet;
   Candidate* cloneJet;
   
   fItJetInputArray->Reset();

   // loop over input jets
   while((originalJet = static_cast<Candidate*>(fItJetInputArray->Next())))
   {
      cloneJet = static_cast<Candidate*>(originalJet->Clone());
      fJetOutputArray->Add(cloneJet);

      const Double_t originalJetPt = originalJet->Momentum.Pt();

      // At best, neutrino estimation can only double the pt of the jet,
      // so only look at jets with at least half the required final pt
      if(originalJetPt >= .5*fMinJetPt)
      {
         // We'll sort jet constituents into two categories 
         // (where taggableMuons pass the pt cut)
         std::vector<Candidate*> taggableMuons;
         std::vector<Candidate const*> everythingElse;
         everythingElse.reserve(32);
         
         Double_t neutrinoPt = 0.;

         // Fill taggableMuons & everythingElse 
         {
            TIter itJetConstituents(cloneJet->GetCandidates());
            Candidate* constituent;

            // Find all muons that pass the pt cuts
            while((constituent = static_cast<Candidate*>(itJetConstituents.Next())))
            {
               if(abs(constituent->PID) == 13)
               {
                  const Double_t muonPt = (constituent->Momentum).Pt();
                  if(muonPt >= fMinMuonPt)
                  {
                     neutrinoPt += muonPt;
                     taggableMuons.push_back(constituent);
                     continue;
                  }
               }
               everythingElse.push_back(constituent);
            }
         }
         
         // Ensure we have at least one taggable muon and possibly enough pT after neutrino estimation
         if((not taggableMuons.empty()) and (originalJetPt + neutrinoPt >= fMinJetPt))
         {
            // Sort muons low to high
            {
               // Sometimes there will be more than 1 muon. In those cases, we will
               // assume that the muons were emitted from highest to lowest pt.
               // This means we'll need to add them back from low to high pt, so that
               // the CM frame of the muon emitted later doesn't include the previously emitted muon.

               // Sort muons low to high
               std::sort(taggableMuons.begin(), taggableMuons.end(), SortCandidatePt_Low2High);
            }

            // Code for finding muon mother/matriarch
            // Store PID inside muon in unused Candidate fields
            // Useful for understanding tag, but useless without 
            // way of getting info out of Delphes 
            // (new output classes, altered TreeWriter).
            
            /*
            std::vector<Candidate const*> taggableMuonsMatriarch;
            // Find the muon matriarch and mother 
            // Store inside each muon
            //       Muon.BTag = abs(PID_mother)
            //       Muon.TauTag = abs(PID_matriarch)
            //       Muon.Tau[2] = xTrue (real x to matriarch)
            {
               // * The "matriarch" is the original boosted primary hadron
               //      We need to look back all the way to the start of hadronization.
               // * The "mother" is the particle which emitted the muon
               //      If the mother is a tau, we consider the grandmother the mother (because we are interested in hadron flavor)
               for(auto itMuon = taggableMuons.begin(); itMuon != taggableMuons.end(); ++itMuon)
               {
                  Candidate const* matriarch;

                  // We need to account for PileUp, which may not have fully de-refernceable mothers
                  if((*itMuon)->IsPU == 1)
                  {
                     matriarch = 0; // Pile-up is assumed to be non-derferenceable, so it has no matriarch
                     // Use default values for BTag, TauTag, and Tau[2]
                  }
                  else
                  {
                     bool noMotherYet = true;
                     // To look back all the way to hadronization, we find the first matriarch with
                     // two seperate mothers (more than 1 mother means color charge still existed)
                     Candidate const* possibleMatriarch = static_cast<Candidate*>(fAllParticles->At((*itMuon)->M1));

                     if(not possibleMatriarch)
                        matriarch = 0;
                     else
                     {
                        while(possibleMatriarch)
                        {
                           const UInt_t absPossibleMatriarchPID = abs(possibleMatriarch->PID);
                           matriarch = possibleMatriarch;

                           // Tau cannot be mother, we want hadrons
                           if(noMotherYet and (absPossibleMatriarchPID not_eq 15))
                           {
                              noMotherYet = false;
                              (*itMuon)->BTag = absPossibleMatriarchPID;
                           }

                           // Check to see if we've found the matriarch
                           if( (not ((possibleMatriarch->M2 == 0) or (possibleMatriarch->M2 == possibleMatriarch->M1))) or
                              ((absPossibleMatriarchPID >= 22) and (absPossibleMatriarchPID <= 24)))
                           {
                              // First check (post hadronic)
                              // Post hadronic decays should have a M2 == 0 or (in the case of interaction
                              // with the vacuum) M1 == M2. Anthing else and we're reaching back into color charge

                              // Second check (primary muons from A, Z, W)
                              // Pythia hadronic decay never specifically invokes a W or Z
                              // so any encountered here must have come from the hard interaction.
                              // Pythia DOES produce A, which can split to muons, but not during decay (only during fragmentation)

                              // We break here because we've found the matriarch.
                              // If we break here, the mother should already be set.
                              break; // out of the do loop
                           }
                           else // Keep tracing back the lineage
                           {
                              possibleMatriarch = static_cast<Candidate*>(fAllParticles->At(possibleMatriarch->M1));
                           }
                        }

                        // Store the matriarch's PID
                        (*itMuon)->TauTag = abs(matriarch->PID);

                        if(noMotherYet)
                        {
                           // If I understand Pythia correctly, we should never reach here
                           // So let me know if we ever reach here, then default the mother to the matriarch
                           cout << "\n\nMotherFail\n\n" << endl;
                           (*itMuon)->BTag = (*itMuon)->TauTag;
                        }

                        // Calculate the muon's true x to the matriarch
                        // Store in Tau[2] (repurposed N-sub-jettyness field)
                        const TVector3 matriarchP3 = (matriarch->Momentum).Vect();
                        const TVector3 muonP3 = ((*itMuon)->Momentum).Vect();

                        (*itMuon)->Tau[2] = (matriarch->Momentum).E() *
                           (muonP3.Cross(matriarchP3)).Mag() /
                           ((matriarch->Momentum).M() * muonP3.Dot(matriarchP3));
                     }
                  }
                
                  // Record the muon's matriarch
                  taggableMuonsMatriarch.push_back(matriarch);
               }
            }
            */

            // Now we find the jet's core by reculstering its jet constituents
            vector<fastjet::PseudoJet> reclusterInput;

            // Add the taggable muons to the reclusterInput
            {
               // Give them a negative user_index ( user_index = index in taggableMuons - taggableMuons.size() )
               // This makes them easy to find after the reclustering.

               const int numMuons = taggableMuons.size();
               for(int iMu = 0; iMu < numMuons; ++iMu)
               {
                  const TLorentzVector& muonMomentum = taggableMuons[iMu]->Momentum;
                  reclusterInput.emplace_back(muonMomentum.Px(), muonMomentum.Py(), muonMomentum.Pz(), muonMomentum.E());
                  reclusterInput.back().set_user_index(iMu - numMuons);
               }
            }

            // Add jet constituents (all tracks, but only towers passing pt cut)
            {
               // Use all tracks, but only use towers/eFlowNeutrals (charge == 0) if they are above the relative pt threshold
               // This is because the angular resolution of towers is much poorer, and we don't
               // want pileup / soft QCD unduly influencing the direction of the core

               // Here, use original pT, not neutrino estimated pT
               const Double_t minTowerPt = fMinTowerPtRatio*originalJetPt;
               for(unsigned int iEverythingElse = 0; iEverythingElse < everythingElse.size(); ++iEverythingElse)
               {
                  Candidate const* const constituent = everythingElse[iEverythingElse];

                  // neutral charge ==> tower (an assumption)
                  if((constituent->Charge == 0) && ((constituent->Momentum).Pt() < minTowerPt))
                     continue;

                  reclusterInput.emplace_back((constituent->Momentum).Px(), (constituent->Momentum).Py(), (constituent->Momentum).Pz(), (constituent->Momentum).E());
                  reclusterInput.back().set_user_index(iEverythingElse);
               }
            }

            // Scope of fastjet ClusterSequence
            {
               // Recluster the jet, to find core candidates
               fastjet::ClusterSequence recluster(reclusterInput, *fCoreDefinition);

               // No minimum pT for reclustered candidates (for now)
               const Double_t minCoreCandidatePt = 0.;
               
               // Get the core candidates (NOT sorted by pT until we remove muons)
               std::vector<fastjet::PseudoJet> coreCandidates = recluster.inclusive_jets(minCoreCandidatePt);

               // Make sure at least one subjet passed the cut
               if(not coreCandidates.empty())
               {
                  fastjet::PseudoJet core;
                  TLorentzVector p4core;
                  
                  // Get fastjets internal jets (because fastjet only has const access to our input jets,
                  // so they have no internal clustering information, which we'll need to use).
                  // The documentation does not specify that the input jets are the first
                  // members of list returned by jets(), so we'll have to assume that
                  // this is the case (this assumption was tested a lot during debug)
                  const std::vector<fastjet::PseudoJet>& internalJets = recluster.jets();

                  // Remove all taggable muons from the core candidates, then sort by pt
                  {
                     const int numMuons = taggableMuons.size();
                     for(int iMu = 0; iMu < numMuons; ++iMu)
                     {
                        const fastjet::PseudoJet& muon = internalJets[iMu];
                        // Debugging code
                        if(muon.user_index() not_eq (iMu - numMuons))
                           throw runtime_error("MuXboostedBTagger: Muon ID Fail!");

                        for(auto itCoreCandidate = coreCandidates.begin(); itCoreCandidate not_eq coreCandidates.end(); ++itCoreCandidate)
                        {
                           if(muon.is_inside(*itCoreCandidate))
                           {
                              *itCoreCandidate -= muon;
                              break; // out of inner loop - the muon can only be in one candidate at a time
                           }
                        }
                     }
                     
                     coreCandidates = fastjet::sorted_by_pt(coreCandidates);
                  }

                  // Find the core 
                  // WARNING: only the hardest muon is used to find the core
                  {
                     // Here me must find the "correct" core by finding the
                     // subjet mass closest to fSubjetMassHypothesis 
                     // (using only the hardest muon)
                     // 
                     // Please refer to /delphes/doc/MuXboostedBTagger.pdf/Code notes/1
                     // for an explanation of the mass math (i.e. g and y)

                     const fastjet::PseudoJet& hardestMuon = reclusterInput[taggableMuons.size()-1];
                     
                     Double_t minDeltaMass = 9e9;
                     auto itCore = coreCandidates.begin(); // Default to hardest core

                     for(auto itCoreCandidate = coreCandidates.begin(); itCoreCandidate not_eq coreCandidates.end(); ++itCoreCandidate)
                     {
                        const Double_t g = Squared(fCoreMassHypothesis / itCoreCandidate->E());
                        // Make sure we still have sufficient boost after muon subtraction
                        if(g*fCoreMinBoostSquared > 1.)
                           break; // Since we're sorted by pt, the next core candidate is also too weak

                        const Double_t y = Tan2(*itCoreCandidate, hardestMuon);

                        // The core has the mass closest to fSubjetMassHypothesis
                        const Double_t deltaMass = abs(sqrt(fCoreMassHypothesisSquared + 
                           (4. * hardestMuon.E() * itCoreCandidate->E() * (g + y)) / (1. + y + sqrt(1. - ((g - y) + g*y)))) - fSubjetMassHypothesis);

                        if(deltaMass < minDeltaMass)
                        {
                           minDeltaMass = deltaMass;
                           itCore = itCoreCandidate;
                        }
                     }
                     
                     core = *itCore;
                  }

                  // From now on, we will keep the original PseudoJet core around,
                  // in case we are interested in constituents. However, all 4-vector math
                  // will now be handled by p4core (i.e. core will no longer have an 
                  // accurate momentum).
                  
                  // Fix core mass to fCoreMassHypothesis
                  {                    
                     // The current mass depends more on the granularity of the CAL more than anything else
                     // To fix the mass, We need to scale the momentum of the core (but not the energy)
                     const Double_t momentumScale = sqrt(((core.E() - fCoreMassHypothesis)*(core.E() + fCoreMassHypothesis))/core.modp2());
                     p4core.SetPxPyPzE(momentumScale*core.px(), momentumScale*core.py(), momentumScale*core.pz(), core.E());
                  }

                  // Reconstruct the subjet, measure x, find the smallest value
                  {
                     // Keep track of the smallest x for any muon
                     Double_t minX = 1024.; // Large enough, simple significand
                     
                     // Iterate through each muon, add it to the core, and calculate the resulting boost invariant
                     for(auto itMuon = taggableMuons.begin(); itMuon not_eq taggableMuons.end(); ++itMuon)
                     {
                        const TLorentzVector& muonP4 = (*itMuon)->Momentum;
                        const TVector3 muonP3 = muonP4.Vect();
                        
                        // Neutrino estimation
                        // A more sophisticated choice will improve dijet mass resolution
                        const TLorentzVector& neutrinoP4 = muonP4;

                        // Add the muon and the neutrino
                        p4core += muonP4;
                        p4core += neutrinoP4;
                        
                        // Now calcuale x
                        const Double_t boostMass = std::min(p4core.M(), fMaxSubjetMass);
                        const Double_t xCore = (p4core.E() * (muonP3.Cross(p4core.Vect())).Mag()) / (muonP3.Dot(p4core.Vect()) * boostMass);

                        // Now add the neutrino to original jet ...
                        // IFF the muon passes the boosted test
                        if(xCore <= fMaxX)
                           (cloneJet->Momentum) += neutrinoP4;

                        minX = std::min(minX, xCore);
                     }// End loop over muons
                     
                     const Double_t jetPt = (cloneJet->Momentum).Pt();

                     // Re-check pt now that we've fully estimated neutrinos
                     if(jetPt >= fMinJetPt)
                     {
                        // Ensure we had a boosted muon emission from
                        // inside a subjet which appears to be a b-hadron
                        if((minX <= fMaxX) and 
                           ((p4core.Pt()/jetPt) >= fMinCorePtRatio))
                        {
                           // Tag the original jet
                           originalJet->BTag |= fTaggedBit;
                           // and the neutrino estimated jet                           
                           cloneJet->BTag |= fTaggedBit;                           
                        }// End tagged
                     }// End final jet pT check
                  }// End subjet reconstruction 
               }// End valid subjets
            }// End reclustering to find core candidates
         }//End muon pt cut and secondary jet pT cut
      }//End jet initial pt cut
   }//End loop over jets
}

//------------------------------------------------------------------------------

// Returns PID of the heaviest quark (for hadrons), 
// otherwise returns 0 (i.e. quarks,leptons, bosons).
// This is not totally perfect, but works for most normal SM particles
Int_t HadronFlavor(Int_t PID)
{
   PID = abs(PID);

   if((PID > 1000) and (PID < 10000)) // Baryons
      PID /= 1000;
   else // Mesons/diquarks/everything else
      PID = (PID%1000)/100;

   return PID;
}

//------------------------------------------------------------------------------

bool SortCandidatePt_Low2High(Candidate const* one, Candidate const* two)
{
   return (one->Momentum).Pt() < (two->Momentum).Pt();
}


//------------------------------------------------------------------------------

inline Double_t Squared(const Double_t arg)
{
   return arg*arg;
}

//------------------------------------------------------------------------------

// Find the tan(theta)**2 between two PseudoJets by calculating (p3>_1 x p3>_2)**2 / (p3>_1 . p3>_2)**2
Double_t Tan2(const fastjet::PseudoJet& one, const fastjet::PseudoJet& two)
{
   // This is more accurate than an alternate form 
   //    (one.one*two.two-(one.two)**2)/(one.two)**2 
   // because the cancellation in the numerator (cross-product) is 
   // handled component-by-component, rather than all at once.
   Double_t
      crossSquared =  Squared(one.px()*two.py() - one.py()*two.px());
      crossSquared += Squared(one.px()*two.pz() - one.pz()*two.px());
      crossSquared += Squared(one.py()*two.pz() - one.pz()*two.py());

   const Double_t dot = one.px()*two.px() + one.py()*two.py() + one.pz()*two.pz();

   return crossSquared / (dot*dot);
}

//------------------------------------------------------------------------------

Double_t VectorAngle(const TVector3& one, const TVector3& two)
{
   // TVector3 finds the angle by acos(one.two / sqrt(one.one*two.two))
   // This has a relative error of 0.5*epsilon / angle**2 (where epsilon is machine epsilon)
   // due to catastrophic cancellation

   // This form is better (from my personal experiments)
   return atan2((one.Cross(two)).Mag(), one.Dot(two));
}
