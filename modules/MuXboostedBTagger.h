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

#ifndef MuXboostedBTagger_h
#define MuXboostedBTagger_h

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

#include "classes/DelphesModule.h"

class TObjArray;
class TIterator;
class Candidate;
class TVector3;
//class TFile;

namespace fastjet
{
   class JetDefinition;
   class PseudoJet;
}

class MuXboostedBTagger: public DelphesModule
{
   public:
      MuXboostedBTagger();
      ~MuXboostedBTagger();

      void Init();
      void Process();
      void Finish();

   private:
      Int_t fBitNumber; // Candidate->BTag bit for tag (0 = default bit used by BTagging)
      Int_t fTaggedBit;
      
      //Main tagging parameters
      Double_t fMaxX; // xMax for the tag
      Double_t fMinCorePtRatio; // fMin for the tag
      
      // pT minimums
      Double_t fMinJetPt;     // Min jet pT to apply for tag
      Double_t fMinMuonPt;    // Min muon pT to be a "taggable" muon
      
      // Parameters for subjet reconstruction
      Double_t fMinTowerPtRatio; // Min pT ratio of tower to entire jet, to qualify for reclustering. Helps reduce pileup sensitivity  
      Double_t fCoreAntiktR; // Reculstering radius
      
      Double_t fCoreMinBoost; // Minimum core boost during core selection (keep this low to reduce fake rate)
      Double_t fCoreMinBoostSquared; // Squared version, to prevent redundant math
      
      Double_t fCoreMassHypothesis; // Mass constraint applied to core
      Double_t fCoreMassHypothesisSquared; // Squared version, to prevent redundant math

      Double_t fSubjetMassHypothesis; // Target subjet mass for core selection

      Double_t fMaxSubjetMass; // Ceiling of reconstructed subjet mass, for sanity

      // List of entire event, for searching for mother/matriarch (currently commented out)
      TObjArray const* fAllParticles;
      
      fastjet::JetDefinition* fCoreDefinition;

      TObjArray const* fJetInputArray; //!
      TIterator* fItJetInputArray; //!
      
      DelphesFactory* fFactory;
      
      TObjArray* fJetOutputArray; //!
      
      ClassDef(MuXboostedBTagger, 1)
};

// Helper functions
Int_t HadronFlavor(Int_t PID);
inline Double_t Squared(const Double_t arg);
Double_t Tan2(const fastjet::PseudoJet& one, const fastjet::PseudoJet& two);
Double_t VectorAngle(const TVector3& one, const TVector3& two);
bool SortCandidatePt_Low2High(Candidate const* one, Candidate const* two);

#endif
