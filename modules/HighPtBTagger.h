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

#ifndef HighPtBTagger_h
#define HighPtBTagger_h

/** \class HighPtBTagger
 *
 *  Tags high-pt jets with a heavy flavor tag.
 *
 *  1. Look for jets passing a pt cut.
 *       1a.    (ptJet > fMinJetPt)
 *  2. Look inside for at least one muon passing a pt cut. This helps ensure
 *     that the the muon is well reconstructed.
 *       2a.    (ptMu > fMinMuonPt)
 *  3. Look for a hard central core carrying a majority of the jet's pt
 *       3a. Cluster the jet's constituents using anti-kt with (R = fCoreAntiktR)
 *       3b. Require "the core" (highest pt sub-jet) to have (pt_subjet / pt_jet >= fMinCoreRatio)
 *           and also still contain the muon.
 *  4. Adjust the core's 4-vector to account for the unobserved neutrino
 *       3a. Add the muon again to estimate the neutrino, BUT only project out the
 *           muon component parallel to the core's axis, so as not to bias the
 *           direction of the jet, and keep the neutrino massless, so the jet mass
 *           doesn't change.
 *  5. Find the mass of the core, given a B quark hypothesis
 *       5a.   coreMass = max(fMinCoreMass, mass(pMu_core))
 *  6. Find the boost (gamma) of the core, using its mass
 *       6a. This "mass" essentially comes from an angular spread of energy in the
 *           calorimeter, and depends only on the hard central components, so it
 *           won't be sensitive to the fragmentation/hadronization treatment
 *           of wide angle components.
 *       6b. Unfortunately, in the highly boosted environment, there is good way
 *           to discriminate between (B -> mu nuMu X C) and (B -> X1 C -> X1 X2 mu nuMu).
 *           Obviously, the muon from the secondary c-decay will be from a less boosted
 *           environment, and will have a slightly higher emission angle. But we can
 *           only calculate one approximate gamma.
 *  7. Find the emission angle "invariant"
 *       7a. Assume z^ is the boost axis of the jet (before emission), (theta) is the
 *           angle the muon makes with the boost axis in the jet's CM frame, and
 *           (theta') is the angle the muon makes with the jet in the lab frame.
 *       7b. The "invariant" (assuming the muon is effectively massless in the CM frame
 *           and gamma >> 1) is
 *               invariant = gamma * tan(theta') = tan(theta/2)
 *  8. If (fSigma = 1 - epsilon**2) is the % of the total solid angle of emission
 *     you want to capture, then (thetaMax = Pi - 2*epsilon).
 *       8a. Thus tag any particle who's (invariant <= tan(thetaMax/2)), which
 *           is equivalent to (invariant <= cot(sqrt(1-fSigma)))
 *       8b. This filters out wide angle components (i.e. from light pi decays in
 *           flight) quite well.
 *
 *  \author K. Pedersen - Illinois Institute of Technology
 *
 */

#include "classes/DelphesModule.h"

class TObjArray;
class TIterator;
class Candidate;
class TFile;
class TH1I;
class TH1F;
class TH2F;

namespace fastjet
{
  class JetDefinition;
}

class HighPtBTagger: public DelphesModule
{
	public:

		HighPtBTagger();
		~HighPtBTagger();

		void Init();
		void Process();
		void Finish();

	private:
		void FillTH1F(const std::vector<std::vector<TH1F*> >& histoVec, const Int_t coreIndex, const Int_t flavor, const Double_t fillValue);
		void FillTH2F(const std::vector<std::vector<TH2F*> >& histoVec, const Int_t coreIndex, const Int_t flavor, const Double_t fillValueX, const Double_t fillValueY);

		void FillMotherPIDHisto(const Int_t PID);
		void WriteMotherPIDHisto();

		Int_t fBitNumber;
		Double_t fMinJetPt;
		Double_t fMinMuonPt;
		Double_t fMinTowerPtRatio;
		Double_t fCoreAntiktR;
		Double_t fMinCorePtRatio;
		Double_t fCoreMassHypothesis;
		Double_t fCoreMassHypothesis2;
		Double_t fMinFinalMass;
		Double_t fMaxFinalMass;
		Double_t fMaxEmissionInvariant;

		TObjArray const* fJetInputArray; //!
		TIterator* fItJetInputArray; //!

		TObjArray const* fAllParticles;

		fastjet::JetDefinition* fCoreDefinition;

		TDirectory* histoDir;

		std::map<Int_t, Int_t> pdgID;

		TH1I
			*pt_Jets,
			*eta_Jets,
			*pt_Muons,
			*motherPID;

		// To make histograms more versatile, we will create vectors of them

		// These histograms are not core dependent
		//     x_True[motherFlavor]
		std::vector<TH1F*>
			x_True,
			pt_MuJets,
			eta_MuJets;

		// These histograms depend on core reconstruction, and the core's matriarch
		//     histoVec[coreIndex][matriarchFlavor]
		std::vector<std::vector<TH1F*> >
			coreRatio,
			mass;

		// These histograms depend on core reconstruction, and the muons mother
		//     histoVec[coreIndex][matriarchFlavor]
		std::vector<std::vector<TH1F*> >
			deltaMass,
			x_Core,
			deltaTrue,
			deltaTrue2Mu,
			pt_Tagged,
			eta_Tagged;

		std::vector<std::vector<TH2F*> >
			x_Core__vs__x_True,
			x_Core__vs__coreRatio;

		// We'll need a map from flavor to histogram index
		std::map<Int_t, Int_t> flavorToHistoIndex;

		ClassDef(HighPtBTagger, 1)
};

bool SortCandidatePt_Low2High(Candidate const* one, Candidate const* two);

#endif
