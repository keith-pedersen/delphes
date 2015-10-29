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


/** \class HighPtBTagger
 *
 *  Tags high-pt jets with a heavy flavor tag
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/BoostedBTagger.h"
#include "classes/DelphesClasses.h"

#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

const Double_t PI = acos(-1.);

// Returns heaviest quark PID for hadrons, otherwise returns actual PID
UInt_t BoostedBTagger::Flavor(Int_t PID)
{
	PID = abs(PID);

	if(PID > 1000 && PID < 10000)
		PID /= 1000;
	else
		PID = (PID%1000)/100;

	return PID;
}

BoostedBTagger::BoostedBTagger():
	fJetInputArray(0), fItJetInputArray(0), fAllParticles(0)
{}

//------------------------------------------------------------------------------

BoostedBTagger::~BoostedBTagger()
{ /* No memory to clean-up */}

//------------------------------------------------------------------------------

inline Double_t Squared(const Double_t arg)
{
	return arg*arg;
}

//~ // Find the tan(theta)**2 between two PseudoJets by calculating (p3>_1 x p3>_2)**2 / (p3>_1 . p3>_2)**2
//~ Double_t Tan2(const fastjet::PseudoJet& one, const fastjet::PseudoJet& two)
//~ {
	//~ Double_t
		//~ crossSquared =  Squared(one.px()*two.py() - one.py()*two.px());
		//~ crossSquared += Squared(one.px()*two.pz() - one.pz()*two.px());
		//~ crossSquared += Squared(one.py()*two.pz() - one.pz()*two.py());

	//~ const Double_t dot = one.px()*two.px() + one.py()*two.py() + one.pz()*two.pz();

	//~ return crossSquared / (dot*dot);
//~ }

//~ Double_t AccurateAngle(const TVector3& one, const TVector3& two)
//~ {
	//~ // TVector3 finds the angle by acos(one.two / sqrt(one.one*two.two))
	//~ // This has a relative error of 0.5*epsilon / angle**2 (where epsilon is machine epsilon)
	//~ // due to catastrophic cancellation

	//~ // This form is better (from my personal experiments)
	//~ const Double_t dot = one.Dot(two);

	//~ return atan2(sqrt((one.Cross(two)).Mag2()/(dot*dot)), std::copysign(1., dot));
//~ }

void BoostedBTagger::Init()
{
 	// read parameters

	fBitNumber = GetInt("BitNumber", 4);

	fMaxJetRank = GetInt("MaxJetRank", 4);

	fMinJetPt = GetDouble("MinJetPt", 100.);
	fMinMuonPt = GetDouble("MinMuonPt", 20.);
        
        fMaxDeltaR = GetDouble("MaxDeltaR", 0.1);

	// import input array(s)
	fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
	fItJetInputArray = fJetInputArray->MakeIterator();
        
        fAllParticles = ImportArray("Delphes/allParticles");
}



//------------------------------------------------------------------------------

void BoostedBTagger::Finish()
{
        delete fItJetInputArray;
}

//------------------------------------------------------------------------------

// After processing:
//     Jet with:
//        good muons (potential for tagging) ........... TauTag == numMuons
//        viable HardCore .............................. FracPt[0] > 0.
//        viable BCore (which isn't the HardCore) .... FracPt[1] > 0.
//     Muon:
//        good muon in jet with interesting pt ......... BTag/TauTag non-zero
//        good muon in jet with interesting core ....... Tau[2] > 0. / Tau[0] > 0. / FracPt[0] < 8.
//        good muon in jet with distinct BCore ....... Tau[1] > 0. / FracPt[1] < 8.

void BoostedBTagger::Process()
{
	// loop over input jets
	Candidate* jet;
	fItJetInputArray->Reset();

        // Here we assume that the jet's have already been sorted by pT (by FastJetFinder, which calls sorted_by_pt on the clustered list)
	Int_t jetRank = 0;

	while((jet = static_cast<Candidate*>(fItJetInputArray->Next())) and (++jetRank <= fMaxJetRank))
	{
		const TLorentzVector& jetMomentum = jet->Momentum; 
                const Double_t jetPt = jetMomentum.Pt();
                
                // Set used fields to NULL values (checked 23.06.2015)
                {
			jet->FracPt[0] = -1.;   // f_subjet of HardCore
			jet->FracPt[1] = -1.;   // f_subjet of BCore
			jet->Tau[0] = -1;       // x_min for HardCore
			jet->Tau[1] = -1;       // x_min for BCore
                        jet->TauTag = 0;        // Maximum flavor of jet (previously Number of Good Muons)
			// Do not clear BTag, since we only set one bit
		}

		// At best, neutrino estimation can only double the pt of the jet,
		// so only look at jets with at least half the required final pt
		if(jetPt >= fMinJetPt)
		{
			// We'll sort jet constituents into two categories (where goodMuons are those which pass the pt cut)
			std::vector<Candidate*> goodMuons;
			{
				TIter itJetConstituents(jet->GetCandidates());
				Candidate* constituent;

				// Find all muons that pass the pt cuts
				while((constituent = static_cast<Candidate*>(itJetConstituents.Next())))
				{
					if(abs(constituent->PID) == 13)
					{
						if((constituent->Momentum).Pt() >= fMinMuonPt)
						{
                                                        // Set used fields to NULL values;
							constituent->Tau[0] = -1.;      // x_HardCore
							constituent->Tau[1] = -1.;      // x_BCore
							constituent->Tau[2] = -1.;      // x_True
							constituent->FracPt[0] = 8.;    // deltaAngle(recoHardCore, matriarch) from adding muon twice
							constituent->FracPt[1] = 8.;    // deltaAngle(recoBCore, matriarch) from adding muon twice
							constituent->BTag = 0;          // abs(PID_mother)
							constituent->TauTag = 0;        // abs(PID_matriarch)
                                                        
							goodMuons.push_back(constituent);
						}
					}
				}
			}

			if(not goodMuons.empty()) // Ensure we have at least one good muon
			{
				// Find the muon matriarch and mother (CHECKED 23.06.2015)
				// Store inside each muon
				//       Muon.BTag = abs(PID_mother)
				//       Muon.TauTag = abs(PID_matriarch)
				//       Muon.Tau[2] = xTrue (real x to matriarch)
				{
					// * The "matriarch" is the original boosted primary hadron
					//      We need to look back all the way to the start of hadronization.
					// * The "mother" is the particle which emitted the muon
					//      If the mother is a tau, we consider the grandmother the mother (because we are interested in hadron flavor)
					for(auto itMuon = goodMuons.begin(); itMuon != goodMuons.end(); ++itMuon)
					{
						Candidate const* matriarch;

						// We need to account for PileUp, which may not have fully derefernceable mothers
						if((*itMuon)->IsPU == 0)
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
						else
						{
							matriarch = 0; // Pile-up is assumed to be non-derferenceable, so it has no matriarch
							// Use default values for BTag, TauTag, and Tau[2]
						}
                                                
                                                if(matriarch not_eq 0)
                                                        jet->TauTag = std::max(jet->TauTag, Flavor(matriarch->PID));
					}
				}
                                
                                for(auto itMuon = goodMuons.begin(); itMuon != goodMuons.end(); ++itMuon)
                                {
                                        if(jetMomentum.DeltaR((*itMuon)->Momentum) <= fMaxDeltaR)
                                        {
                                                jet->BTag |= (1 << fBitNumber);
                                        }
                                }
			}//End muon pt cut
		}//End jet initial pt cut
	}//End loop over jets
}
