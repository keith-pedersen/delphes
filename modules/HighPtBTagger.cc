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

#include "modules/HighPtBTagger.h"
#include "classes/DelphesClasses.h"

#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

const Double_t PI = acos(-1.);

// Returns heaviest quark PID for hadrons, otherwise returns actual PID
Int_t flavor(Int_t PID)
{
	PID = abs(PID);

	if(PID > 1000 && PID < 10000)
		PID /= 1000;
	else
		PID = (PID%1000)/100;

	return PID;
}

//------------------------------------------------------------------------------

HighPtBTagger::HighPtBTagger():
	fJetInputArray(0), fItJetInputArray(0), fAllParticles(0), fCoreDefinition(0)
{}

//------------------------------------------------------------------------------

HighPtBTagger::~HighPtBTagger()
{ /* No memory to clean-up */}

//------------------------------------------------------------------------------

inline Double_t Squared(const Double_t arg)
{
	return arg*arg;
}

// Find the tan(theta)**2 between two PseudoJets by calculating (p3>_1 x p3>_2)**2 / (p3>_1 . p3>_2)**2
Double_t Tan2(const fastjet::PseudoJet& one, const fastjet::PseudoJet& two)
{
	Double_t
		crossSquared =  Squared(one.px()*two.py() - one.py()*two.px());
		crossSquared += Squared(one.px()*two.pz() - one.pz()*two.px());
		crossSquared += Squared(one.py()*two.pz() - one.pz()*two.py());

	const Double_t dot = one.px()*two.px() + one.py()*two.py() + one.pz()*two.pz();

	return crossSquared / (dot*dot);
}

Double_t AccurateAngle(const TVector3& one, const TVector3& two)
{
	// TVector3 finds the angle by acos(one.two / sqrt(one.one*two.two))
	// This has a relative error of 0.5*epsilon / angle**2 (where epsilon is machine epsilon)
	// due to catastrophic cancellation

	// A form suggested by WK (in "How Mindless"),  accurate over all domains
	const TVector3
		v1 = one*two.Mag(),
		v2 = two*one.Mag();

	return 2*atan(sqrt((v1-v2).Mag2()/(v1+v2).Mag2()));
}

void HighPtBTagger::Init()
{
 	// read parameters

	fBitNumber = GetInt("BitNumber", 3);

	fMinJetPt = GetDouble("MinJetPt", 500.);
	fMinMuonPt = GetDouble("MinMuonPt", 10.);

	fMinTowerPtRatio = GetDouble("MinTowerPtRatio", .05);
	fCoreAntiktR = GetDouble("CoreAntiktR", .04);
	fCorePtRatioMin = GetDouble("CorePtRatioMin", .5);

	fMinCoreMinBoost = GetDouble("MinCoreMinBoost", 1.);
	fMinCoreMinBoost2 = Squared(fMinCoreMinBoost);

	fCoreMassHypothesis = GetDouble("CoreMassHypothesis", 2.0);
	fCoreMassHypothesis2 = Squared(fCoreMassHypothesis);

	fMinFinalMass = GetDouble("MinFinalMass", 5.3);
	fMaxFinalMass = GetDouble("MaxFinalMass", fMinFinalMass*2);

	fMaxEmissionInvariant = GetDouble("MaxX", 3.);
	cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Emission Invariant: " << fMaxEmissionInvariant << endl;

	fCoreDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm, fCoreAntiktR);

	// import input array(s)
	fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
	fItJetInputArray = fJetInputArray->MakeIterator();

	fAllParticles = ImportArray("Delphes/allParticles");
}

//------------------------------------------------------------------------------

void HighPtBTagger::Finish()
{
	delete fCoreDefinition;
	delete fItJetInputArray;
}

//------------------------------------------------------------------------------

// After processing:
//     Jet with:
//        good muons (potential for tagging) ........... TauTag == numMuons
//        viable HardCore .............................. FracPt[0] > 0.
//        viable MinCore (which isn't the HardCore) .... FracPt[1] > 0.
//     Muon:
//        good muon in jet with interesting pt ......... BTag/TauTag non-zero
//        good muon in jet with interesting core ....... Tau[2] > 0. / Tau[0] > 0. / FracPt[0] < 8.
//        good muon in jet with distinct MinCore ....... Tau[1] > 0. / FracPt[1] < 8.

void HighPtBTagger::Process()
{
	// loop over input jets
	Candidate* jet;
	fItJetInputArray->Reset();

	while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
	{
		const Double_t jetOriginalPt = jet->Momentum.Pt();

		// Set used fields to NULL values
		jet->FracPt[0] = -1.;
		jet->FracPt[1] = -1.;
		jet->Tau[0] = -1;
		jet->Tau[1] = -1;
		jet->TauTag = 0;

		// At best, neutrino estimation can only double the pt of the jet,
		// so only look at jets with at least half the required final pt
		if(jetOriginalPt >= .5*fMinJetPt)
		{
			// We'll sort jet constituents into two categories (where goodMuons are those which pass the pt cut)
			std::vector<Candidate*> goodMuons;
			std::vector<Candidate const*> everythingElse;

			// Fill goodMuons & everythingElse (CHECKED 14.04.2015)
			{
				TIter itJetConstituents(jet->GetCandidates());
				Candidate* constituent;

				everythingElse.reserve(16);

				// Find all muons that pass the pt cuts
				while((constituent = static_cast<Candidate*>(itJetConstituents.Next())))
				{
					if(abs(constituent->PID) == 13)
					{
						if((constituent->Momentum).Pt() >= fMinMuonPt)
						{
							// Set used fields to NULL values;
							constituent->Tau[0] = -1.;
							constituent->Tau[1] = -1.;
							constituent->Tau[2] = -1.;
							constituent->FracPt[0] = 8.;
							constituent->FracPt[1] = 8.;
							goodMuons.push_back(constituent);
							continue;
						}
					}
					everythingElse.push_back(constituent);
				}
			}

			if(not goodMuons.empty()) // Ensure we have at least one good muon
			{
				jet->TauTag = goodMuons.size();

				// Sort muons low to high (CHECKED 14.04.2015)
				{
					// Sometimes there will be more than 1 muon. In those cases, we will
					// assume that the muons were emitted from highest to lowest pt.
					// This means we'll need to add them back from low to high pt, so that
					// the CM frame of the muon emitted later doesn't include the previously emitted muon.

					// Sort muons low to high
					std::sort(goodMuons.begin(), goodMuons.end(), SortCandidatePt_Low2High);
				}

				std::vector<Candidate const*> goodMuonsMatriarch;
				// Find the muon matriarch and mother (CHECKED 14.04.2015)
				// Store inside each muon
				//       Muon.BTag = abs(motherPID)
				//       Muon.TauTag = abs(matriarchPID)
				//       Muon.Tau[2] = xTrue (real x to matriarch)
				{
					// * The "matriarch" is the original boosted primary hadron
					//      We need to look back all the way to the start of hadronization.
					// * The "mother" is the particle which emitted the muon
					//      If the mother is a tau, we consider the grandmother the mother (because we are interested in hadron flavor)
					for(std::vector<Candidate*>::iterator itMuon = goodMuons.begin(); itMuon != goodMuons.end(); ++itMuon)
					{
						bool noMotherYet = true;

						Candidate const* matriarch;
						{
							// To look back all the way to hadronization, we find the first matriarch with
							// two seperate mothers (more than 1 mother means color charge still existed)
							Candidate const* possibleMatriarch = static_cast<Candidate*>(fAllParticles->At((*itMuon)->M1));

							do
							{
								const UInt_t absPossibleMatriarchPID = abs(possibleMatriarch->PID);
								matriarch = possibleMatriarch;

								// Tau cannot be mother, we want hadrons
								if(noMotherYet and (absPossibleMatriarchPID not_eq 15))
								{
									noMotherYet = false;
									(*itMuon)->BTag = absPossibleMatriarchPID;
								}

								if( (not ((possibleMatriarch->M2 == 0) or (possibleMatriarch->M2 == possibleMatriarch->M1))) or
									((absPossibleMatriarchPID >= 22) and (absPossibleMatriarchPID <= 24)))
								{
									// First check (post hadronic)
									// Post hadronic decays should have a M2 == 0 or (in the case of interaction
									// with the vacuum) M1 == M2. Anthing else and we've reached back into color charge

									// Second check (primary muons from A, Z, W)
									// Pythia hadronic decay never specifically invokes a W or Z
									// so any encountered here must have come from the hard interaction.
									// Pythia DOES produce A which can split to muons.

									// However, if we break here, the mother should already be set
									break; // out of the do loop
								}
								else // Keep tracing back the lineage
								{
									possibleMatriarch = static_cast<Candidate*>(fAllParticles->At(possibleMatriarch->M1));
								}
							}while(possibleMatriarch);

							// Record the muon's matriarch
							goodMuonsMatriarch.push_back(matriarch);
							(*itMuon)->TauTag = abs(matriarch->PID);

							if(noMotherYet)
							{
								// If I understand Pythia correctly, we should never reach here
								// So let me know if we ever reach here.
								cout << "\n\nMotherFail\n\n" << endl;
								(*itMuon)->BTag = (*itMuon)->TauTag;
							}

							// Calculate the muon's true x to the matriarch
							// Store in Tau[2] (repurposed N-sub-jettyness)
							const TVector3 matriarchP3 = (matriarch->Momentum).Vect();
							const TVector3 muonP3 = ((*itMuon)->Momentum).Vect();

							(*itMuon)->Tau[2] = (matriarch->Momentum).E() *
								(muonP3.Cross(matriarchP3)).Mag() /
								((matriarch->Momentum).M() * muonP3.Dot(matriarchP3));
						}
					}
				}

				// Now we find the jet's core by reculstering its jet constituents
				vector<fastjet::PseudoJet> reclusterInput;

				// Add the good muons to the reclusterInput (CHECKED 14.04.2015)
				{
					// Give them a negative user_index ( user_index = index in goodMuons - goodMuons.size() )
					// This makes them easy to find after the reclustering.

					const int numMuons = goodMuons.size();
					for(int iMu = 0; iMu < numMuons; ++iMu)
					{
						const TLorentzVector& muonMomentum = goodMuons[iMu]->Momentum;
						reclusterInput.emplace_back(muonMomentum.Px(), muonMomentum.Py(), muonMomentum.Pz(), muonMomentum.E());
						reclusterInput.back().set_user_index(iMu - numMuons);
					}
				}

				// Add jet constituents (all tracks towers passing pt cut) (CHECKED 14.04.2015)
				{
					// Use all tracks, but only use towers/eFlowNeutrals (charge == 0) if they are above the relative pt threshold
					// This is because the angular resolution of towers is much poorer, and we don't
					// want soft radiation unduly influencing the direction of the core

					const Double_t minTowerPt = fMinTowerPtRatio*jetOriginalPt;
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

				// We'll be working with 2 core definitions. Keeping them in a vector allows us to re-use the code more easily
				//    index 0: The Hardcore (the subjet with the highest pt)
				//    index 1: The MinCore (the core which minimizes the reconstructed mass, when adding the muon, if different than the HardCore)
				std::vector<TLorentzVector> coreP4Vec;
				// We'll also keep track of the minimum x we get for each core definition, for the final tag
				std::vector<Double_t> coreMinX;

				{// Scope of fastjet ClusterSequence
					// Only consider subjets if they have a boost greater than fMinCoreMinBoost
					// pt not_eq E, let's turn this off for now
					// const Double_t minSubjetPt = fMinCoreMinBoost*fCoreMassHypothesis;
					const Double_t minSubjetPt = 0.;

					// Recluster the jet, to find subjets, then ensure the pt_subjet > 10 * coreMassHypothesis
					fastjet::ClusterSequence recluster(reclusterInput, *fCoreDefinition);
					const std::vector<fastjet::PseudoJet> subJets = fastjet::sorted_by_pt(recluster.inclusive_jets(minSubjetPt));

					// Make sure at least one subjet passed the cut
					if(not subJets.empty())
					{
						// Keep pseudojets in a vector that parallels coreP4Vec
						std::vector<fastjet::PseudoJet> coreVec;

						// The hardCore is the core with the highest pt
						coreVec.push_back(subJets.front());
						coreMinX.push_back(9e9);

						// Get fastjets internal jets (because it can't alter our input jets,
						// so they have no internal information, which we'll need to use).
						// The documentation does not specify that the input jets are the first
						// members of list returned by jets(), so we'll have to assume that
						// this is the case (though we will check this assumption later)
						const std::vector<fastjet::PseudoJet>& internalJets = recluster.jets();

						// Find the minCore (WARNING, only the hardest muon is used to define the MinCore) (CHECKED 14.04.2015)
						{
							// Unfortunately, we can't just add the muon to the core and calculate
							// the mass, because the core currently has the wrong mass (from the granularity of
							// the cal, not from the actual particles). And we don't really want to recompute
							// the 4-momentum of every subjet with a new mass.

							// However, from a rather trivial calculation, one can show that, if the mass of
							// the subjet is constrained to CoreMassHypothesis, the M**2 after adding the muon twice is:
							//
							// M**2 = CoreMassHypothesis**2 + 4*Emu*Esub*(g + y)/(1 + (y + sqrt(1 - ((g-y) + g*y))))
							//
							// where (g = (CoreMassHypothesis/Esub)**2) and (y = tan(theta)**2 = (p3>_mu x p3>_sub)**2 / (p3>_mu.p3>_sub)**2
							//
							// Luckily, we can calculate y without applying the CoreMassHypothesis to each subjet,
							// since their direction won't be altered by a new mass.

							const fastjet::PseudoJet& hardestMuon = internalJets[goodMuons.size() - 1];
							// With the indexing scheme used when filling reclusterInput, the
							// hardest muon should have a user_index of -1. Verify.
							if(hardestMuon.user_index() not_eq -1)
								throw runtime_error("HighPtBTagger: Can't access original muon from vector returned by ClusterSequence::jets()!");

							// The hardest muon can only be inside one sub-jet, so we can
							// stop looking to see if it's inside subjets after it's
							// already found and subtracted from one of them
							bool hardestMuonUnsubtracted = true;

							// By minimizing the squared mass, we'll minimize the mass
							Double_t minMass2 = 9e9;
							unsigned int minCoreIndex = 0;

							for(unsigned int iSub = 0; iSub < subJets.size(); ++iSub)
							{
								// Make a copy of the subjet, so we can operate on its 4-vector
								fastjet::PseudoJet subjet = subJets[iSub];

								// If the muon is already inside the subjet, subtract it's momentum
								if(hardestMuonUnsubtracted and hardestMuon.is_inside(subjet))
								{
									subjet -= hardestMuon;
									hardestMuonUnsubtracted = false;
								}

								const Double_t g = Squared(fCoreMassHypothesis / subjet.E());
								// Make sure we still have sufficient boost after muon subtraction;
								if(g*fMinCoreMinBoost2 > 1.)
									continue;

								const Double_t y = Tan2(subjet, hardestMuon);

								// Technically, every squared mass has CoreMassHypothesis**2 (and is multiplied by 4),
								// so there is no point in doing these to everything
								// const Double_t mass2 = fCoreMassHypothesis2 + 4. * hardestMuon.E() * subjet.E() * (g + y) / (1. + (y + sqrt(1. - ((g - y) + g*y))));

								//WARNING: Currently this implementation does not account for cores with two muons of almost equal energy (small percentage of cores)
								const Double_t mass2 = hardestMuon.E() * subjet.E() * (g + y) / (1. + (y + sqrt(1. - ((g - y) + g*y))));

								if(mass2 < minMass2)
								{
									minMass2 = mass2;
									minCoreIndex = iSub;
								}
							}

							// Check if the minCore is different than the hardCore
							if(minCoreIndex > 0)
							{
								coreVec.push_back(subJets[minCoreIndex]);
								coreMinX.push_back(9e9);
							}
						}
						// If the minCore is the hardcore, coreVec will only have 1 member

						// Subtract muons from the cores (CHECKED 14.04.2015)
						{
							for(unsigned int iMu = 0; iMu < goodMuons.size(); ++iMu)
							{
								const fastjet::PseudoJet& muon = internalJets[iMu];
								if(muon.user_index() >= 0)
									throw runtime_error("HighPtBTagger: Muon ID Fail!");

								for(std::vector<fastjet::PseudoJet>::iterator itCore = coreVec.begin(); itCore not_eq coreVec.end(); ++itCore)
								{
									if(muon.is_inside(*itCore))
									{
										*itCore -= muon;
										break; // out of inner loop - the muon can only be in one core at a time
									}
								}
							}
						}

						// From now on, we will keep the original PseudoJet cores around (i.e. coreVec)
						// in case we are interested in constituents. However, all 4-vector math
						// will now be handled by TLorentzVector, since it has more/better math.

						// Fix core mass to fCoreMassHypothesis (CHECKED 14.04.2015)
						{
							for(std::vector<fastjet::PseudoJet>::iterator itCore = coreVec.begin(); itCore not_eq coreVec.end(); ++itCore)
							{
								// The current mass depends more on the granularity of the CAL than anything else
								// To fix the mass, We need to scale the momentum of the core (but not the energy)
								const Double_t momentumScale = sqrt((Squared(itCore->E()) - fCoreMassHypothesis2)/itCore->modp2());
								coreP4Vec.emplace_back(momentumScale*itCore->px(), momentumScale*itCore->py(), momentumScale*itCore->pz(), itCore->E()); //xyzt
							}
						}

						//const bool minCoreIsHardCore = (coreVec.size() == 1);

						// Iterate through the cores and analyze
						for(unsigned int iCore = 0; iCore < coreP4Vec.size(); ++iCore)
						{
							TLorentzVector& core = coreP4Vec[iCore];

							// Iterate through each muon, adding it back to the core, and calculate the resulting boost invariant
							for(unsigned int iMu = 0; iMu < goodMuons.size(); ++iMu)
							{
								Candidate* muon = goodMuons[iMu];
								const TLorentzVector& muonP4 = muon->Momentum;
								const TVector3 muonP3 = muonP4.Vect();
								//const Double_t originalMass = core.M();

								// Add the muon back to the core
								core += muonP4;

								// Now estimate the neutrino
								{
									// Only the HardCore adds the "neutrino" to the original jet
									if(iCore == 0)
										jet->Momentum += muonP4;

									// Figure out if adding the muon twice helps orient the subjet
									{
										const TVector3 matriarchP3 = (goodMuonsMatriarch[iMu]->Momentum).Vect();

										const Double_t originalAngle = AccurateAngle(matriarchP3, core.Vect());
										// Add the "neutrino" to the core
										core += muonP4;
										// Store the change in angle
										muon->FracPt[iCore] = AccurateAngle(matriarchP3, core.Vect()) - originalAngle;
									}
								}

								// Find the deltaMass
								// const Double_t delta_Mass = core.M() - originalMass;

								// Now calcuale x
								const Double_t boostMass = std::min(std::max(core.M(), fMinFinalMass), fMaxFinalMass);
								const Double_t	xCore = (core.E() * (muonP3.Cross(core.Vect())).Mag()) / (muonP3.Dot(core.Vect()) * boostMass);

								muon->Tau[iCore] = xCore;
								coreMinX[iCore] = std::min(coreMinX[iCore], xCore);
							}// End loop over muons
						}// End loop over cores
					}// End valid subjets
				}// End Reclustering

				const Double_t	jetPt = (jet->Momentum).Pt();

				// Re-check pt now that we've estimated neutrinos
				if(jetPt > fMinJetPt)
				{
					for(unsigned int iCore = 0; iCore < coreP4Vec.size(); ++iCore)
					{
						jet->FracPt[iCore] = coreP4Vec[iCore].Pt() / jetPt; // Store core ratio
						jet->Tau[iCore] = coreMinX[iCore]; // Store core's smallest x
					}//End loop through cores

					// If the MinCore is the HardCore, copy its values
					if (coreP4Vec.size() == 1)
					{
						jet->FracPt[1] = jet->FracPt[0];
						jet->Tau[1] = jet->Tau[0];
					}

					// Tag based on the MinCore (make a flag maybe)
					const unsigned int taggingCore = 1;

					if(jet->FracPt[taggingCore] >= fCorePtRatioMin)
					{
						if(jet->Tau[taggingCore] <= fMaxEmissionInvariant)
						{
							jet->BTag |= (1 << fBitNumber);
						}//End tagged
					}//End core is hard enough
				}
			}//End muon pt cut
		}//End jet initial pt cut
	}//End loop over jets
}

//------------------------------------------------------------------------------

bool SortCandidatePt_Low2High(Candidate const* one, Candidate const* two)
{
	return (one->Momentum).Pt() < (two->Momentum).Pt();
}

/* Old Neutrino Estimation Code
{
	TVector3 neutrino;
	{
		// Since the angle between the core and the muon is extremely important,
		// we don't want to bias the core towards the muon. Thus, instead of
		// adding the muon twice, we give the neutrino the direction of the core,
		// projecting out the magnitude of the muon momentum.
		//     neutrinoP3> = (muonP3>.coreP3> / coreP3>.coreP3>) * coreP3>
		const TVector3 coreP3 = coreRaw.Vect(); // raw and fixed should point in same direction
		neutrino = coreP3;
		neutrino *= coreP3.Dot(muonP3)/coreP3.Mag2();
	}

	{
		TLorentzVector neutrinoMomentum(neutrino.Px(), neutrino.Py(), neutrino.Pz(), neutrino.Mag()); //xyzt

		// Add the nuetrino to the core AND to the jet
		coreRaw += neutrinoMomentum;
		coreFixed += neutrinoMomentum;
		jet->Momentum += neutrinoMomentum;
	}
}
*/
