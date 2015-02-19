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

//------------------------------------------------------------------------------

HighPtBTagger::HighPtBTagger():
	fMinJetPt(0.), fMinMuonPt(0.), fCoreAntiktR(0.), fMinCoreRatio(0.), fMinCoreMass(0.), fMaxEmissionInvariant(0.),
	fJetInputArray(0), fItJetInputArray(0) {}

//------------------------------------------------------------------------------

HighPtBTagger::~HighPtBTagger()
{ /* No memory to clean-up */}

//------------------------------------------------------------------------------

void HighPtBTagger::Init()
{
 	// read parameters

	fMinJetPt = GetDouble("MinJetPt", 1000.);
	fMinMuonPt = GetDouble("MinMuonPt", 15.);
	fCoreAntiktR = GetDouble("CoreAntiktR", .2);
	fMinCoreRatio = GetDouble("MinCoreRatio", .8);
	fMinCoreMass = GetDouble("MinCoreMass", 5.3);
	
	{
		const Double_t sigma = abs(GetDouble("Sigma", .9));
		
		if(sigma >= 1.)
		{	
			throw runtime_error("HightPtBTagger. Supplied sigma is out of bounds (must be in [0,1) )");
		}
		
		fMaxEmissionInvariant = 1./tan(sqrt(1. - sigma));
		
		cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Emission Invariant: " << fMaxEmissionInvariant << endl;
	}
	
	fJetDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm, fCoreAntiktR);
	
	// import input array(s)
	fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
	fItJetInputArray = fJetInputArray->MakeIterator();
	
	jetsAboveThreshold = 0;
	jetsWithGoodMuons = 0;
	jetsTagged = 0.;
	invariants.clear();
}

//------------------------------------------------------------------------------

void HighPtBTagger::Finish()
{
	if(not invariants.empty())
	{
		cout << std::fixed;
		cout << "\n\n\n{ " << invariants[0];
	
		for(unsigned int i = 1; i < invariants.size(); ++i)
		{
			cout << ", " << invariants[i];
			if(i%10 == 9)
				cout << "\n";
		}
		cout << "}\n\n\n";
	}
	
	delete fJetDefinition;
	delete fItJetInputArray;
}

//------------------------------------------------------------------------------

void HighPtBTagger::Process()
{
	vector<fastjet::PseudoJet> inputList;
	
	// loop over input jets
	Candidate* jet;
	fItJetInputArray->Reset();
	
	stringstream debugOut;
	
	const int granularity = 4; // granularity of ATLAS middle sample, 4x4 ECAL cells per HCAL cell
	const Double_t eCalThreshold = 20.;
	const Double_t eCalRadiusFactor = .5;
	const Double_t invariantScaleFactor = 2./3.;
	
	const fastjet::JetDefinition eCalJetDefinition(fastjet::antikt_algorithm, fCoreAntiktR * eCalRadiusFactor);
		
	while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
	{
		// Reset the output stream 
		debugOut.str("");
		debugOut << setprecision(4);
	
		const Double_t jetPt = (jet->Momentum).Pt();
		
		if(jetPt >= fMinJetPt) // Ensure the jet is above the pt cut
		{
			++jetsAboveThreshold;
			debugOut << "pt/eta/phi:         " << jetPt << "\t| " << (jet->Momentum).Eta() << "\t| " << (jet->Momentum).Phi() << endl;
		
			TObjArray const* const jetConstituents = jet->GetCandidates();
			TIterator* const itJetConstituents = jetConstituents->MakeIterator();
			std::vector<Candidate const*> goodMuons;
			
			// Find all muons that pass the pt cuts
			{
				Candidate const* constituent;
				
				while((constituent = static_cast<Candidate const*>(itJetConstituents->Next())))
				{
					if((abs(constituent->PID) == 13) and ((constituent->Momentum).Pt() >= fMinMuonPt))
						goodMuons.push_back(constituent);
				}
			}
			
			if(not goodMuons.empty()) // Ensure we have at least one muon with enough pt
			{
				++jetsWithGoodMuons;
				// To find the hard core, we could cluster HCAL cells. But the problem
				// is that we're looking for a hard core whose radius is likely smaller
				// than the size of HCAL cells. To improve the angular resolution, we 
				// can make a "reasonable" assumption that a highly boosted jet will have
				// at least one pi0, and that this energy will get boosted into significant
				// energy deposits in the ECAL. Thus, ECAL hits can be used for pointing 
				// information for all the energy of the tower.
				// Initialy, this algorithm was implemented on a cell-by-cell basis, 
				// using each tower's ECAL hits to direct it's HCAL energy.
				// In the second iteration, we will allow the ECAL energy from 
				// neighboring towers to cluster. This ECAL energy will then pull HCAL
				// energy along with it (in an amount proportional to its ECAL propotion 
				// in the tower).
				// To reduce noise, we will require the energy in an ECAL cell to 
				// exceed a given threshold before it can participate in the clustering.
								
				std::vector<ECalTower> eCalTowers;
				std::vector<Candidate const*> standaloneHCalTowers;
				{
					Candidate* tower;			
					
					itJetConstituents->Reset();
					while((tower = static_cast<Candidate*>(itJetConstituents->Next())))
					{
						// Make sure it has meaningful energy and is not a track
						if(((tower->Momentum).Pt() >= .01*jetPt) and (tower->Charge == 0))
						{
							// Find the edges of the tower (HCAL cell)
							const Double_t
								etaMin = tower->Edges[0],
								etaMax = tower->Edges[1],
								phiMin = tower->Edges[2],
								phiMax = tower->Edges[3];
						
							// Find the width of each ECAL cell
							const Double_t
								etaWidth = (etaMax - etaMin)/granularity,
								phiWidth = (phiMax - phiMin)/granularity;
								 
							// Create the ECAL grid
							// Linear grid would be more efficient, implement later
							Double_t eCalGrid[granularity][granularity];
							
							// Zero the grid
							for(int eta = 0; eta < granularity; ++eta)
								for(int phi = 0; phi < granularity; ++phi)
									eCalGrid[eta][phi] = 0.;
							
							// Prepare to loop through the tower's constituents
							TIterator* fItTowerConstituents = (tower->GetCandidates())->MakeIterator();
							Candidate const* towerConstituent;
													
							// Loop over all particles and bin them into ECAL grids
							while((towerConstituent = static_cast<Candidate*>(fItTowerConstituents->Next())))
							{
								const int absID = abs(towerConstituent->PID);
								
								if((absID == 11) or (absID == 22)) // electrons and photons
								{
									// Cal eta/phi depends on position, not momentum
									const TLorentzVector& thisPosition = towerConstituent->Position;
									const int 
										thisEtaIndex = int((thisPosition.Eta() - etaMin)/etaWidth),
										thisPhiIndex = int((thisPosition.Phi() - phiMin)/phiWidth);
								
									eCalGrid[thisEtaIndex][thisPhiIndex] += (towerConstituent->Momentum).E();
									
									if(absID == 11)
										debugOut << "Electron:           " << (towerConstituent->Momentum).E() << "\n";
								}
							}
						
							TLorentzVector eCalMomentum;
							Double_t qualifyingEnergy = 0.;
							const unsigned int firstETowerIndex = eCalTowers.size();
						
							// Now loop over all ECAL grids and, if they are above the 
							// energy threshold, creat an ECalTower
							for(int eta = 0; eta < granularity; ++eta)
							{
								for(int phi = 0; phi < granularity; ++phi)
								{
									const Double_t eCalEnergy = eCalGrid[eta][phi];
									
									if(eCalEnergy >= eCalThreshold)
									{
										const Double_t newEta = etaMin + (0.5 + eta)*etaWidth;
										eCalMomentum.SetPtEtaPhiM(eCalEnergy/cosh(newEta), newEta, phiMin + (0.5 + phi)*phiWidth, 0.);
										eCalTowers.emplace_back(eCalMomentum, eCalEnergy*(tower->Ehad));
										qualifyingEnergy += eCalEnergy;
									}
								}
							}
							
							const unsigned int pastLastETowerIndex = eCalTowers.size();
							
							if(pastLastETowerIndex == firstETowerIndex)
							{
								standaloneHCalTowers.push_back(tower);
								debugOut << "Standalone:         " << tower->Momentum.E() << "\n";
							}
							else
							{
								// Now loop over any newly created ECalTowers and scale their hCalEnergy by 
								// the amount of ECalEnergy that qualified
								for(unsigned int eTowerIndex = firstETowerIndex; eTowerIndex < pastLastETowerIndex; ++eTowerIndex)
									eCalTowers[eTowerIndex].hCalEnergy /= qualifyingEnergy;
								debugOut << "EM pointed:         " << "{ " << tower->Eem << " , " << tower->Ehad << " }\n";
							}
						
							delete fItTowerConstituents;
						}
					}
				}
				
				{
					inputList.clear();
					int numConstituents = 0;
	
					// Now feed all ECalTowers into fastjet
					for(std::vector<ECalTower>::const_iterator itETower = eCalTowers.begin(); itETower not_eq eCalTowers.end(); ++itETower)
					{
						const TLorentzVector& eTowerP4 = itETower->eCalMomentum;
						
						inputList.emplace_back(eTowerP4.Px(), eTowerP4.Py(), eTowerP4.Pz(), eTowerP4.E());						
						inputList.back().set_user_index( numConstituents++ );
					}
				}
				
				{	
					// Now cluster the ECalTowers
					fastjet::ClusterSequence clusterSequence(inputList, eCalJetDefinition);
					std::vector<fastjet::PseudoJet> outputList = fastjet::sorted_by_pt(clusterSequence.inclusive_jets(0.));
					
					if(not outputList.empty())
						debugOut << "highest energy EM: " << outputList.front().E() << "\n";
					
					// Now loop through each output jet
					for(std::vector<fastjet::PseudoJet>::iterator itOutJet = outputList.begin(); itOutJet not_eq outputList.end(); ++itOutJet)
					{
						// Find its constituents
						std::vector<fastjet::PseudoJet> inJet = clusterSequence.constituents(*itOutJet);
						
						Double_t hCalAddition = 0.;
						// Now loop through those constituents and add up all the HCalEnergy we need to add
						for(std::vector<fastjet::PseudoJet>::const_iterator itInJet = inJet.begin(); itInJet not_eq inJet.end(); ++itInJet)
							hCalAddition += eCalTowers[itInJet->user_index()].hCalEnergy;
						
						// Now scale the 4-vectors to bring their HCAL energy along with them
						// We'll do this by discarding any aqcuired mass and use the momentum
						// as a massless 4-vector
						
						const Double_t p3Mag = sqrt(itOutJet->modp2());
						const Double_t newScale = (1. + hCalAddition/p3Mag);
						itOutJet->reset_momentum(itOutJet->px()*newScale, itOutJet->py()*newScale, itOutJet->pz()*newScale, p3Mag*newScale);						
					}
					
					// Now clear the inputlist and re-fill it
					inputList.clear();
					int numConstituents = 0;
					// Start with the muons, so that their user_index is easy to find (less than goodMuons.size())
					for(std::vector<Candidate const*>::const_iterator itMuons = goodMuons.begin(); itMuons not_eq goodMuons.end(); ++itMuons)
					{
						const TLorentzVector& thisMuon = (*itMuons)->Momentum;
						inputList.emplace_back(thisMuon.Px(), thisMuon.Py(), thisMuon.Pz(), thisMuon.E());
						inputList.back().set_user_index( numConstituents++ );
					}
					debugOut << "good muons:         " << goodMuons.size() << "\n";
					
					// Now add the ECAL directed HCAL towers
					for(std::vector<fastjet::PseudoJet>::const_iterator itOutJet = outputList.begin(); itOutJet not_eq outputList.end(); ++itOutJet)
					{
						inputList.emplace_back(itOutJet->px(), itOutJet->py(), itOutJet->pz(), itOutJet->E());
						inputList.back().set_user_index( numConstituents++ );
					}
					debugOut << "ECAL jets:          " << outputList.size() << "\n";
					
					// Don't add standalone HCAL towers (using the center of the tower)
					/*
					for(std::vector<Candidate const*>::const_iterator itStandalone = standaloneHCalTowers.begin(); itStandalone not_eq standaloneHCalTowers.end(); ++itStandalone)
					{
						const TLorentzVector& thisTower = (*itStandalone)->Momentum;
						inputList.emplace_back(thisTower.Px(), thisTower.Py(), thisTower.Pz(), thisTower.E());
						inputList.back().set_user_index( numConstituents++ );
					}
					*/
					debugOut << "standalone HCAL:    " << standaloneHCalTowers.size() << "\n";
				}
				
				// Now cluster the new list
				fastjet::ClusterSequence clusterSequence(inputList, *fJetDefinition);
				std::vector<fastjet::PseudoJet> outputList = fastjet::sorted_by_pt(clusterSequence.inclusive_jets( jetPt * fMinCoreRatio ));
				
				if(not outputList.empty()) // Make sure at least one subjet passed the pt cut
				{
					fastjet::PseudoJet coreCopy = outputList.front(); // The highest pt jet is the core
										
					debugOut << "core Pt:            " << sqrt(coreCopy.kt2()) << endl;
					debugOut << "core mass:          " << coreCopy.m() << endl;
					
					{
						// Find the number of core constituents
						std::vector<fastjet::PseudoJet> coreConstituents = clusterSequence.constituents(outputList.front());
						
						debugOut << "coreConstituents:   " << coreConstituents.size() << endl;
						
						// Find any muons in the core and subtract their 4-momentum
						for(std::vector<fastjet::PseudoJet>::iterator itCore = coreConstituents.begin(); itCore != coreConstituents.end(); ++itCore)
						{
							if(itCore->user_index() < goodMuons.size())
								coreCopy -= *itCore;
						}
					}
					
					TLorentzVector coreMomentum(coreCopy.px(), coreCopy.py(), coreCopy.pz(), coreCopy.E()); //xyzt
						
					// Sometimes there will be more than 1 muon. In those cases, we will
					// assume that the muons were emitted from highest to lowest pt.
											
					// Sort muons low to high
					std::sort(goodMuons.begin(), goodMuons.end(), SortCandidatePt_Low2High);
					
					debugOut << "\n";
					
					// Now add back each muon's 4-momentum, estimate its neutrino, and calculate the resulting boost invariant
					for(std::vector<Candidate const*>::iterator itMuon = goodMuons.begin(); itMuon != goodMuons.end(); ++itMuon)
					{
						TVector3 muonP3;
						{
							const TLorentzVector& muonMomentum = (*itMuon)->Momentum;
							muonP3 = muonMomentum.Vect();
						
							// Add the muon back to the core
							coreMomentum += muonMomentum;
							
							debugOut << "muon Pt:            " << muonMomentum.Pt() << endl;
						}
						
						TVector3 neutrino;
						{
							// We do not want to bias the core momentum towards the muon emission direction
							// Thus, we give the neutrino the direction of the core, projecting out the 
							// magnitude of the muon momentum.
							//     neutrinoP3> = (muonP3>.coreP3> / coreP3>.coreP3>) * coreP3>
							const TVector3 coreP3 = coreMomentum.Vect();
							neutrino = coreP3;							
							neutrino *= coreP3.Dot(muonP3)/coreP3.Mag2();
						}
						
						{
							TLorentzVector neutrinoMomentum(neutrino.Px(), neutrino.Py(), neutrino.Pz(), neutrino.Mag());
							
							// Don't add the neutrino to the core, since it won't change its direction							
							// coreMomentum += neutrinoMomentum;
							// Also add the neutrino the jet;
							jet->Momentum += neutrinoMomentum;
							
							debugOut << "neutrino Pt:        " << neutrinoMomentum.Pt() << endl;
						}
						
						// Find the boost of the core, disallowing any mass lower than (fMinCoreMass)
						const Double_t coreBoostMass = std::max(coreMomentum.M(), fMinCoreMass);
						const Double_t coreBoost = coreMomentum.E()/coreBoostMass;						
													
						debugOut << "core Pt:            " << coreMomentum.Pt() << endl;
						debugOut << "core Mass:          " << coreBoostMass << endl;							
						
						const TVector3 coreP3 = coreMomentum.Vect();
						const Double_t CMInvariant = invariantScaleFactor * coreBoost * (muonP3.Cross(coreP3)).Mag() / muonP3.Dot(coreP3);
						
						if(CMInvariant <= fMaxEmissionInvariant)
						{
							++jetsTagged;
						}
						
						debugOut << "CM invaraint:       " << CMInvariant << endl;
						debugOut << "\n";
						
						invariants.push_back(CMInvariant);
						
					}// End loop over muons
					debugOut << "w/goodMuon efficiency...." << jetsTagged/jetsWithGoodMuons << "\n";
					debugOut << "Total efficiency........." << jetsTagged/jetsAboveThreshold << "\n\n"; 
					
					cout << debugOut.str();
					cout << "\n---------------------------------------------------------\n";				
				}// End core pt cut					
			}//End Muon pt cut			
			
			delete itJetConstituents;			
		}//End jet pt cut  	
	}//End loop over jets
}

//------------------------------------------------------------------------------

ECalTower::ECalTower(const TLorentzVector& eCalMomentum_in, const Double_t hCalEnergy_in):
	eCalMomentum(eCalMomentum_in), hCalEnergy(hCalEnergy_in) {}

//------------------------------------------------------------------------------

bool SortCandidatePt_Low2High(Candidate const* one, Candidate const* two)
{
	return (one->Momentum).Pt() < (two->Momentum).Pt();
}

//------------------------------------------------------------------------------

/* Old algorithm, use EM cells to point each hadronic cell individually.


void HighPtBTagger::Process()
{
	vector<fastjet::PseudoJet> inputList;
	
	// loop over input jets
	Candidate* jet;
	fItJetInputArray->Reset();
	
	stringstream debugOut;
		
	while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
	{
		// Reset the output stream 
		debugOut.str("");
		debugOut << setprecision(4);
	
		const Double_t jetPt = (jet->Momentum).Pt();
		
		if(jetPt >= fMinJetPt) // Ensure the jet is above the pt cut
		{
			debugOut << "pt/eta/phi:         " << jetPt << "\t| " << (jet->Momentum).Eta() << "\t| " << (jet->Momentum).Phi() << endl;
		
			TObjArray const* const jetConstituents = jet->GetCandidates();
			TIterator* const itJetConstituents = jetConstituents->MakeIterator();
			std::vector<Candidate const*> goodMuons;
			
			// Find all muons that pass the pt cuts
			{
				Candidate const* constituent;
				
				while((constituent = static_cast<Candidate const*>(itJetConstituents->Next())))
				{
					if((abs(constituent->PID) == 13) and ((constituent->Momentum).Pt() >= fMinMuonPt))
						goodMuons.push_back(constituent);
				}
			}
			
			if(not goodMuons.empty()) // Ensure we have at least one muon with enough pt
			{
				// Now loop over all jet constituents again. If they're pt > 1% of jet pt,
				// bin their ECal hits and use the enhanced angular information to reset
				// their direction and mass.
				itJetConstituents->Reset();
				
				{
					Candidate* tower;
					const int granularity = 4; // granularity of ATLAS middle sample
				
					while((tower = static_cast<Candidate*>(itJetConstituents->Next())))
					{
						// Make sure it's not a track
						if(((tower->Momentum).Pt() >= .01*jetPt) and (tower->Charge == 0))
						{
							const Double_t
								etaMin = tower->Edges[0],
								etaMax = tower->Edges[1],
								phiMin = tower->Edges[2],
								phiMax = tower->Edges[3];
						
							const Double_t
								etaWidth = (etaMax - etaMin)/granularity,
								phiWidth = (phiMax - phiMin)/granularity;
								 
							Double_t eCalGrid[granularity][granularity];
					
							for(int eta = 0; eta < granularity; ++eta)
								for(int phi = 0; phi < granularity; ++phi)
									eCalGrid[eta][phi] = 0.;
					
							TIterator* fItTowerConstituents = (tower->GetCandidates())->MakeIterator();
							Candidate const* towerConstituent;
						
							// Loop over all particles and bin them into ECAL grids
							while((towerConstituent = static_cast<Candidate*>(fItTowerConstituents->Next())))
							{
								const int absID = abs(towerConstituent->PID);
								
								if((absID == 11) or (absID == 22)) // electrons and photons
								{
									// Cal eta/phi depends on position, not momentum
									const TLorentzVector& thisPosition = towerConstituent->Position;
									const int 
										thisEtaIndex = int((thisPosition.Eta() - etaMin)/etaWidth),
										thisPhiIndex = int((thisPosition.Phi() - phiMin)/phiWidth);
								
									eCalGrid[thisEtaIndex][thisPhiIndex] += (towerConstituent->Momentum).E();
								}
							}
						
							TLorentzVector eCalMomentum;
						
							// Now loop over all ECAL grids and sum their 4-momenta
							for(int eta = 0; eta < granularity; ++eta)
							{
								TLorentzVector	gridMomentum;
								for(int phi = 0; phi < granularity; ++phi)
								{
									const Double_t newEta = etaMin + (0.5 + eta)*etaWidth;
									gridMomentum.SetPtEtaPhiM(eCalGrid[eta][phi]/cosh(newEta), newEta, phiMin + (0.5 + phi)*phiWidth, 0.);
									eCalMomentum += gridMomentum;
								}
							}

							debugOut << "eta/phi/Ehad/Eem:   " << (tower->Position).Eta() << "\t" << (tower->Position).Phi() << "\t" << (tower->Ehad) << "\t" << (tower->Eem) << "\n";
						
							// Now use the total 4-momentum of the ECAL to define the 
							// mass and direction of the HCAL
							// Store the original momentum in the area field
							(tower->Area) = (tower->Momentum);
							(tower->Momentum).SetPtEtaPhiM((tower->Momentum).Pt(), eCalMomentum.Eta(), eCalMomentum.Phi(), (tower->Momentum).E()/(eCalMomentum.E()/eCalMomentum.M()));
						
							delete fItTowerConstituents;
						}
					}
				}
				
				// Now that we've increased the angular resolution of the towers, 
				// re-cluster the constituents
				
				// Now feed all constituents into fastjet
				{
					inputList.clear();
					Candidate const* constituent;
					int numConstituents = 0;
		
					itJetConstituents->Reset();
					while((constituent = static_cast<Candidate const*>(itJetConstituents->Next())))
					{
						const TLorentzVector& constituentP4 = constituent->Momentum;
						inputList.emplace_back(constituentP4.Px(), constituentP4.Py(), constituentP4.Pz(), constituentP4.E());
						inputList.back().set_user_index( numConstituents++ );
					}
				
					debugOut << "numConstituents:    " << numConstituents << endl;
					debugOut << "numGoodMuons:       " << goodMuons.size() << endl;
				}
				
				// re-cluster the constituents				
				fastjet::ClusterSequence clusterSequence(inputList, *fJetDefinition);
				std::vector<fastjet::PseudoJet> outputList = fastjet::sorted_by_pt(clusterSequence.inclusive_jets( jetPt * fMinCoreRatio ));
				
				if(not outputList.empty()) // Make sure at least one subjet passed the pt cut
				{
					fastjet::PseudoJet coreCopy = outputList.front(); // The highest pt jet is the core
										
					debugOut << "core Pt:            " << sqrt(coreCopy.kt2()) << endl;
					debugOut << "core mass:          " << coreCopy.m() << endl;
					
					{
						// Find the number of core constituents
						std::vector<fastjet::PseudoJet> coreConstituents = clusterSequence.constituents(outputList.front());
						
						debugOut << "coreConstituents:   " << coreConstituents.size() << endl;
						
						// Find any muons in the core and subtract their 4-momentum
						for(std::vector<fastjet::PseudoJet>::iterator itCore = coreConstituents.begin(); itCore != coreConstituents.end(); ++itCore)
						{
							Candidate const* const coreConstituent = static_cast<Candidate*>(jetConstituents->At(itCore->user_index()));
						
							if(abs(coreConstituent->PID) == 13)
								coreCopy -= *itCore;
						}
					}
					
					TLorentzVector coreMomentum(coreCopy.px(), coreCopy.py(), coreCopy.pz(), coreCopy.E()); //xyzt
						
					// Sometimes there will be more than 1 muon. In those cases, we will
					// assume that the muons were emitted from highest to lowest pt.
											
					// Sort muons low to high
					std::sort(goodMuons.begin(), goodMuons.end(), SortCandidatePt_Low2High);
					
					debugOut << "\n";
					
					// Now add back each muon's 4-momentum, estimate its neutrino, and calculate the resulting boost invariant
					for(std::vector<Candidate const*>::iterator itMuon = goodMuons.begin(); itMuon != goodMuons.end(); ++itMuon)
					{
						TVector3 muonP3;
						{
							const TLorentzVector& muonMomentum = (*itMuon)->Momentum;
							muonP3 = muonMomentum.Vect();
						
							// Add the muon back to the core
							coreMomentum += muonMomentum;
							
							debugOut << "muon Pt:            " << muonMomentum.Pt() << endl;
						}
						
						TVector3 neutrino;
						{
							// We do not want to bias the core momentum towards the muon emission direction
							// Thus, we give the neutrino the direction of the core, projecting out the 
							// magnitude of the muon momentum.
							//     neutrinoP3> = (muonP3>.coreP3> / coreP3>.coreP3>) * coreP3>
							const TVector3 coreP3 = coreMomentum.Vect();
							neutrino = coreP3;							
							neutrino *= coreP3.Dot(muonP3)/coreP3.Mag2();
						}
						
						{
							TLorentzVector neutrinoMomentum(neutrino.Px(), neutrino.Py(), neutrino.Pz(), neutrino.Mag());
							
							// Add the neutrino to the core							
							coreMomentum += neutrinoMomentum;
							// Also add the neutrino the jet;
							jet->Momentum += neutrinoMomentum;
							
							debugOut << "neutrino Pt:        " << neutrinoMomentum.Pt() << endl;
						}
						
						// Find the boost of the core, disallowing any mass lower than (fMinCoreMass)
						const Double_t coreBoostMass = fMinCoreMass;//fMinCoreMass std::max(coreMomentum.M(), fMinCoreMass);
						const Double_t coreBoost = coreMomentum.E()/coreBoostMass;						
													
						debugOut << "core Pt:            " << coreMomentum.Pt() << endl;
						debugOut << "core Mass:          " << coreBoostMass << endl;							
						
						const TVector3 coreP3 = coreMomentum.Vect();
						const Double_t CMInvariant = coreBoost * (muonP3.Cross(coreP3)).Mag() / muonP3.Dot(coreP3);
						
						if(CMInvariant <= fMaxEmissionInvariant)
						{
							// Congratulations, you've just tagged a b-jet
						}
						
						debugOut << "CM invaraint:       " << CMInvariant << endl;
						debugOut << "\n";
						
					}// End loop over muons
					
					cout << debugOut.str();
					cout << "\n---------------------------------------------------------\n";				
				}// End core pt cut					
			}//End Muon pt cut			
			
			delete itJetConstituents;			
		}//End jet pt cut  	
	}//End loop over jets
}

*/

