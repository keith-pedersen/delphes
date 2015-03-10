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
int flavor(int PID)
{
	PID = abs(PID);
	
	if(PID > 1000 && PID < 10000)
		PID /= 1000;
	else if(PID > 100)
		PID = (PID%1000)/100;
	
	return PID;
}

//------------------------------------------------------------------------------

HighPtBTagger::HighPtBTagger():
	fJetInputArray(0), fItJetInputArray(0), fCoreDefinition(0),
	histoFile(0) {}

//------------------------------------------------------------------------------

HighPtBTagger::~HighPtBTagger()
{ /* No memory to clean-up */}

//------------------------------------------------------------------------------

void HighPtBTagger::Init()
{
 	// read parameters

	fBitNumber = GetInt("BitNumber", 3);
	
	fMinJetPt = GetDouble("MinJetPt", 1000.);
	fMinMuonPt = GetDouble("MinMuonPt", 15.);
	
	fMinTowerPtRatio = GetDouble("MinTowerPtRatio", .01);
	fCoreAntiktR = GetDouble("CoreAntiktR", .1);
	fMinCorePtRatio = GetDouble("MinCorePtRatio", .8);
	
	fCoreMassHypothesis = GetDouble("CoreMassHypothesis", 1.9);
	fCoreMassHypothesis2 = fCoreMassHypothesis * fCoreMassHypothesis;
	fMinFinalMass = GetDouble("MinFinalMass", 5.3);
	fMaxFinalMass = GetDouble("MaxFinalMass", fMinFinalMass*2);	
	
	{
		const Double_t sigma = abs(GetDouble("PercentSolidAngle", .9));
		
		if(sigma >= 1.)
		{	
			throw runtime_error("HightPtBTagger. Supplied sigma is out of bounds (must be in [0,1) )");
		}
		
		fMaxEmissionInvariant = 1./tan(sqrt(1. - sigma));
		
		cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Emission Invariant: " << fMaxEmissionInvariant << endl;
	}
	
	fCoreDefinition = new fastjet::JetDefinition(fastjet::antikt_algorithm, fCoreAntiktR);
	
	// import input array(s)
	fJetInputArray = ImportArray(GetString("JetInputArray", "FastJetFinder/jets"));
	fItJetInputArray = fJetInputArray->MakeIterator();
	
	fAllParticles = ImportArray("Delphes/allParticles");
	
	// Restore current gDirectory after done initializing. 
	//TDirectory::TContext ctxt(0);
	
	// Append to existing root file
	//histoFile = new TFile("HighPtBTagger_histos.root", "recreate");
	
	const Double_t
		ptMax_jet = 2000.,
		ptWidth_jet = 20.,
		
		ptMax_muon = 1000.,
		ptWidth_muon = 5.,
		
		massMin = 1e-1, 
		massMax = 100.,
		
		invariantMin = 1e-2,
		invariantMax = 10.*fMaxEmissionInvariant,
		
		angleMin = 4e-4,
		angleMax = .4;
	
	const Int_t 
		numBins_jetPt = (ptMax_jet - fMinJetPt)/ptWidth_jet,
		numBins_muonPt = (ptMax_muon - fMinMuonPt)/ptWidth_muon,
		numBins_mass = 100,
		numBins_invariant = 100,
		numBins_angle = 50;
		
	// Array of lower edges (+ total upper edge) for mass and invariant
	// This allows an intrinsic log scale for the binning
	Double_t 
		invariantEdges[numBins_invariant + 1],
		massEdges[numBins_mass + 1],
		angleEdges[numBins_angle +1];
		
	// Fill the lower edges
	{	
		Double_t lowerEdge = log(invariantMin);	
		const Double_t logStep = (log(invariantMax) - lowerEdge)/numBins_invariant;
		const Double_t* lastBin = invariantEdges + numBins_invariant;
		
		for(Double_t* bin = invariantEdges; bin <= lastBin; ++bin)
		{
			*bin = exp(lowerEdge);
			lowerEdge += logStep; // Cumulative rounding error, but very small
		}
	}
	
	{	
		Double_t lowerEdge = log(massMin);	
		const Double_t logStep = (log(massMax) - lowerEdge)/numBins_mass;
		const Double_t* lastBin = massEdges + numBins_mass;
		
		for(Double_t* bin = massEdges; bin <= lastBin; ++bin)
		{
			*bin = exp(lowerEdge);
			lowerEdge += logStep; // Cumulative rounding error, but very small
		}
	}
	
	{	
		Double_t lowerEdge = log(angleMin);	
		const Double_t logStep = (log(angleMax) - lowerEdge)/numBins_angle;
		const Double_t* lastBin = angleEdges + numBins_angle;
		
		for(Double_t* bin = angleEdges; bin <= lastBin; ++bin)
		{
			*bin = exp(lowerEdge);
			lowerEdge += logStep; // Cumulative rounding error, but very small
		}
	}
	
	// The following histograms will be automatically owned by histoFile.
	// Thus we need not delete them, only histoFile.
	
	// create ID
	
	pt_Jets = new TH1I("pt_Jets", "Jets vs. pt_Jet (GeV)",
		numBins_jetPt, fMinJetPt, ptMax_jet);
	pt_Jets->Sumw2();
	
	pt_Muons = new TH1I("pt_Muons", "Muons (in jets above threshold) vs. pt_Muon (GeV)",
		numBins_muonPt, fMinMuonPt, ptMax_muon);
	
	
	pt_MuJets = new TH1F("pt_MuJets", "% MuJet (Jet with muons) vs. pt_Jet (GeV)",
		numBins_jetPt, fMinJetPt, ptMax_jet);
	pt_MuJets->Sumw2();
	
	coreRatioHisto = new TH1F("coreRatio", "% Cores vs. (pt_Core/pt_MuJet)",
		100, 0., 1.);
	coreRatioHisto->Sumw2();
		
	coreMass_Raw = new TH1F("coreMass_Raw", "Raw HardCore (w/ muons subtracted) vs. Their mass (GeV)",
		numBins_mass, massEdges);
		
	deltaMass_Raw = new TH1F("deltaMass_Raw", "Raw HardCores vs. Their deltaMass from adding muon (GeV)",
		numBins_mass, massEdges);
		
	deltaMass_Fixed = new TH1F("deltaMass_Fixed", "Fixed HardCores vs. Their deltaMass from adding muon (GeV)",
		numBins_mass, massEdges);
		
	pt_TaggedJets_Raw = new TH1F("pt_TaggedJets_Raw", "% Jets Tagged with Raw Core Mass (wrt MuJets) vs. pt_Jet (GeV)",
		numBins_jetPt, fMinJetPt, ptMax_jet);
	pt_TaggedJets_Raw->Sumw2();
	
	pt_TaggedJets_Fixed = new TH1F("pt_TaggedJets_Fixed", "% Jets Tagged with Fixed Core Mass (wrt MuJets) vs. pt_Jet (GeV)",
		numBins_jetPt, fMinJetPt, ptMax_jet);
	pt_TaggedJets_Fixed ->Sumw2();
	
	
	xTrueHisto = new TH1F("xTrue", "Muons vs. x(muon, mother)",
		numBins_invariant, invariantEdges);
	
	xRawHisto = new TH1F("x(mu, HardCore_Raw)", "Muons vs. x(muon, HardCore_Raw)",
		numBins_invariant, invariantEdges);
		
	xFixedHisto = new TH1F("x(mu, HardCore_Fixed)", "Muons vs. x(muon, HardCore_Fixed)",
		numBins_invariant, invariantEdges);
	
	xRawRatio = new TH1F("x_Raw/x_True", "Muons vs. x_Raw/x_True",
		numBins_invariant, invariantEdges);
	
	xFixedRatio = new TH1F("x_Fixed/x_True", "Muons vs. x_Fixed/x_True",
		numBins_invariant, invariantEdges);
		
 	trueCoreError = new TH1F("trueCoreError", "difference in angle from core to actual mother",
		numBins_angle, angleEdges);
		
	trueCoreErrorMuonTwice = new TH1F("trueCoreErrorMuonTwicer", "difference in angle from core to actual mother\n(when muon is added twice)",
		numBins_angle, angleEdges);
	
	// create 2D
	
	coreMass_Raw_X_coreSpread = new TH2F("coreMass_Raw vs. coreSpread", "coreMass_Raw vs. coreSpread (rel. std. dev. of constituent E)",
		numBins_mass, massEdges,
		100, 0., 5.);
	
	deltaMass_Raw_X_coreMass_Raw = new TH2F("deltaMass_Raw vs. coreMass_Raw", "deltaMass_Raw (coreAfterMuon - coreBeforeMuon) vs. coreMass_Raw",
		numBins_mass, massEdges,
		numBins_mass, massEdges);
	
	xRaw_X_coreRatio = new TH2F("x(mu, HardCore_Raw) vs. coreRatio", "Muons vs. x(mu, HardCore_Raw) vs. coreRatio",
		numBins_invariant, invariantEdges, 
		100, 0., 1.);
		
	xFixed_X_coreRatio = new TH2F("x(mu, HardCore_Fixed) vs. coreRatio", "Muons vs. x(mu, HardCore_Fixed) vs. coreRatio",
		numBins_invariant, invariantEdges, 
		100, 0., 1.);
	 
	
	/*
	Not needed for now, just do it right before we divide
			
	// Tell all the histograms to keep track of squared weights. This will help
	// for those we decide to divide
	TIterator activeHistograms(histoFile->GetList()->MakeIterator());
	TH1* histogram;
	
	while((histogram = static_cast<TH1*>(activeHistograms.Next())))
	{
		histogram->Sumw2();
	}
	*/
}

//------------------------------------------------------------------------------

void HighPtBTagger::Finish()
{
	// Turn the percent histos into percent histos
	pt_TaggedJets_Raw->Divide(pt_MuJets);
	pt_TaggedJets_Fixed->Divide(pt_MuJets);
	
	pt_MuJets->Divide(pt_Jets);

	delete fCoreDefinition;
	delete fItJetInputArray;
	//histoFile->Write();
	delete histoFile; // Close() is first line of dtor
}

//------------------------------------------------------------------------------

void HighPtBTagger::Process()
{
	// loop over input jets
	Candidate* jet;
	fItJetInputArray->Reset();
	
	while((jet = static_cast<Candidate*>(fItJetInputArray->Next())))
	{
		// At best, neutrino estimation can only double the pt of the jet, 
		// so only look at jets with at least half the required final pt
		if((jet->Momentum.Pt()) >= .5*fMinJetPt) 
		{
			// We'll sort jet constituents into two categories (where goodMuons are those which pass the pt cut)
			std::vector<Candidate const*> goodMuons, everythingElse;
			
			{
				TObjArray const* const jetConstituents = jet->GetCandidates();
				TIterator* const itJetConstituents = jetConstituents->MakeIterator();
				Candidate const* constituent;
				
				everythingElse.reserve(jetConstituents->GetEntriesFast());
				
				// Find all muons that pass the pt cuts
				while((constituent = static_cast<Candidate const*>(itJetConstituents->Next())))
				{
					if(abs(constituent->PID) == 13)
					{
						const Double_t muonPt = (constituent->Momentum).Pt();
					
						if(muonPt >= fMinMuonPt)
						{
							pt_Muons->Fill(muonPt);
							goodMuons.push_back(constituent);
							continue;
						}
					}
					everythingElse.push_back(constituent);
				}
				
				delete itJetConstituents;
			}
			
			if(not goodMuons.empty()) // Ensure we have at least one good muon
			{
				// Sometimes there will be more than 1 muon. In those cases, we will
				// assume that the muons were emitted from highest to lowest pt.
				// This means we'll need to add them back from low to high pt, so that 
				// the CM frame of the muon emitted 2nd doesn't include the 1st muon.
									
				// Sort muons low to high 
				std::sort(goodMuons.begin(), goodMuons.end(), SortCandidatePt_Low2High);
			
				// To find the hard core, we recluster the jet constituents
				vector<fastjet::PseudoJet> reclusterInput;
				
				// Add the good muons to the reclusterInput
				// Give them a negative user_index ( user_index = index in goodMuons - goodMuons.size() )
				// This makes them easy to find after the reclustering.
				{
					const int numMuons = goodMuons.size();
					for(int iMu = 0; iMu < numMuons; ++iMu)
					{
						const TLorentzVector& muonMomentum = goodMuons[iMu]->Momentum;
						reclusterInput.emplace_back(muonMomentum.Px(), muonMomentum.Py(), muonMomentum.Pz(), muonMomentum.E());
						reclusterInput.back().set_user_index(iMu - numMuons);
					}
				}
				
				// Now put the jet constituents into reclusterInput
				// Use all tracks, but only use towers/eFlowNeutrals (charge == 0) if they are above the relative pt threshold
				{
					const Double_t minTowerPt = fMinTowerPtRatio*(jet->Momentum).Pt();
					for(unsigned int iEverythingElse = 0; iEverythingElse < everythingElse.size(); ++iEverythingElse)
					{
						Candidate const* const constituent = everythingElse[iEverythingElse];
						
						if((constituent->Charge == 0) && ((constituent->Momentum).Pt() < minTowerPt))
							continue;
						
						reclusterInput.emplace_back((constituent->Momentum).Px(), (constituent->Momentum).Py(), (constituent->Momentum).Pz(), (constituent->Momentum).E());
						reclusterInput.back().set_user_index(iEverythingElse);
					}
				}
				
				// We'll be working with 2 cores, one with a raw mass, the other
				// where the mass has been fixed
				TLorentzVector coreRaw, coreFixed;				
				
				Double_t 
					minRawX = 99.,
					minFixedX = 99.;
				// Now recluster
				{
					fastjet::ClusterSequence recluster(reclusterInput, *fCoreDefinition);
					std::vector<fastjet::PseudoJet> subJets = fastjet::sorted_by_pt(recluster.inclusive_jets(20.));
					
					// Old choice, 
					fastjet::PseudoJet coreBig = subJets.front(); 
					
					/*
					fastjet::PseudoJet core;					
					
					// The core should be the subjet which will intrinsically minimize the CM mass
					// Thus, we want to minimize massFactor: sqrt(pt_core)*deltaR(muon, core)
					{
						const fastjet::PseudoJet& hardestMuon = reclusterInput[goodMuons.size() - 1];
						Double_t minMassFactor2 = 9e9;
						unsigned int coreIndex = 0;
						
						for(unsigned int iSub = 0; iSub < subJets.size(); ++iSub)
						{
							// By minimizing the squared massFactor, we'll minimize the massFactor
							const Double_t massFactor2 = subJets[iSub].kt2()*hardestMuon.squared_distance(subJets[iSub]);
							
							if(massFactor2 < minMassFactor2)
							{
								minMassFactor2 = massFactor2;
								coreIndex = iSub;
							}
						}
						
						core = subJets[coreIndex];
					}
					
					*/
												
					// Find the core's constituents
					std::vector<fastjet::PseudoJet> coreConstituents = recluster.constituents(subJets.front());
					
					{
						// Find any good muons in the core and subtract their 4-momentum
						// For the rest of the core constituents, keep a running tally of their energy and energy**2
						Double_t 
							sumCoreE = 0.,
							sumCoreE2 = 0.;
						unsigned int muonsInCore = 0;
						for(std::vector<fastjet::PseudoJet>::const_iterator itCore = coreConstituents.begin(); itCore != coreConstituents.end(); ++itCore)
						{
							if(itCore->user_index() < 0)
							{
								core -= *itCore;
								++muonsInCore;
							}
							else
							{
								const Double_t thisE = itCore->E();
								sumCoreE += thisE;
								sumCoreE2 += thisE*thisE;
							}
						}
					
						coreRaw.SetPxPyPzE(core.px(), core.py(), core.pz(), core.E());
					
						// Now assume that the reconstructed hadronic mass is horribly incorrect,
						// created almost entirely by the geometry/granularity of the Cal cells, 
						// as well as the magnetic spreading in phi
						// Thus, apply a mass hypothesis to the core.
						{
							TVector3 coreP3 = coreRaw.Vect();
							Double_t coreE = coreRaw.E();
							// Rescale the momentum to the new mass
							coreP3 *= sqrt(coreE*coreE - fCoreMassHypothesis2) / coreP3.Mag();

							coreFixed.SetPxPyPzE(coreP3.Px(), coreP3.Py(), coreP3.Pz(), coreE);
						}
					
						// Study the raw core mass, now that the muons have been subtracted
						{
							const Double_t coreMassRaw = coreRaw.M();
							coreMass_Raw->Fill(coreMassRaw);
							// Here, coreSpread is the relative standard deviation of the energy (sqrt(N * sum(E**2)/sum(E)**2  - 1))
							coreMass_Raw_X_coreSpread->Fill(coreMassRaw, sqrt((coreConstituents.size() - muonsInCore)*sumCoreE2/(sumCoreE*sumCoreE) - 1.));
						}
					}
										
					
					
					// Now add back each muon's 4-momentum, estimate its neutrino, and calculate the resulting boost invariant
					for(std::vector<Candidate const*>::const_iterator itMuon = goodMuons.begin(); itMuon != goodMuons.end(); ++itMuon)
					{
						const Double_t 
							originalMassRaw = coreRaw.M(),
							originalMassFixed = coreFixed.M();
						
						const TVector3 muonP3 = ((*itMuon)->Momentum).Vect();
						
						// Add the muon back to the core
						coreRaw += (*itMuon)->Momentum;
						coreFixed += (*itMuon)->Momentum;
						
						Candidate const* matriarch = static_cast<Candidate*>(fAllParticles->At((*itMuon)->M1));
						{
							// We are interested in the original boosted mother (the matriarch). This muon
							// could have come from a tau which came from something else.
							// We need to look back all the way to the start of hadronization.
						
							// To look back all the way to hadronization, we find the first mother with more than 1 mother (when color charge still existed)
							Candidate const* mother = matriarch;
							while(mother->M2 == 0)
							{
								const int absMotherPID = abs(mother->PID);
							
								if((absMotherPID >= 22) or (absMotherPID <= 24))
								{
									// Safety check for primary muons (A, Z, W)
									// Pythia hadronic decay never specifically invokes a W; hence,
									// if this is a W, it must have come from the hard interaction.
									// By breaking here, we treat the A/Z/W as the original mother
									break;
								}
								else
								{
									matriarch = mother;
									mother = static_cast<Candidate*>(fAllParticles->At(mother->M1));
								
									if(not mother)
										break;
								}
							}
						}
						
						const TVector3 motherP3 = (matriarch->Momentum).Vect();
						
						const Double_t xTrue = (matriarch->Momentum).E() / matriarch->Mass * 
								(muonP3.Cross(motherP3).Mag() / muonP3.Dot(motherP3));
						
						// Study difference to core for 
						{
							const int matriarchFlavor = flavor(matriarch->PID);
							
							if((matriarchFlavor > 3) and (matriarchFlavor < 6))
							{
								xTrueHisto->Fill(xTrue);						
								trueCoreError->Fill(motherP3.Angle(coreFixed.Vect()));
						
								// Does adding the muon twice get us closer to the actual mother?
								{
									TLorentzVector muonTwice(coreFixed);
									muonTwice += (*itMuon)->Momentum;
						
									trueCoreErrorMuonTwice->Fill(motherP3.Angle(muonTwice.Vect()));
								}	
								// Yes, it does. This helps us in two regards. The core is closer to
								// the mother, and the muon is closer to the core. 
							}
						}						
						
						// Simulate the neutrino by adding the muon twice
						coreRaw += (*itMuon)->Momentum;
						coreFixed += (*itMuon)->Momentum;	
						jet->Momentum += (*itMuon)->Momentum;
						
						/*						
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
						
						// Study original mass of core, vs. mass when muon (and neutrino) are added
						{
							const Double_t deltaMassRaw = coreRaw.M() - originalMassRaw;
							deltaMass_Raw->Fill(deltaMassRaw);
							deltaMass_Raw_X_coreMass_Raw->Fill(deltaMassRaw, originalMassRaw);
							
							deltaMass_Fixed->Fill(coreFixed.M() - originalMassFixed);
						}
						
						// Find the boost of the core, forcing the mass to be in [fMinFinalMass, fMaxFinalMass]
						
						const Double_t 
							boostMassRaw = std::min(std::max(coreRaw.M(), fMinFinalMass), fMaxFinalMass),
							boostMassFixed = std::min(std::max(coreFixed.M(), fMinFinalMass), fMaxFinalMass);
							
						const Double_t 
							coreBoostRaw = coreRaw.E()/boostMassRaw,
							coreBoostFixed = coreFixed.E()/boostMassFixed;
														
						const Double_t 
							xRaw = coreBoostRaw * (muonP3.Cross(coreRaw.Vect()).Mag() / muonP3.Dot(coreRaw.Vect())),
							xFixed = coreBoostFixed * (muonP3.Cross(coreFixed.Vect()).Mag() / muonP3.Dot(coreFixed.Vect()));
							
						xRawHisto->Fill(xRaw);
						xFixedHisto->Fill(xFixed);
						
						xRawRatio->Fill(xRaw/xTrue);
						xFixedRatio->Fill(xFixed/xTrue);
					
						minRawX = std::min(minRawX, xRaw);
						minFixedX = std::min(minFixedX, xFixed);
					}// End loop over muons
				}// End Reclustering
				
				const Double_t jetPt = (jet->Momentum).Pt();
				pt_MuJets->Fill(jetPt);		
				
				const Double_t coreRatio = (coreRaw.Pt() / jetPt);
				coreRatioHisto->Fill(coreRatio);
				
				xRaw_X_coreRatio->Fill(minRawX, coreRatio);
				xFixed_X_coreRatio->Fill(minFixedX, coreRatio);
				
					
				// Make sure it's a hard core
				if(coreRatio >= fMinCorePtRatio)
				{
					if(minRawX <= fMaxEmissionInvariant)
					{
						jet->BTag |= (1 << fBitNumber);
						pt_TaggedJets_Raw->Fill(jetPt);
					}
					
					if(minFixedX < fMaxEmissionInvariant)
					{
						jet->BTag |= (1 << fBitNumber);
						pt_TaggedJets_Fixed->Fill(jetPt);
					}
				}
			}//End muon pt cut
			
			pt_Jets->Fill((jet->Momentum).Pt());			
		}//End jet initial pt cut  	
	}//End loop over jets
}

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
