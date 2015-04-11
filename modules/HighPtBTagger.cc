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
	fJetInputArray(0), fItJetInputArray(0), fAllParticles(0), fCoreDefinition(0),
	histoDir(0) {}

//------------------------------------------------------------------------------

HighPtBTagger::~HighPtBTagger()
{ /* No memory to clean-up */}

//------------------------------------------------------------------------------

// Fill array[numBins+1] from [minValue, maxValue] with lower edges (plus the upper edge of the last bin)
void FillLogEdges(Double_t* const array, const int numBins, const Double_t minValue, const Double_t maxValue)
{
	Double_t lowerEdge = log(minValue);
	const Double_t logStep = (log(maxValue) - lowerEdge)/numBins;
	Double_t* const lastBin = array + numBins;

	for(Double_t* bin = array; bin <= lastBin; ++bin)
	{
		*bin = exp(lowerEdge);
		lowerEdge += logStep; // Cumulative rounding error, but very small
	}
}

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

	// A form suggested by WK (in "How Mindless") is accurate over all domains
	const TVector3
		v1 = one*two.Mag(),
		v2 = two*one.Mag();

	return 2*atan(sqrt((v1-v2).Mag2()/(v1+v2).Mag2()));
}


// This templated function allows new root histograms to be created with formatted suffixes (damn those C-strings, can't add em)
template<class TH, typename... Args>
TH* NewTH(const std::string& label, const std::string& namePrefix, const std::string& titlePrefix, Args&&... args)
{
	return new TH((namePrefix + "_" + label).c_str(), (titlePrefix + " (" + label + ")").c_str(), std::forward<Args>(args)...);
}

// This templated function allows a set of root histograms (of same type and bin specs,
// but with different labels and directories)
template<class TH, typename... Args>
std::vector<TH*> NewTHVec(const std::vector<std::string>& labelVec, const std::vector<TDirectory*>& dirVec,
	const std::string& namePrefix, const std::string& titlePrefix, Args&&... args)
{
	std::vector<TH*> returnVec;

	for(unsigned int i = 0; i < dirVec.size(); ++i)
	{
		dirVec[i]->cd();
		returnVec.push_back(NewTH<TH>(labelVec[i], namePrefix, titlePrefix, std::forward<Args>(args)...));
	}

	return returnVec;
}

void HighPtBTagger::FillTH1F(const std::vector<std::vector<TH1F*> >& histoVec, const Int_t coreIndex, const Int_t flavor, const Double_t fillValue)
{
	// Fill the specific and general histograms
	histoVec[coreIndex][0]->Fill(fillValue);
	histoVec[coreIndex][flavorToHistoIndex[flavor]]->Fill(fillValue);
}

void HighPtBTagger::FillTH2F(const std::vector<std::vector<TH2F*> >& histoVec, const Int_t coreIndex, const Int_t flavor,
	const Double_t fillValueX, const Double_t fillValueY)
{
	// Fill the specific and general histograms
	histoVec[coreIndex][0]->Fill(fillValueX, fillValueY);
	histoVec[coreIndex][flavorToHistoIndex[flavor]]->Fill(fillValueX, fillValueY);
}

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

		coreRatioMin = 0.,
		coreRatioMax = 1.,

		angleMin = 4e-4,
		angleMax = .4;

	const Int_t
		numBins_jetPt = (ptMax_jet - fMinJetPt)/ptWidth_jet,
		numBins_muonPt = (ptMax_muon - fMinMuonPt)/ptWidth_muon,
		numBins_mass = 100,
		numBins_invariant = 100,
		numBins_coreRatio = 50,
		numBins_angle = 50;

	// Array of lower edges (+ total upper edge) for mass and invariant
	// This allows an intrinsic log scale for the binning
	Double_t
		invariantEdges[numBins_invariant + 1],
		massEdges[numBins_mass + 1],
		angleEdges[numBins_angle +1];

	FillLogEdges(invariantEdges, numBins_invariant, invariantMin, invariantMax);
	FillLogEdges(massEdges, numBins_mass, massMin, massMax);
	FillLogEdges(angleEdges, numBins_angle, angleMin, angleMax);

	// The following histograms will be automatically owned by the current
	// gDirectory. So as long as that diretory is deleted, these histograms
	// will be as well

	// We'll have 3 sets of histos; all, light, and heavy
	const int
		allIndex = 0,
		lightIndex = 1,
		heavyIndex = 2;

	// Map each quark flavor to one of those histos
	flavorToHistoIndex[0] = lightIndex;
	flavorToHistoIndex[1] = lightIndex;
	flavorToHistoIndex[2] = lightIndex;
	flavorToHistoIndex[3] = lightIndex;
	flavorToHistoIndex[4] = heavyIndex;
	flavorToHistoIndex[5] = heavyIndex;
	flavorToHistoIndex[6] = heavyIndex;

	// Create a vector of labels (suffixes) for the histograms
	const std::vector<std::string> labels = {"all", "light", "heavy"};

	{
		// We'll keep the histograms inside a special sub-directory
		histoDir = gDirectory->mkdir("HighPtBTagger_Histograms");
		TDirectory::TContext cdHisto(histoDir);

		// Create 1D Histograms

		// TH1I

		pt_Jets = new TH1I("pt_Jets", "Jets vs. pt_Jet (GeV)",
			numBins_jetPt, fMinJetPt, ptMax_jet);
		pt_Jets->Sumw2();

		pt_Muons = new TH1I("pt_Muons", "Muons (in jets above threshold) vs. pt_Muon (GeV)",
			numBins_muonPt, fMinMuonPt, ptMax_muon);


		// motherPID will be created and filled during Finish().
		// Reset the PID map, however
		pdgID = std::map<int, int>();

		{
			TDirectory* x_TrueDir = gDirectory->mkdir("xTrue");
			// Change to the x_True directory. When the brackets close, we'll return
			TDirectory::TContext cd_xTrue(x_TrueDir);

			for(std::vector<string>::const_iterator itLabel = labels.begin(); itLabel not_eq labels.end(); ++itLabel)
				x_True.push_back(NewTH<TH1F>(*itLabel, "x_True", "Muons vs. x(muon, boosted matriarch)", numBins_invariant, invariantEdges));
		}

		{
			TDirectory* pt_muJetsDir = gDirectory->mkdir("pt_muJets");
			TDirectory::TContext cd_xTrue(pt_muJetsDir);

			for(std::vector<string>::const_iterator itLabel = labels.begin(); itLabel not_eq labels.end(); ++itLabel)
			{
				pt_MuJets.push_back(NewTH<TH1F>(*itLabel, "pt_MuJets", "% MuJet (Jet with muons) vs. pt_Jet (GeV)  Matriarch: ", numBins_jetPt, fMinJetPt, ptMax_jet));
				pt_MuJets.back()->Sumw2();
			}
		}

		// Create the HardCore and MinCore histograms in seperate directories;
		TDirectory* hardcoreDir = gDirectory->mkdir("HardCore");
		{
			// Change to the HardCore directory. When the brackets close, we'll return
			TDirectory::TContext cdHardCore(hardcoreDir);

			TDirectory
				*hardcoreLight = hardcoreDir->mkdir("Light"),
				*hardcoreHeavy = hardcoreDir->mkdir("Heavy");

			std::vector<TDirectory*> hardcoreDirVec = {hardcoreDir, hardcoreLight, hardcoreHeavy};

			coreRatio.push_back(
				NewTHVec<TH1F>(labels, hardcoreDirVec, "coreRatio_HardCore", "HardCores vs. (pt_HardCore/pt_MuJet) Matriarch: ",
				numBins_coreRatio, coreRatioMin, coreRatioMax));

			mass.push_back(
				NewTHVec<TH1F>(labels, hardcoreDirVec, "mass_HardCore", "HardCores (w/ muons subtracted) vs. Their mass (GeV) Matriarch: ",
				numBins_mass, massEdges));

			deltaMass.push_back(
				NewTHVec<TH1F>(labels, hardcoreDirVec, "deltaMass_HardCore", "HardCores vs. Their deltaMass from adding muon (GeV)",
				numBins_mass, massEdges));

			x_Core.push_back(
				NewTHVec<TH1F>(labels, hardcoreDirVec, "x_HardCore", "Muons vs. x(muon, HardCore)",
				numBins_invariant, invariantEdges));

			deltaTrue.push_back(
				NewTHVec<TH1F>(labels, hardcoreDirVec, "deltaTrue_HardCore", "Difference in angle from HardCore to actual boosted matriarch",
				numBins_angle, angleEdges));

			deltaTrue2Mu.push_back(
				NewTHVec<TH1F>(labels, hardcoreDirVec, "deltaTrue2Mu_HardCore", "Difference in angle from HardCore (after adding muon twice) to actual boosted matriarch",
				numBins_angle, angleEdges));

			pt_Tagged.push_back(
				NewTHVec<TH1F>(labels, hardcoreDirVec, "pt_Tagged_HardCore", "% MuJets (HardCore) Tagged vs. pt_MuJet (GeV)",
				numBins_jetPt, fMinJetPt, ptMax_jet));

			for(std::vector<TH1F*>::iterator itTagged = pt_Tagged.back().begin(); itTagged not_eq pt_Tagged.back().end(); ++itTagged)
				(*itTagged)->Sumw2();

			// 2D Histograms
			x_Core__vs__x_True.push_back(
				NewTHVec<TH2F>(labels, hardcoreDirVec, "x_HardCore__vs__x_True", "Muons vs. x(muon, HardCore) vs. x(muon, boosted matriarch)",
				numBins_invariant, invariantEdges,
				numBins_invariant, invariantEdges));

			x_Core__vs__coreRatio.push_back(
				NewTHVec<TH2F>(labels, hardcoreDirVec, "x_HardCore__vs__coreRatio", "Muons vs. x(muon, HardCore) vs. (pt_HardCore/pt_MuJet)",
				numBins_invariant, invariantEdges,
				numBins_coreRatio, coreRatioMin, coreRatioMax));
		}

		TDirectory* mincoreDir = gDirectory->mkdir("MinCore");
		{
			TDirectory::TContext cdMinCore(mincoreDir);

			TDirectory
				*mincoreLight = mincoreDir->mkdir("Light"),
				*mincoreHeavy = mincoreDir->mkdir("Heavy");

			std::vector<TDirectory*> mincoreDirVec = {mincoreDir, mincoreLight, mincoreHeavy};

			coreRatio.push_back(
				NewTHVec<TH1F>(labels, mincoreDirVec, "coreRatio_MinCore", "MinCores vs. (pt_MinCore/pt_MuJet) Matriarch: ",
				numBins_coreRatio, coreRatioMin, coreRatioMax));

			mass.push_back(
				NewTHVec<TH1F>(labels, mincoreDirVec, "mass_MinCore", "MinCores (w/ muons subtracted) vs. Their mass (GeV) Matriarch: ",
				numBins_mass, massEdges));

			deltaMass.push_back(
				NewTHVec<TH1F>(labels, mincoreDirVec, "deltaMass_MinCore", "MinCores vs. Their deltaMass from adding muon (GeV)",
				numBins_mass, massEdges));

			x_Core.push_back(
				NewTHVec<TH1F>(labels, mincoreDirVec, "x_MinCore", "Muons vs. x(muon, MinCore)",
				numBins_invariant, invariantEdges));

			deltaTrue.push_back(
				NewTHVec<TH1F>(labels, mincoreDirVec, "deltaTrue_MinCore", "Difference in angle from MinCore to actual boosted matriarch",
				numBins_angle, angleEdges));

			deltaTrue2Mu.push_back(
				NewTHVec<TH1F>(labels, mincoreDirVec, "deltaTrue2Mu_MinCore", "Difference in angle from MinCore (after adding muon twice) to actual boosted matriarch",
				numBins_angle, angleEdges));

			pt_Tagged.push_back(
				NewTHVec<TH1F>(labels, mincoreDirVec, "pt_Tagged_MinCore", "% MuJets (MinCore) Tagged vs. pt_MuJet (GeV)",
				numBins_jetPt, fMinJetPt, ptMax_jet));

			for(std::vector<TH1F*>::iterator itTagged = pt_Tagged.back().begin(); itTagged not_eq pt_Tagged.back().end(); ++itTagged)
				(*itTagged)->Sumw2();

			x_Core__vs__x_True.push_back(
				NewTHVec<TH2F>(labels, mincoreDirVec, "x_MinCore__vs__x_True", "Muons vs. x(muon, MinCore) vs. x(muon, boosted matriarch)",
				numBins_invariant, invariantEdges,
				numBins_invariant, invariantEdges));

			x_Core__vs__coreRatio.push_back(
				NewTHVec<TH2F>(labels, mincoreDirVec, "x_MinCore__vs__coreRatio", "Muons vs. x(muon, MinCore) vs. (pt_MinCore/pt_MuJet)",
				numBins_invariant, invariantEdges,
				numBins_coreRatio, coreRatioMin, coreRatioMax));
		}
	}
}

//------------------------------------------------------------------------------

void HighPtBTagger::FillMotherPIDHisto(const Int_t PID)
{
	++(pdgID[abs(PID)]);
}

//------------------------------------------------------------------------------

void HighPtBTagger::WriteMotherPIDHisto()
{
	TDirectory::TContext histo(histoDir);

	motherPID = new TH1I("motherPID", "Count of muons vs. their mother's abs(PDG ID)",
		pdgID.size(), 0., 1.);

	TAxis* const xAxis = motherPID->GetXaxis();

	int binNumber = 1; //  bin zero is underflow
	for(std::map<Int_t, Int_t>::const_iterator itPID = pdgID.begin(); itPID not_eq pdgID.end(); ++itPID)
	{
		motherPID->SetBinContent(binNumber, itPID->second);
		xAxis->SetBinLabel(binNumber, to_string(itPID->first).c_str());
		++binNumber;
	}
}

//------------------------------------------------------------------------------

void HighPtBTagger::Finish()
{
	WriteMotherPIDHisto();

	// Turn the percent histos into percent histos
	for(std::vector<std::vector<TH1F*> >::iterator itCore = pt_Tagged.begin(); itCore not_eq pt_Tagged.end(); ++itCore)
	{
		for(unsigned int iFlavor = 0; iFlavor < itCore->size(); ++iFlavor)
		{
			(*itCore)[iFlavor]->Divide(pt_MuJets[iFlavor]);
		}
	}

	for(unsigned int iFlavor = 0; iFlavor < pt_MuJets.size(); ++iFlavor)
	{
		pt_MuJets[iFlavor]->Divide(pt_Jets);
	}

	delete fCoreDefinition;
	delete fItJetInputArray;
	//histoFile->Write();
	//delete histoFile; // Close() is first line of dtor
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
		const Double_t jetOriginalPt = jet->Momentum.Pt();

		if(jetOriginalPt >= .5*fMinJetPt)
		{
			// We'll sort jet constituents into two categories (where goodMuons are those which pass the pt cut)
			std::vector<Candidate const*> goodMuons, everythingElse;

			// Fill goodMuons & everythingElse (CHECKED 3/24/15)
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
				// Sort muons low to high (CHECKED 3/24/15)
				{
					// Sometimes there will be more than 1 muon. In those cases, we will
					// assume that the muons were emitted from highest to lowest pt.
					// This means we'll need to add them back from low to high pt, so that
					// the CM frame of the muon emitted later doesn't include the already emitted muon.

					// Sort muons low to high
					std::sort(goodMuons.begin(), goodMuons.end(), SortCandidatePt_Low2High);
				}

				// To find the hard core, we recluster the jet constituents
				vector<fastjet::PseudoJet> reclusterInput;

				// Add the good muons to the reclusterInput (CHECKED 3/24/15)
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

				// Add jet constituents (all tracks towers passing pt cut) (CHECKED 3/24/15)
				{
					// Use all tracks, but only use towers/eFlowNeutrals (charge == 0) if they are above the relative pt threshold
					// This is because the angular resolution of towers is much poorer, and we don't
					// want soft radiation unduly influencing the direction of the core

					const Double_t minTowerPt = fMinTowerPtRatio*jetOriginalPt;
					for(unsigned int iEverythingElse = 0; iEverythingElse < everythingElse.size(); ++iEverythingElse)
					{
						Candidate const* const constituent = everythingElse[iEverythingElse];

						// neutral means tower
						if((constituent->Charge == 0) && ((constituent->Momentum).Pt() < minTowerPt))
							continue;

						reclusterInput.emplace_back((constituent->Momentum).Px(), (constituent->Momentum).Py(), (constituent->Momentum).Pz(), (constituent->Momentum).E());
						reclusterInput.back().set_user_index(iEverythingElse);
					}
				}

				// We'll be working with 2 cores. Keeping them in a vector allows us to re-use the code more easily
				//    index 0: the hardest subjet (will control the global histograms)
				//    index 1: the core which minimizes the reconstructed mass (when adding the muon), if different than 0
				std::vector<TLorentzVector> coreP4Vec;
				std::vector<Double_t> coreMinX = {9e99, 9e99};
				Int_t coreFlavor = 0; // The flavor of the highest energy matriarch

				{// Scope of ClusterSequence
					const Double_t minSubjetPt = 10.*fCoreMassHypothesis;

					// Recluster the jet, to find subjets, then ensure the pt_subjet > 10 * coreMassHypothesis
					fastjet::ClusterSequence recluster(reclusterInput, *fCoreDefinition);
					std::vector<fastjet::PseudoJet> subJets = fastjet::sorted_by_pt(recluster.inclusive_jets(minSubjetPt));

					// Make sure at least one subjet passed the cut
					if(not subJets.empty())
					{
						std::vector<Candidate const*> goodMuonsMatriarch;
						std::vector<Int_t> goodMuonsMotherFlavor;
						// Find the muon matriarch and fill the mother histogram (CHECKED 3/26/15)
						{
							// * The "matriarch" is the original boosted mother
							//      We need to look back all the way to the start of hadronization.
							// * The "mother" is the particle which emitted the muon
							//      Only if the real mother is a tau do we consider the grandmother the mother
							for(std::vector<Candidate const*>::const_iterator itMuon = goodMuons.begin(); itMuon != goodMuons.end(); ++itMuon)
							{
								bool noMotherYet = true;

								Candidate const* matriarch;
								{
									// To look back all the way to hadronization, we find the first matriarch with more than 1 mother
									// (more than 1 mother means color charge still existed)
									Candidate const* possibleMatriarch = static_cast<Candidate*>(fAllParticles->At((*itMuon)->M1));

									do
									{
										const int absPossibleMatriarchPID = abs(possibleMatriarch->PID);
										matriarch = possibleMatriarch;

										if(noMotherYet && (absPossibleMatriarchPID not_eq 15))
										{
											noMotherYet = false;
											FillMotherPIDHisto(absPossibleMatriarchPID);
											goodMuonsMotherFlavor.push_back(flavor(absPossibleMatriarchPID));
										}

										if((possibleMatriarch->M2 not_eq 0) or
											((absPossibleMatriarchPID >= 22) and (absPossibleMatriarchPID <= 24)))
										{
											// Safety check for primary muons (A, Z, W)
											// Pythia hadronic decay never specifically invokes a W; hence,
											// if this is a W, it must have come from the hard interaction.
											// By breaking here, we treat the A/Z/W as the original mother
											break;
										}
										else
										{
											possibleMatriarch = static_cast<Candidate*>(fAllParticles->At(possibleMatriarch->M1));
										}
									}while(possibleMatriarch);

									const Int_t matriarchFlavor = flavor(matriarch->PID);

									if(noMotherYet)
									{
										cout << "MotherFail" << endl;
										goodMuonsMotherFlavor.push_back(matriarchFlavor);
									}
									goodMuonsMatriarch.push_back(matriarch);

									coreFlavor = std::max(matriarchFlavor, coreFlavor);
								}
							}
						}

						// Keep pseudojets in a parallel vector to coreP4Vec
						std::vector<fastjet::PseudoJet> coreVec;

						// The hardCore is the core with the highest pt
						coreVec.push_back(subJets.front());

						// Bin the mass of the HardCore
						FillTH1F(mass, 0, coreFlavor, subJets.front().m());

						// Get fastjets internal jets (because it couldn't alter our input jets,
						// so they have no internal information). I guess we have to assume
						// that the input jets are the first in the vector
						const std::vector<fastjet::PseudoJet>& internalJets = recluster.jets();

						// Find the minCore (WARNING, only hardest muon use for core finding) (CHECKED 3/26/15)
						{
							// Unfortunately, we can't just add the muon to the core and calculate
							// the mass, because the core currently has the wrong mass (from the granularity of
							// the cal, not from the actual particles). And we don't really want to recompute
							// the 4-momentum of every subjet with a new mass.

							// However, from a rather trivial calculation, one can show that, if the mass of
							// the subjet is constrained to CoreMassHypothesis, the M**2 after adding the muon is:
							//
							// M**2 = fCoreMassHypothesis2 + 2*Emu*Esub*(g + y)/(1 + (y + sqrt(1 - ((g-y) + g*y))))
							//
							// where (g = (fCoreMassHypothesis/Esub)**2) and (y = tan(theta)**2 = |(p3>_mu x p3>_sub)|**2 / (p3>_mu.p3>_sub)**2
							//
							// Luckily, we can calculate y without applying the CoreMassHypothesis to each subjet,
							// since their direction won't be altered by a new mass.

							const fastjet::PseudoJet& hardestMuon = internalJets[goodMuons.size() - 1];
							if(hardestMuon.user_index() >= 0)
								cout << "FAIL" << endl;

							// By minimizing the squared mass, we'll minimize the mass
							Double_t minMass2 = 9e9;
							unsigned int minCoreIndex = 0;

							for(unsigned int iSub = 0; iSub < subJets.size(); ++iSub)
							{
								// Make a copy of the subjet, so we can operate on its 4-vector
								fastjet::PseudoJet subjet = subJets[iSub];

								if(hardestMuon.is_inside(subjet))
								{
									// If the muon is already inside the subjet, subtract it first
									subjet -= hardestMuon;
								}

								const Double_t g = Squared(fCoreMassHypothesis / subjet.E());
								const Double_t y = Tan2(subjet, hardestMuon);

								// Technically, every squared mass has fCoreMassHypothesis2 (and is multiplied by 2),
								// so there is no point in doing these to everything
								// const Double_t mass2 = fCoreMassHypothesis2 + 2. * hardestMuon.E() * subjet.E() * (g + y) / (1. + (y + sqrt(1. - ((g - y) + g*y))));

								//WARNING: Currently this implementation does not account for cores with two muons of almost equal energy (small percentage of cores)
								const Double_t mass2 = hardestMuon.E() * subjet.E() * (g + y) / (1. + (y + sqrt(1. - ((g - y) + g*y))));

								if(mass2 < minMass2)
								{
									minMass2 = mass2;
									minCoreIndex = iSub;
								}
							}

							if(minCoreIndex > 0)
								coreVec.push_back(subJets[minCoreIndex]);

							// Bin the mass of the MinCore (which could be the HardCore)
							FillTH1F(mass, 1, coreFlavor, coreVec.back().m());
						}
						// If the minCore is the hardcore, coreVec will only have 1 member

						// Subtract muons from the cores (CHECKED 3/26/15)
						{
							for(unsigned int iMu = 0; iMu < goodMuons.size(); ++iMu)
							{
								const fastjet::PseudoJet& muon = internalJets[iMu];
								if(muon.user_index() >= 0)
									cout << "FAIL" << endl;

								for(std::vector<fastjet::PseudoJet>::iterator itCore = coreVec.begin(); itCore not_eq coreVec.end(); ++itCore)
								{
									if(muon.is_inside(*itCore))
									{
										*itCore -= muon;
										break; // The muon can only be in one core at a time
									}
								}
							}
						}

						// From now on, we will keep the original PseudoJet cores around in
						// case we are interested in constituents. However, their 4-vectors
						// will now be handled by TLorentzVector, since it has more/better math.

						// Fix core mass to fCoreMassHypothesis (CHECKED 3/26/15)
						{
							for(std::vector<fastjet::PseudoJet>::iterator itCore = coreVec.begin(); itCore not_eq coreVec.end(); ++itCore)
							{
								// The current mass depends more on the granularity of the CAL than anything else
								// To fix the mass, We need to scale the momentum of the core
								const Double_t momentumScale = sqrt((Squared(itCore->E()) - fCoreMassHypothesis2)/itCore->modp2());

								coreP4Vec.emplace_back(momentumScale*itCore->px(), momentumScale*itCore->py(), momentumScale*itCore->pz(), itCore->E()); //xyzt
							}
						}

						bool minCoreIsHardCore = (coreP4Vec.size() == 1);

						// Iterate through the cores and bin core/muon properties
						for(unsigned int iCore = 0; iCore < coreP4Vec.size(); ++iCore)
						{
							TLorentzVector& core = coreP4Vec[iCore];

							// Iterate through each muon, adding it back to the core, and calculate the resulting boost invariant
							for(unsigned int iMuon = 0; iMuon < goodMuons.size(); ++iMuon)
							{
								const TLorentzVector& muonP4 = goodMuons[iMuon]->Momentum;
								const TVector3 muonP3 = muonP4.Vect();
								const Double_t originalMass = core.M();
								const Int_t motherFlavor = goodMuonsMotherFlavor[iMuon];

								// Add the muon back to the core
								core += muonP4;

								// Only the HardCore adds the "neutrino" to the original jet
								if(iCore == 0)
									jet->Momentum += muonP4;

								const TLorentzVector& matriarchP4 = goodMuonsMatriarch[iMuon]->Momentum;
								const TVector3 matriarchP3 = matriarchP4.Vect();

								// Bin the angle from the raw core to the matriarch
								{
									const Double_t coreAngle2True = AccurateAngle(matriarchP3, core.Vect());

									FillTH1F(deltaTrue, iCore, motherFlavor, coreAngle2True);
									if(minCoreIsHardCore)
										FillTH1F(deltaTrue, 1, motherFlavor, coreAngle2True);
								}

								// Add the "neutrino" to the core
								core += muonP4;

								// Find the deltaMass
								{
									const Double_t delta_Mass = core.M() - originalMass;

									FillTH1F(deltaMass, iCore, motherFlavor, delta_Mass);
									if(minCoreIsHardCore)
										FillTH1F(deltaMass, 1, motherFlavor, delta_Mass);
								}

								// Bin the angle from the core (with muon twice) to matriarch
								{
									const Double_t coreAngle2True2Mu = AccurateAngle(matriarchP3, core.Vect());

									FillTH1F(deltaTrue2Mu, iCore, motherFlavor, coreAngle2True2Mu);
									if(minCoreIsHardCore)
										FillTH1F(deltaTrue2Mu, 1, motherFlavor, coreAngle2True2Mu);
								}

								const Double_t xTrue = (matriarchP4.E() * (muonP3.Cross(matriarchP3)).Mag()) /
										(matriarchP4.M() * muonP3.Dot(matriarchP3));

								// Bin xTrue
								{
									if(iCore == 0) // Only HardCore bins xTrue (i.e. once per muon), no core dependence
									{
										x_True[0]->Fill(xTrue);
										x_True[flavorToHistoIndex[motherFlavor]]->Fill(xTrue);
									}
								}


								/* Old Neutrino Code
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

								const Double_t boostMass = std::min(std::max(core.M(), fMinFinalMass), fMaxFinalMass);
								const Double_t	xCore = (core.E() * (muonP3.Cross(core.Vect())).Mag()) / (muonP3.Dot(core.Vect()) * boostMass);

								coreMinX[iCore] = std::min(coreMinX[iCore], xCore);

								// Bin x
								{
									FillTH1F(x_Core, iCore, motherFlavor, xCore);
									if(minCoreIsHardCore)
										FillTH1F(x_Core, 1, motherFlavor, xCore);
								}

								// Bin x vs. xTrue
								{
									FillTH2F(x_Core__vs__x_True, iCore, motherFlavor, xCore, xTrue);
									if(minCoreIsHardCore)
										FillTH2F(x_Core__vs__x_True, 1, motherFlavor, xCore, xTrue);
								}
							}// End loop over muons
						}// End loop over cores
					}// End valid subjets
				}// End Reclustering

				const Double_t jetPt = (jet->Momentum).Pt();
				bool minCoreIsHardCore = (coreP4Vec.size() == 1);

				for(unsigned int iCore = 0; iCore < coreP4Vec.size(); ++iCore)
				{
					if(iCore == 0)
					{
						pt_MuJets[0]->Fill(jetPt);
						pt_MuJets[flavorToHistoIndex[coreFlavor]]->Fill(jetPt);
					}

					const Double_t core_Ratio = coreP4Vec[iCore].Pt() / jetPt;

					FillTH1F(coreRatio, iCore, coreFlavor, core_Ratio);
					FillTH2F(x_Core__vs__coreRatio, iCore, coreFlavor, coreMinX[iCore], core_Ratio);
					if(minCoreIsHardCore)
					{
						FillTH1F(coreRatio, 1, coreFlavor, core_Ratio);
						FillTH2F(x_Core__vs__coreRatio, 1, coreFlavor, coreMinX[iCore], core_Ratio);
					}

					if(core_Ratio >= fMinCorePtRatio)
					{
						if(coreMinX[iCore] <= fMaxEmissionInvariant)
						{
							// Only tag jets from the minCore
							if(minCoreIsHardCore or iCore == 1)
								jet->BTag |= (1 << fBitNumber);

							FillTH1F(pt_Tagged, iCore, coreFlavor, jetPt);
							if(minCoreIsHardCore)
								FillTH1F(pt_Tagged, 1, coreFlavor, jetPt);
						}//End tagged
					}//End core is hard enough
				}//End loop through cores
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
