#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include <stdlib.h>
#include <signal.h>
#include <stdio.h>

#include "TROOT.h"
//#include "TDirectory.h"
#include "TString.h"
#include "TMath.h"
//#include "TApplication.h"
#include "TLorentzVector.h"

//#include "TFile.h"
//#include "TClonesArray.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "TNamed.h"

#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "LHCOWriterDirect.h"
//#include "ExRootAnalysis/ExRootTreeBranch.h"
//#include "ExRootAnalysis/ExRootProgressBar.h"

//using namespace std;

/*

LHCO+ format

See the header file for a detailed explanation of LHCO+

*/
char const* const LHCOWriterDirect::BASE_NAME = "BaseName";

const char* LHCOWriterDirect::BANNER = "   #  typ     eta     phi       pt      jmas     ntrk     btag   had/em     dum1     dum2\n";
const char* LHCOWriterDirect::EVENT_FIRST_LINE_FORMAT = "%4d %13u %8d\n";
const char* LHCOWriterDirect::EVENT_PARTICLE_LINE_FORMAT = "%4d %4d %7.3f %7.3f %8.3f %8.3f %7.3f %8.3f %7.3f %8.3f %8.5f\n";

const Float_t LHCOWriterDirect::PHOTON_MIN_PT_DEFAULT = 10.0;
const Float_t LHCOWriterDirect::ELECTRON_MIN_PT_DEFAULT = 5.0;
const Float_t LHCOWriterDirect::MUON_MIN_PT_DEFAULT = 5.0;
const Float_t LHCOWriterDirect::TAU_MIN_PT_DEFAULT = 10.0;
const Float_t LHCOWriterDirect::JET_MIN_PT_DEFAULT = 10.0;

// WARNING: Both of these values depend on the precision in EVENT_PARTICLE_LINE_FORMAT
const Float_t LHCOWriterDirect::MAX_HADEM = 999.999;
const Float_t LHCOWriterDirect::MAX_TANH_DELTA_R = .999;

const Int_t LHCOWriterDirect::EFLOW_PHOTON = 0;
const Int_t LHCOWriterDirect::ELECTRON = 1;
const Int_t LHCOWriterDirect::MUON = 2;
const Int_t LHCOWriterDirect::EFLOW_TAU = 3;
const Int_t LHCOWriterDirect::EFLOW_JET = 4;
const Int_t LHCOWriterDirect::MISSING_ET = 6;
const Int_t LHCOWriterDirect::SCALAR_HT = 7;
const Int_t LHCOWriterDirect::TOWER_JET = 14;
const Int_t LHCOWriterDirect::PYTHIA_PARTON_JET = 34;
const Int_t LHCOWriterDirect::PYTHIA_SINGLET_JET = 35;

const Int_t LHCOWriterDirect::NON_ISOLATED_FLAG =  -99999;
const Int_t LHCOWriterDirect::ISO_PHOTON_FLAG   = -100000;
const Int_t LHCOWriterDirect::ISO_ELECTRON_FLAG = -100001;
const Int_t LHCOWriterDirect::ISO_MUON_FLAG     = -100002;

const size_t LHCOWriterDirect::BUFFER_SIZE = 8*1024*1024; // (8 MB) This is to reduce I/O conflicts between reading the LHE and outputing the LHCO

const Float_t TWO_PI = 2. * TMath::Pi();
const Float_t PI = TMath::Pi();

LHCOParticleLine::LHCOParticleLine(Int_t type_in, Float_t eta_in, Float_t phi_in, Float_t pt_in, Float_t jmas_in) :
type (type_in), eta(eta_in), phi(phi_in), pt(pt_in), jmas(jmas_in), ntrk(0.), btag(0.), hadEm(0.), dum1(0.), dum2(0.) {}

LHCOWriterDirect::LHCOWriterDirect() :
	fOutputFile(0),
	fAllParticles(0), fEFlowPhotons(0), fElectrons(0), fMuons(0),
	fEFlowJets(0), fTowerJets(0), fMissingEt(0), fScalarHt(0), fPythiaPartonJets(0), fPythiaSingletJets(0),
	fItEFlowPhotons(0), fItElectrons(0), fItMuons(0),
	fItEFlowJets(0), fItTowerJets(0), fItPythiaPartonJets(0), fItPythiaSingletJets(0),
	photonMinPt(0.), electronMinPt(0.), muonMinPt(0.), tauMinPt(0.), jetMinPt(0.),
	subtractMuons(false),
	eventRecord(), jetToConstituent(),
	eventCounter(0), lineCounter(0), jetBegin(0), jetEnd(0)
{}

LHCOWriterDirect::~LHCOWriterDirect()
{}

void LHCOWriterDirect::Init()
{
	// Open the output file
	{
		// The main which initiates the Delphes module must create a TNamed(LHCOWriterDirect::BASE_NAME, baseNameOfLHCOFile)
		// and add the TNamed object to gROOT. This allows the TNamed to be searchable by its "Name" (LHCOWriterDirect::BASE_NAME)
		// so that the TName's "Title" (baseNameOfLHCOFile) can be used by this module.
		TNamed const* const baseName = static_cast<TNamed*>(gROOT->Get(LHCOWriterDirect::BASE_NAME));

		if(not baseName)
		{
			std::stringstream message;
			message << "Can't locate TNamed named (\"" << BASE_NAME << "\") which must be set in the main() \n that initiates Delphes! This is how LHCOWriterDirect knows how to name its file.\n";
			throw std::runtime_error(message.str()); // An unhandled exception - oops!
		}

		std::string fileName(baseName->GetTitle());
		fileName = fileName.substr(0, fileName.find_first_of('.')) + ".lhco";

		fOutputFile = fopen(fileName.c_str(), "w"); // Open the file for replacement

		if(fOutputFile == NULL) // Check that the file is open
		{
			std::stringstream message;
			message << "LHCOWriterDirect: Can't open  <" << fileName << ">\n";
			throw std::runtime_error(message.str()); // An unhandled exception - oops!
		}
	}

	setvbuf(fOutputFile, 0, _IOFBF, BUFFER_SIZE); // Set full buffering, with suggested size (the null pointer forces automatic allocation)

	fprintf(fOutputFile, LHCOWriterDirect::BANNER); // Write the first line of the file
	eventCounter = 0; // There are no events yet

	// Check for supplied parameters, hook up the arrays and initialize the iterators (when appropriate)

	const char* arrayName = GetString("AllParticles", "Delphes/allParticles");
	if(arrayName[0] != '\0')
		fAllParticles = ImportArray(arrayName);

	arrayName = GetString("EFlowPhotons", "");
	if(arrayName[0] != '\0')
	{
		fEFlowPhotons = ImportArray(arrayName);
		fItEFlowPhotons = fEFlowPhotons->MakeIterator();
	}

	arrayName = GetString("Electrons", "");
	if(arrayName[0] != '\0')
	{
		fElectrons = ImportArray(arrayName);
		fItElectrons = fElectrons->MakeIterator();
	}

	arrayName = GetString("Muons", "");
	if(arrayName[0] != '\0')
	{
		fMuons = ImportArray(arrayName);
		fItMuons = fMuons->MakeIterator();
	}

	arrayName = GetString("EFlowJets", "");
	if(arrayName[0] != '\0')
	{
		fEFlowJets = ImportArray(arrayName);
		fItEFlowJets = fEFlowJets->MakeIterator();
	}

	arrayName = GetString("MissingEt", "");
	if(arrayName[0] != '\0')
		fMissingEt = ImportArray(arrayName);

	arrayName = GetString("ScalarHt", "");
	if(arrayName[0] != '\0')
		fScalarHt = ImportArray(arrayName);

	arrayName = GetString("TowerJets", "");
	if(arrayName[0] != '\0')
	{
		fTowerJets = ImportArray(arrayName);
		fItTowerJets = fTowerJets->MakeIterator();
	}

	arrayName = GetString("PythiaPartonJets", "");
	if(arrayName[0] != '\0')
	{
		fPythiaPartonJets = ImportArray(arrayName);
		fItPythiaPartonJets = fPythiaPartonJets->MakeIterator();
	}

	arrayName = GetString("PythiaSingletJets", "");
	if(arrayName[0] != '\0')
	{
		fPythiaSingletJets = ImportArray(arrayName);
		fItPythiaSingletJets = fPythiaSingletJets->MakeIterator();
	}

	photonMinPt = GetDouble("PhotonMinPt", LHCOWriterDirect::PHOTON_MIN_PT_DEFAULT);
	electronMinPt = GetDouble("ElectronMinPt", LHCOWriterDirect::ELECTRON_MIN_PT_DEFAULT);
	muonMinPt = GetDouble("MuonMinPt", LHCOWriterDirect::MUON_MIN_PT_DEFAULT);
	tauMinPt = GetDouble("TauMinPt", LHCOWriterDirect::TAU_MIN_PT_DEFAULT);
	jetMinPt = GetDouble("JetMinPt", LHCOWriterDirect::JET_MIN_PT_DEFAULT);

	subtractMuons = GetBool("SubtractMuons", false);
}

void LHCOWriterDirect::Process()
{
	eventRecord.clear(); // Clear the record
	jetToConstituent.clear(); // Clear the mapper
	lineCounter = 0; // reset the counter

	// Process the potentially isolated candidates
	// The @nearJet and @jetIn will be filled once we start processing jets
	ProcessEFlowPhotons();
	ProcessElectrons();
	ProcessMuons();

	// Process the jets (and taus)
	ProcessEFlowJets();

	// Process the energy statistics
	ProcessMissingEt();
	//ProcessScalarHt(); // Not supported yet

	ProcessTowerJets();
	ProcessPythiaPartonJets();
	ProcessPythiaSingletJets();

	WriteLines(); // Write out the event
}

void LHCOWriterDirect::Finish()
{
	fclose(fOutputFile);
	// Delete iterators
	// No need to check if they're non-zezo, delete will not delete the nullpointer
	delete fItEFlowPhotons;
	delete fItElectrons;
	delete fItMuons;
	delete fItEFlowJets;
	delete fItTowerJets;
	delete fItPythiaPartonJets;
	delete fItPythiaSingletJets;
}

void LHCOWriterDirect::WriteLines()
{
	if(eventRecord.empty())
		return;
	else
	{
		fprintf(fOutputFile, LHCOWriterDirect::EVENT_FIRST_LINE_FORMAT, 0, eventCounter++, 0); // Write the event line. Currently there is no trigger word (the last argument).

		Int_t lineNumber = 0;

		// Iterate through the lines and write them out
		for(std::vector<LHCOParticleLine>::iterator itLine = eventRecord.begin(); itLine != eventRecord.end(); ++itLine)
		{
			LHCOParticleLine& line = *itLine;
			fprintf(fOutputFile, LHCOWriterDirect::EVENT_PARTICLE_LINE_FORMAT,
			++lineNumber, line.type, line.eta, line.phi, line.pt, line.jmas, line.ntrk, line.btag, line.hadEm, line.dum1, line.dum2);
		}

		// Explicitly prevent double writing
		eventRecord.clear(); // Clear the record
		lineCounter = 0; // reset the counter
	}
}

// Inserts a particles 4-momentum into the LHCO event record
// If lineNumberBTag, puts the line number of the event in the BTag field
bool LHCOWriterDirect::FillKinematicsAndPushBack(Candidate* candidate, Int_t type, Float_t ptMin)
{
	TLorentzVector& momentum = candidate->Momentum; // Get the candidates 4-momentum
	Float_t pt = momentum.Pt();

	if(pt >= ptMin)
	{
		// Here we would prefer to use emplace_back, but alas, cannot trust that c++11 will be fully supported
		eventRecord.push_back(LHCOParticleLine(type, momentum.Eta(), momentum.Phi(), pt, momentum.M()));
		// momentum.M() is better than Candidate->Mass, as it informs us which 4-momentum was used inside Delphes itself
		// If Delphes is using massless 4-momentum when it adds 4-vectors, we should not pretend it isn't

		++lineCounter; // We have a new line in the event record

		return true;
	}
	else return false;
}

// Used to find the original mother of a lepton
// Will trace back until the mother is not another lepton
Float_t LHCOWriterDirect::MotherFlavor(Candidate* candidate)
{
	Int_t flavor = -99; // An unsupported flavor, indicating that fAllParticles was not supplied

	if(fAllParticles and (candidate->IsPU == 0)) // Cannot assume that pileup is dereferenceable
	{
		Int_t pid = 0; // Another unsupported flavor, indicating that something went wrong

		while((candidate = static_cast<Candidate*>(fAllParticles->At(candidate->M1))))  // Try to get the mother1
		{
			pid = TMath::Abs(candidate->PID);

			if(!(pid == 13 || pid == 15)) // If the particle is NOT a muon or a tau (a lepton from a lepton), then we've found the mother
				break;
		}

		flavor = pid;

		/* temporary print out entire mother ID
		if(pid < 100) // The mother was not a hadron, so it has no flavor: output mom's absPDGID as a negative number
			flavor = -pid;
		else if(pid >= 1000 && pid < 10000) // Baryon, heaviest quark in 1000s place
			flavor = pid / 1000;
		else // Meson, heaviest quark in 100s place
			flavor = (pid%1000)/100;
		*/
	}

	return Float_t(flavor);
}

void LHCOWriterDirect::SubtractIsolated_CalculateJetStats(Candidate* jet)
{
	// Zero the internal counters
	Float_t EhadSum = 0., EemSum = 0.;
	Int_t tracks = 0, chargeSum = 0;

	{
		TIterator* fItConstituents = jet->GetCandidates()->MakeIterator(); // Make an iterator for the jet's constituents

		Candidate* constituent;

		// Loop through the jet constituents
		while((constituent = static_cast<Candidate*>(fItConstituents->Next())))
		{
			Int_t& constituentTag = constituent->IsConstituent;
			if(constituentTag < LHCOWriterDirect::NON_ISOLATED_FLAG && constituentTag >= ISO_MUON_FLAG)
			{
				// A jet constituent should not have an IsConstituent flag. The only way this should be non-zero
				// is if we've already placed this particle in the event record and set the BTag ourselves

				// Associate the constituent pointer with the jet
				jetToConstituent[jet].push_back(constituent);

				// Subtract the particle from the jet and don't allow it to contribute to jet statistics
				// That is, unless it's a muon and we're not subtracting muons
				if(!(constituentTag == ISO_MUON_FLAG && !subtractMuons))
				{
					jet->Momentum -= constituent->Momentum;
					constituentTag = -1; // Record that we've subtracted the constituent
					continue; // Don't sum the jet statistics
				}
				else
					constituentTag = 1; // Record that we haven't subtracted the constituent
			}

			// Sum the EM and HAD energy of all constituents
			// This will work perfectly assuming that tracks have their Ehad and Eem stored inside them
			// (which is not the Delphes default but has been manuall added by KDP in Calorimeter)
			// Otherwise, tracks will not contribute to this ratio.
			// Muons will not be summed, but EHAD is supposed to reflect the calorimeter state
			EhadSum += constituent->Ehad;
			EemSum += constituent->Eem;

			if(constituent->Charge != 0) // If it's a track
			{
				tracks++;
				chargeSum += constituent->Charge;
			}
		}

		delete fItConstituents;
	}

	// These fields are unused for jets
	// Actually, Charge is used for Taus, but it is determined clairvoyently from the tau, so it isn't very meaningful
	jet->Status = tracks;
	jet->Charge = chargeSum;
	jet->Ehad = EhadSum;
	jet->Eem = EemSum;
}

void LHCOWriterDirect::AssociateConstituents(Candidate* jet, UInt_t jetLine)
{
	std::map<Candidate*, std::vector<Candidate*> >::iterator itMap(jetToConstituent.find(jet));

	// Make sure this jet has constituents
	if(itMap != jetToConstituent.end())
	{
		std::vector<Candidate*>& constituents = itMap->second;

		// iterate through the constituent pointers
		for(std::vector<Candidate*>::iterator itConstituents = constituents.begin(); itConstituents != constituents.end(); ++itConstituents)
			eventRecord[(*itConstituents)->BTag - 1].dum1 = ((*itConstituents)->IsConstituent)*Int_t(jetLine); // the lineNum of the constituent is stored in BTag
	}
}

Float_t LHCOWriterDirect::DeltaRSquared(std::vector<LHCOParticleLine>::const_iterator one, std::vector<LHCOParticleLine>::const_iterator two)
{
	// Retutn the DeltaR of the two particles in the event record (used in finding the closest jet)
	// If we needed to calculate phi directly, it would be faster to have a optional third argument which gives the current min
	// That way, if the new deltaR is already larger from eta, don't proceed.
	// But, since we already have access to eta and phi already, just calculate it (the number of conditional jumps are less with this way)

	Float_t absDPhi = TMath::Abs(one->phi - two->phi);
	if(absDPhi > PI)// Check that we're using the interior angle, not the exterior one
		absDPhi -= TWO_PI;

	Float_t dEta = (one->eta - two->eta);

	return absDPhi*absDPhi + dEta*dEta;
}

void LHCOWriterDirect::FindClosestJets(
	std::vector<LHCOParticleLine>::iterator beginCandidate, // inclusive
	std::vector<LHCOParticleLine>::iterator endCandidate, //exclusive
	std::vector<LHCOParticleLine>::iterator beginJet, //inclusive
	std::vector<LHCOParticleLine>::iterator endJet) // exclusive
{
	// For photons, electrons, muons and taus, find the line# of the closest jet and store it in BTag field
	// WARNING: For now, this function assumes that it will only be called after jets but before anything else

	// Loop through everything before jets
	for(std::vector<LHCOParticleLine>::iterator itCan = beginCandidate; itCan != endCandidate; ++itCan)
	{
		// We'll keep track of the smallest delDeltaRSquaredtaR2 and the corresponding iterator
		Float_t deltaR2Min = 9e99; // No jet (that we can detect) can be farther away than this
		std::vector<LHCOParticleLine>::iterator closestJet = eventRecord.end();

		// Loop through the jets (assuming nothing is after jets)
		for(std::vector<LHCOParticleLine>::iterator itJet = beginJet; itJet != endJet; ++itJet)
		{
			Float_t deltaR2 = DeltaRSquared(itCan, itJet);

			// If the new jet is closer, hooray!
			if(deltaR2 < deltaR2Min)
			{
				deltaR2Min = deltaR2;
				closestJet = itJet;
			}
		}

		if(closestJet != eventRecord.end())
			itCan->dum2 = Float_t(closestJet - eventRecord.begin() + 1) + TMath::Min(Float_t(TMath::TanH(TMath::Sqrt(deltaR2Min))), MAX_TANH_DELTA_R);
		// Store the closest jet AND the distance
		// int(line# of closest jet) + min(tanh(deltaR), .999)
	}
}

void LHCOWriterDirect::ProcessEFlowPhotons()
{
	if(fEFlowPhotons)
	{
		Candidate* photon;

		fEFlowPhotons->Sort(); // These are Candidates, so they will be sorted by Pt
		fItEFlowPhotons->Reset();

		while((photon = static_cast<Candidate*>(fItEFlowPhotons->Next())))
		{
			if(FillKinematicsAndPushBack(photon, LHCOWriterDirect::EFLOW_PHOTON, photonMinPt))
			{
				photon->IsConstituent = LHCOWriterDirect::ISO_PHOTON_FLAG; //Record that this is a photon
				photon->BTag = lineCounter; //Record the photon's line

				LHCOParticleLine& photonLine = eventRecord.back();

				photonLine.ntrk = photon->Tau[0]; // tracks in Cal cell
				photonLine.btag = TMath::Min(Float_t(photon->Momentum.E()) / photon->Tau[1], LHCOWriterDirect::MAX_HADEM); // EFlow_E / E_subtracted
				photonLine.hadEm = photon->Tau[2]; // ptiso
			}
			else
			{
				photon->IsConstituent = LHCOWriterDirect::NON_ISOLATED_FLAG;
				photon->BTag = 0;
			}
		}
	}
}

void LHCOWriterDirect::ProcessElectrons()
{
	if(fElectrons)
	{
		Candidate* electron;

		fElectrons->Sort(); // These are Candidates, so they will be sorted by Pt
		fItElectrons->Reset();

		while((electron = static_cast<Candidate*>(fItElectrons->Next())))
		{
			if(FillKinematicsAndPushBack(electron, LHCOWriterDirect::ELECTRON, electronMinPt))
			{
				electron->IsConstituent = LHCOWriterDirect::ISO_ELECTRON_FLAG;
				electron->BTag = lineCounter;

				LHCOParticleLine& electronLine = eventRecord.back();

				electronLine.ntrk = electron->Charge;
				electronLine.btag = MotherFlavor(electron);
				electronLine.hadEm = electron->Tau[2]; // ptiso
			}
			else
			{
				electron->IsConstituent = LHCOWriterDirect::NON_ISOLATED_FLAG;
				electron->BTag = 0;
			}
		}
	}
}

void LHCOWriterDirect::ProcessMuons()
{
	if(fMuons)
	{
		Candidate* muon;

		fMuons->Sort(); // These are Candidates, so they will be sorted by Pt
		fItMuons->Reset();

		while((muon = static_cast<Candidate*>(fItMuons->Next())))
		{
			if(FillKinematicsAndPushBack(muon, LHCOWriterDirect::MUON, muonMinPt))
			{
				muon->IsConstituent = LHCOWriterDirect::ISO_MUON_FLAG;
				muon->BTag = lineCounter;

				LHCOParticleLine& muonLine = eventRecord.back();

				muonLine.ntrk = muon->Charge;
				muonLine.btag = MotherFlavor(muon);
				muonLine.hadEm = muon->Tau[2]; // ptiso
			}
			else
			{
				muon->IsConstituent = LHCOWriterDirect::NON_ISOLATED_FLAG;
				muon->BTag = 0;
			}
		}
	}
}

void LHCOWriterDirect::ProcessEFlowJets()
{
	if(fEFlowJets)
	{
		Candidate* candidate;

		fItEFlowJets->Reset();

		jetToConstituent.clear();
		// Find isolated leptons and photons and subtract them (only subtract muons if subtractMuons == true)
		while((candidate = static_cast<Candidate*>(fItEFlowJets->Next())))
			SubtractIsolated_CalculateJetStats(candidate);

		// Now sort the jets
		fEFlowJets->Sort();
		fItEFlowJets->Reset();

		std::vector<Candidate*> jetVec, tauVec;

		// We need to seperate the tauJets from the regular jets, and proces them sequentially
		while((candidate = static_cast<Candidate*>(fItEFlowJets->Next())))
		{
			if(candidate->TauTag > 0)
				tauVec.push_back(candidate);
			else
				jetVec.push_back(candidate);
		}

		// Process the taus
		for(std::vector<Candidate*>::iterator itCan = tauVec.begin(); itCan != tauVec.end(); ++itCan)
		{
			candidate = *itCan;

			if(FillKinematicsAndPushBack(candidate, LHCOWriterDirect::EFLOW_TAU, tauMinPt))
			{
				LHCOParticleLine& tauLine = eventRecord.back();

				tauLine.ntrk = candidate->Status*candidate->Charge; // tracks * charge
				// Currently, delphes TauTag is either 0 or 1, so there is no point in recording the TauTag
				tauLine.hadEm = TMath::Min(candidate->Ehad / candidate->Eem, LHCOWriterDirect::MAX_HADEM);

				tauLine.dum1 = (candidate->Area).Pt(); // The scalar area from the jet algorithm (will be 0 if run without area finding)

				// The detected area is approx the same as the widest spread particles: sqrt(DeltaEta^2 + DeltaPhi^2)
				Float_t tmp = candidate->DeltaEta;
				tmp*=tmp;
				tauLine.dum2 = tmp;
				tmp = candidate->DeltaPhi;
				tmp*=tmp;
				tauLine.dum2 += tmp;
				tauLine.dum2 *= PI;

				AssociateConstituents(candidate, lineCounter);
			}
			else break;// These are sorted by Pt, once the first one fails the pt cut, all subsequent ones will too
		}

		if(jetVec.empty())
			return;

		jetBegin = eventRecord.size(); // It would be preferrable to store an iterator, but eventRecord may still grow

		// Fill the jets
		for(std::vector<Candidate*>::iterator itCan = jetVec.begin(); itCan != jetVec.end(); ++itCan)
		{
			candidate = *itCan;

			if(FillKinematicsAndPushBack(candidate, LHCOWriterDirect::EFLOW_JET, jetMinPt))
			{
				LHCOParticleLine& jetLine = eventRecord.back();

				Float_t dTracks(candidate->Status);

				// Here we have the number of tracks AND the charge embedded in the same number
				if(candidate->Charge != 0) // but only if chargeSum != 0, otherwise will have a divide by 0
					jetLine.ntrk = candidate->Charge * (dTracks / TMath::Abs(candidate->Charge) + 1 / dTracks); // WARNING: there remains an ambiguity if all tracks have the same charge (vs. if they perfeclty balance in charge)
				else
					jetLine.ntrk = dTracks;

				jetLine.hadEm = TMath::Min(candidate->Ehad / candidate->Eem, LHCOWriterDirect::MAX_HADEM);

				jetLine.btag = candidate->BTag;
				jetLine.dum1 = (candidate->Area).Pt(); // The scalar area from the jet algorithm (will be 0 if run without area finding)

				// The detected area is approx the same as the widest spread particles: sqrt(DeltaEta^2 + DeltaPhi^2)
				Float_t tmp = candidate->DeltaEta;
				tmp*=tmp;
				jetLine.dum2 = tmp;
				tmp = candidate->DeltaPhi;
				tmp*=tmp;
				jetLine.dum2 += tmp;
				jetLine.dum2 *= PI;

				AssociateConstituents(candidate, lineCounter);
			}
			else break;// These are sorted by Pt, once the first one fails the pt cut, all subsequent ones will too
		}

		jetEnd = eventRecord.size();

		jetToConstituent.clear();

		// Now that all the jets are in place, find the closest jet for each of the non-jets
		// All photons, electrons, muons, and taus have the nearest jet's line # in their btag field
		FindClosestJets(eventRecord.begin(), eventRecord.begin() + jetBegin, eventRecord.begin() + jetBegin, eventRecord.end());
	}
}

void LHCOWriterDirect::ProcessMissingEt()
{
	if(fMissingEt)
	{
		Candidate* Et = static_cast<Candidate*>(fMissingEt->At(0)); // There can be only one (i.e. no need to iterate)

		if(Et)
		{
			FillKinematicsAndPushBack(Et, LHCOWriterDirect::MISSING_ET, 0.);
			LHCOParticleLine& EtLine = eventRecord.back();

			// Delphes calculates the sum of all 4-momenta for its MissingEt, we must reverse eta and phi
			EtLine.eta = -EtLine.eta;

			// We're using phi in [-Pi, Pi], so reversing phi is easy
			if(EtLine.phi > 0.)
				EtLine.phi -= PI;
			else if(EtLine.phi < 0.)
				EtLine.phi += PI;
		}
	}
}

void LHCOWriterDirect::ProcessTowerJets()
{
	// Really need to have tower photons as well

	if(fTowerJets)
	{
		Candidate* towerJet;

		fTowerJets->Sort(); // These are Candidates, so they will be sorted by Pt
		fItTowerJets->Reset();

		int towerBegin = eventRecord.size();

		while((towerJet = static_cast<Candidate*>(fItTowerJets->Next())))
		{
			if(FillKinematicsAndPushBack(towerJet, LHCOWriterDirect::TOWER_JET, jetMinPt))
			{
				LHCOParticleLine& towerLine = eventRecord.back();

				Float_t EhadSum = 0., EemSum = 0.;
				// Zero the counters passed by reference
				Int_t tracks = 0;

				{
					TIterator* fItConstituents = towerJet->GetCandidates()->MakeIterator(); // Make an iterator for the jet's constituents

					Candidate* constituent;

					// Loop through the jet constituents
					while((constituent = static_cast<Candidate*>(fItConstituents->Next())))
					{
						// Sum the EM and HAD energy of all constituents
						EhadSum += constituent->Ehad;
						EemSum += constituent->Eem;

						tracks += constituent->Tau[0];
					}

					delete fItConstituents;
				}

				towerLine.ntrk = tracks;
				towerLine.btag = towerJet->BTag;
				towerLine.hadEm = TMath::Min(EhadSum/EemSum, LHCOWriterDirect::MAX_HADEM); // Return Had/Em, limiting the maximum ratio
			}
			else break;// These are sorted by Pt, once the first one fails the pt cut, all subsequent ones will too
		}

		int towerEnd = eventRecord.size();

		FindClosestJets(eventRecord.begin() + towerBegin, eventRecord.begin() + towerEnd,
			eventRecord.begin() + jetBegin, eventRecord.begin() + jetEnd);
	}
}

void LHCOWriterDirect::ProcessPythiaPartonJets()
{
	if(fPythiaPartonJets)
	{
		Candidate* partonJet;

		fPythiaPartonJets->Sort(); // These are Candidates, so they will be sorted by Pt
		fItPythiaPartonJets->Reset();

		int partonBegin = eventRecord.size();

		while((partonJet = static_cast<Candidate*>(fItPythiaPartonJets->Next())))
		{
			if(FillKinematicsAndPushBack(partonJet, LHCOWriterDirect::PYTHIA_PARTON_JET, jetMinPt))
			{
				LHCOParticleLine& partonLine = eventRecord.back();

				Int_t quarkCount = 0, heaviestQuarkPID = 0, hardestPartonPID = 21;
				Float_t maxPt = 0.;
				TIterator* ItConstituent = partonJet->GetCandidates()->MakeIterator();
				ItConstituent->Reset();

				Candidate* subJet;

				while((subJet = static_cast<Candidate*>(ItConstituent->Next())))
				{
					Int_t absPID = TMath::Abs(subJet->PID);

					if(absPID != 21)
					{
						++quarkCount;

						if(absPID > 6) // diquark
						{
							++quarkCount;
							absPID /= 1000; // convert to heaviest quark in diquark
						}

						if(absPID > heaviestQuarkPID)
							heaviestQuarkPID = absPID;
					}

					Float_t newPt = subJet->Momentum.Pt();

					if(newPt > maxPt)
					{
						maxPt = newPt;
						hardestPartonPID = absPID;
					}
				}

				partonLine.ntrk = quarkCount;
				partonLine.btag = heaviestQuarkPID;
				partonLine.dum1 = hardestPartonPID;

				delete ItConstituent;
			}
			else break;// These are sorted by Pt, once the first one fails the pt cut, all subsequent ones will too
		}

		int partonEnd = eventRecord.size();

		FindClosestJets(eventRecord.begin() + partonBegin, eventRecord.begin() + partonEnd,
			eventRecord.begin() + jetBegin, eventRecord.begin() + jetEnd);
	}
}

void LHCOWriterDirect::ProcessPythiaSingletJets()
{
	if(fPythiaSingletJets)
	{
		Candidate* singletJet;

		fPythiaSingletJets->Sort(); // These are Candidates, so they will be sorted by Pt
		fItPythiaSingletJets->Reset();

		int singletBegin = eventRecord.size();

		while((singletJet = static_cast<Candidate*>(fItPythiaSingletJets->Next())))
		{
			if(!FillKinematicsAndPushBack(singletJet, LHCOWriterDirect::PYTHIA_SINGLET_JET, false))
				throw std::runtime_error("wtf");

			LHCOParticleLine& singletLine = eventRecord.back();

			Int_t contentCount = 0, heaviestQuarkPID = 1, hardestPartonPID = 21;
			Float_t maxPt = 0.;
			TIterator* ItConstituent = singletJet->GetCandidates()->MakeIterator();
			ItConstituent->Reset();

			Candidate* subJet;

			while((subJet = static_cast<Candidate*>(ItConstituent->Next())))
			{
				++contentCount;
				{
					Float_t newPt = subJet->Position.Pt();

					if(newPt > maxPt)
					{
						maxPt = newPt;
						hardestPartonPID = subJet->PID;
					}
				}
				{
					Int_t quarkPID = subJet->Charge;

					if(quarkPID > 1000)
						quarkPID /= 1000;

					if(quarkPID > heaviestQuarkPID)
						heaviestQuarkPID = quarkPID;
				}
			}

			singletLine.ntrk = contentCount;
			singletLine.btag = heaviestQuarkPID;
			singletLine.dum1 = hardestPartonPID;

			delete ItConstituent;
		}

		int singletEnd = eventRecord.size();

		FindClosestJets(eventRecord.begin() + singletBegin, eventRecord.begin() + singletEnd,
			eventRecord.begin() + jetBegin, eventRecord.begin() + jetEnd);
	}
}
