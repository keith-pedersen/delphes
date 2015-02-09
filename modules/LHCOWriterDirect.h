#ifndef LHCOWriterDirect_h
#define LHCOWriterDirect_h

/** \class LHCOWriterDirect
 *
 *  Writes physics objects to a human-readable LHCO file. 
 *
 *  \author K. Pedersen - Illinois Institute of Technology (https://github.com/keith-pedersen/)
 *
 */

#include "classes/DelphesModule.h"
#include <map>

// This class assumes that it is the last class to access the Candidates
// It may create side affects, as it alters the following fields:
//      Candidate::BTag (this is only altered for photon, electrons, and muons, so there should be no real side effects)

/*

LHCO+ format

All energy in GeV
===========
 0-9  | Detector observables with EFlow
10-19 | Detector observables without EFlow
20-29 | Detector observables in a perfect detector (that can see neutirinos (and thus taus))
20-29 | Unobservables (colored particles)

type# | type string        | eta | phi | pt | jmas  |  ntrk            |  btag     | had/em | dum1     | dum2     |
======|====================|=====|=====|====|=======|==================|===========|========|==========|==========|
    0 | EFlowPhoton        | eta | phi | pt | @mass |  @trkTower       | @E/Esub   | @ptiso | @jetIn   | @nearJet |
    1 | Electron           | eta | phi | pt | @mass |   charge         | @maFlavor | @ptiso | @jetIn   | @nearJet |
    2 | Muon               | eta | phi | pt | @mass |   charge         | @maFlavor | @ptiso | @jetIn   | @nearJet |
    3 | EFlowTau           | eta | phi | pt | @mass | charge * prongs  |!!@maFlavor| had/em | @algArea | @nearJet |
    4 | EFlowJet           | eta | phi | pt | @mass |     @trackQ      | BTag (SV) | had/em | @algArea | @detArea |
    6 | MissingEt          | eta | phi | pt | @mass | !!total tracks   |     ?     |    ?   |    ?     |     ?    |
  !!7 | ScalarHT           |  0  |  0  | pt |   0   |        ?         |     ?     |    ?   |    ?     |     ?    |
  !!8 | Charged Hadron     | eta | phi | pt | @mass |   charge         |     ?     |    ?   |    ?     | @nearJet |
  !!9 | EFlowNeutralHadron | eta | phi | pt | @mass |        ?         |     ?     |    ?   |    ?     |     ?    |
======|====================|=====|=====|====|=======|==================|===========|========|==========|==========|
 !!10 | TowerPhoton        | eta | phi | pt | @mass |        ?         |     ?     | had/em |    ?     | @nearJet |
 !!13 | TowerTau           | eta | phi | pt | @mass | charge * prongs  | @maFlavor | had/em |  TauTag  | @nearJet |
 !!14 | TowerJet           | eta | phi | pt | @mass |     @trkTower    |!!BTag (SV)| had/em |    ?     | @nearJet |
======|====================|=====|=====|====|=======|==================|===========|========|==========|==========|
 !!24 | GenJet             | eta | phi | pt | @mass |     @trackQ      | @myFlavor | had/em |    ?     |     ?    |
 !!25 | Neutrino           | eta | phi | pt | @mass |        ?         | @myBuddy  |    ?   |    ?     |     ?    |
======|====================|=====|=====|====|=======|==================|===========|========|==========|==========|
 !!34 | PythiaPartonJet    | eta | phi | pt | @mass |  quark count     | @quarkID  |    ?   | @partonID| @nearJet |
 !!35 | PythiaSingletJet   | eta | phi | pt | @mass |  parton count    | @quarkID  |    ?   | @partonID| @nearJet |
 !!36?| Parton             | eta | phi | pt | @mass |        ?         |  PDGID    |    ?   |    ?     |     ?    |

@mass     =  Not Candidate.Mass but Candidate.Momentum.M(). This is a safety check, since 4-momenta are being added inside Delphes.
@trkTower =  The total number of tracks that flowed through the tower of the EFlowPhoton. Stored in Candidate.Tau[0]
@nearJet  =  int(line# of nearest jet) + min(tanh(deltaR to jet), .999)
@ptiso    =  ptSum/pt     [ptSum = sum of pt in circle of deltaR = 0.5]  This is stored in Candidate.Tau[2]
@jetIn    =  Line# of jet the particle was absorbed into
@E/Esub   =  max(EFlowPhoton energy / (energy subtracted from tracks), 999.99) energy subtracted is stored in Candidate.Tau[1]
@maFlavor =  If positive: the flavor of the heaviest quark from the mother hadron. If negative: the generation of the mother lepton. If -99, broken. 
@trackQ   =  sumQ*(trks/abs(sumQ) + 1/trks) (i.e. the number of tracks as well as how well they balance charge)
@detArea  =  sqrt(dEta^2 + dPhi^2) (where dEta is the distance from the centroid to the widest constituent in eta and dPhi is the same for phi)
@algArea  =  The pt of the area^mu 4-vector returned by FastJet 
@myFlavor =  The flavor of my heaviest constituent (not neccessarily my most energetic). WARNING: Might be gluon splitting
@myBuddy  =  The line# of my lepton sibling
@quarkID  =  The absPID of the highest pt quark
@partonID =  The absPID of the highest pt parton


    ?     =  0 (i.e. it might hold something in the future besides 0, but for now it is 0)
    $     =  may not be possible to determine in current incarnation of Delphes
   !!     =  not yet implemented
    !     =  only partially implemented

*/

class TIterator;
class TObjArray;
class DelphesFormula;
class string;
class Candidate;

struct LHCOParticleLine
{
	Int_t type; 
	Float_t eta, phi, pt, jmas, ntrk, btag, hadEm, dum1, dum2;
	
	LHCOParticleLine(Int_t type_in = -1, Float_t eta_in = 0., Float_t phi_in = 0., Float_t pt_in = 0., Float_t jmas_in = 0.);
};

class LHCOWriterDirect: public DelphesModule
{
	public:
		LHCOWriterDirect();
		~LHCOWriterDirect();

		void Init();
		void Process();
		void Finish(); 
		
		static const char* BANNER;
		static const char* EVENT_FIRST_LINE_FORMAT;
		static const char* EVENT_PARTICLE_LINE_FORMAT;
		
		static const Float_t MAX_HADEM;
		static const Float_t MAX_TANH_DELTA_R;
		
		// enum class preferred, but not available without c++11 support
		static const Int_t EFLOW_PHOTON;
		static const Int_t ELECTRON;
		static const Int_t MUON;
		static const Int_t EFLOW_TAU;
		static const Int_t EFLOW_JET;
		static const Int_t MISSING_ET;
		static const Int_t SCALAR_HT;
		static const Int_t PYTHIA_SINGLET_JET;
		static const Int_t PYTHIA_PARTON_JET;
		static const Int_t TOWER_JET;
		
		static const Float_t PHOTON_MIN_PT_DEFAULT;
		static const Float_t ELECTRON_MIN_PT_DEFAULT;
		static const Float_t MUON_MIN_PT_DEFAULT;
		static const Float_t TAU_MIN_PT_DEFAULT;
		static const Float_t JET_MIN_PT_DEFAULT;
			
		// This module uses the IsConstituent field to record information for jet constituents (muons, isolated photons/electrons)
		// As such, we should create a relatively unique flag that we set the field to 
		static const Int_t ISO_PHOTON_FLAG;
		static const Int_t ISO_ELECTRON_FLAG;
		static const Int_t ISO_MUON_FLAG;
		static const Int_t NON_ISOLATED_FLAG;
		
		static const size_t BUFFER_SIZE; // The file buffer size in bytes
				
	private:
		FILE *fOutputFile;
		
		TObjArray* fAllParticles; //!
		TObjArray* fEFlowPhotons; //!
		TObjArray* fElectrons; //!
		TObjArray* fMuons; //!
		TObjArray* fEFlowJets; //!
		TObjArray* fTowerJets; //!
		TObjArray* fMissingEt; //!
		TObjArray* fScalarHt; //!
		TObjArray* fPythiaPartonJets; //!
		TObjArray* fPythiaSingletJets; //!
		
		TIterator* fItEFlowPhotons; //!
		TIterator* fItElectrons; //!
		TIterator* fItMuons; //!
		TIterator* fItEFlowJets; //!
		TIterator* fItTowerJets; //!
		TIterator* fItPythiaPartonJets; //!
		TIterator* fItPythiaSingletJets; //!

		Float_t photonMinPt;		
		Float_t electronMinPt;
		Float_t muonMinPt;
		Float_t tauMinPt;
		Float_t jetMinPt;
				
		bool subtractMuons;
		
		std::vector<LHCOParticleLine> eventRecord;
		std::map<Candidate*, std::vector<Candidate*> > jetToConstituent;
		
		UInt_t eventCounter, lineCounter, jetBegin, jetEnd;
		
		bool FillKinematicsAndPushBack(Candidate* candidate, Int_t type, Float_t ptMin = 0.);
		Float_t MotherFlavor(Candidate* candidate);
		void SubtractIsolated_CalculateJetStats(Candidate* jet);
		void AssociateConstituents(Candidate* jet, UInt_t jetLine);
		void FindClosestJets( 
			std::vector<LHCOParticleLine>::iterator beginCandidate,
			std::vector<LHCOParticleLine>::iterator endCandidate, 
			std::vector<LHCOParticleLine>::iterator beginJet, 
			std::vector<LHCOParticleLine>::iterator endJet);
				
		void ProcessEFlowPhotons();
		void ProcessElectrons();
		void ProcessMuons();
		void ProcessEFlowJets();
		void ProcessMissingEt();
		void ProcessScalarHt();
		void ProcessTowerJets();
		void ProcessPythiaPartonJets();
		void ProcessPythiaSingletJets();
		
		void WriteLines();
		
		static Float_t DeltaRSquared(std::vector<LHCOParticleLine>::const_iterator one, std::vector<LHCOParticleLine>::const_iterator two);
					
		ClassDef(LHCOWriterDirect, 1)
};

#endif
