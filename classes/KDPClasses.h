// Created this file here to re-use the Delphes machinery to generate the Root dictionary

#ifndef KDP_CLASSES
#define KDP_CLASSES

#include "TObject.h"
#include "classes/SortableObject.h"

class TLorentzVector;
class CompBase;
class ExRootTreeBranch;
class ExRootTreeWriter;
class ExRootTreeReader;
namespace Pythia8
{
	class Pythia;
}
class DelphesFactory;
class TObjArray;
class TClonesArray;
class Candidate;

// KDP
// This class is designed to store information from a Pythia run
// This allows the Pythia run to take place once, not every time
// the detector is run.
class PythiaParticle : public TObject
{
	public:
		typedef Float_t   Kinematic_t;
		typedef UShort_t  Index_t;
		typedef Int_t     ID_t;
		typedef Short_t   Status_t;
		typedef Bool_t    PU_t;
		typedef Char_t    Charge_t;

		ID_t PID; // PDG ID number

		Status_t PythiaStatus; // Pythia status
		PU_t IsPU; // 1 for particles from pile-up interactions (0 for non-PU)
		Charge_t Charge; // particle charge

		// Indices of mother and daughters in the event record
		Index_t M1; // particle 1st mother | hepevt.jmohep[number][0] - 1
		Index_t M2; // particle 2nd mother | hepevt.jmohep[number][1] - 1

		Index_t D1; // particle 1st daughter | hepevt.jdahep[number][0] - 1
		Index_t D2; // particle last daughter | hepevt.jdahep[number][1] - 1

		// Mass, Energy, Mometnum stored in GeV
		Kinematic_t E; // particle energy | hepevt.phep[number][3]
		Kinematic_t Px; // particle momentum vector (x component) | hepevt.phep[number][0]
		Kinematic_t Py; // particle momentum vector (y component) | hepevt.phep[number][1]
		Kinematic_t Pz; // particle momentum vector (z component) | hepevt.phep[number][2]

		// Position and CTau stored in millimeters
		Kinematic_t T; // particle vertex position (t component) | hepevt.vhep[number][3]
		Kinematic_t X; // particle vertex position (x component) | hepevt.vhep[number][0]
		Kinematic_t Y; // particle vertex position (y component) | hepevt.vhep[number][1]
		Kinematic_t Z; // particle vertex position (z component) | hepevt.vhep[number][2]

		Kinematic_t Mass; // particle mass
		Kinematic_t CTau; // The actual decay lifetime (in the rest frame)

		static const char* const PYTHIA_TREE_NAME;
		static const char* const PYTHIA_EVENT_INFO_BRANCH_NAME;
		static const char* const PYTHIA_EVENT_RECORD_BRANCH_NAME;

		static ExRootTreeWriter* CreatePythiaOutputTree(TFile* const pythiaTreeFile);
		static ExRootTreeBranch* CreateInfoBranch(ExRootTreeWriter* pythiaWriter);
		static ExRootTreeBranch* CreateEventBranch(ExRootTreeWriter* pythiaWriter);

		static ExRootTreeReader* LoadPythiaInputTree(TFile* const pythiaTreeFile);
		static TClonesArray* GetInfoBranch(ExRootTreeReader* pythiaReader);
		static TClonesArray* GetEventBranch(ExRootTreeReader* pythiaReader);

		void WritePythiaToTree(const Long64_t eventCounter, Pythia8::Pythia const* const pythia,
			ExRootTreeBranch* const eventInfoBranch = 0, ExRootTreeBranch* const eventParticleBranch = 0, const Bool_t isPU = kFALSE);

		void FillDelphesFactoryCandidatesFromPythiaTree(DelphesFactory* const factory,
			TObjArray* const allParticleOutputArray, TClonesArray const* const pythiaEventRecordArray,
			ExRootTreeBranch* const eventInfo = 0, TClonesArray const* const pythiaEventInfoArray = 0, const Float_t unitWeight = 0.);

		void FillExternalCandidatesFromPythiaTree(ExRootTreeBranch* const newCandidateBranch,
			TObjArray* const allParticleOutputArray, TClonesArray const* const pythiaEventRecordArray,
			ExRootTreeBranch* const eventInfo = 0, TClonesArray const* const pythiaEventInfoArray = 0, const Float_t unitWeight = 0.);

		static void FillInfo(ExRootTreeBranch* const eventInfo, TClonesArray const* const pythiaEventInfoArray, const Float_t unitWeight);
		static void FillCandidate(Candidate* const candidate, PythiaParticle const* const pythiaParticle);

		ClassDef(PythiaParticle, 2)
};

//------------------------------------------------------------------------------

class TaggingEfficiencyJet : public SortableObject
{
	public:
		typedef Float_t Kinematic_t;
		typedef UChar_t Tag_t;

		Kinematic_t PT; // jet transverse momentum
		Kinematic_t Eta; // jet pseudorapidity
		Kinematic_t Phi; // jet azimuthal angle

		Kinematic_t Mass; // jet invariant mass
		Kinematic_t EhadOverEem; // ratio of the hadronic versus electromagnetic energy deposited in the calorimeter

		Kinematic_t HardCoreRatio; // (Stored in Candidate::FracPt[0])
		Kinematic_t MinCoreRatio;  // (Stored in Candidate::FracPt[1])

		Kinematic_t xMinHardCore; // The smallest x for a muon with the HardCore (Stored in Candidate::Tau[0])
		Kinematic_t xMinMinCore;  // The smallest x for a muon with the MinCore (Stored in Candidate::Tau[1])

		// Candidate::TauTag is used to indicate the presence of goodMuons

		Tag_t BTag; // >0 for a jet that has been tagged as containing a heavy quark
		//Tag_t MatriarchFlavor; // 0 or 1 for a jet that has been tagged as containing a heavy quark

		TRefArray Muons; // references to muons

		static CompBase *fgCompare; //!
		const CompBase *GetCompare() const { return fgCompare; }

		TLorentzVector P4() const;

		ClassDef(TaggingEfficiencyJet, 2)
};



class TaggingEfficiencyMuon: public SortableObject
{
	public:
		typedef Float_t Kinematic_t;
		typedef Char_t  Charge_t;
		typedef UInt_t ID_t;

		Kinematic_t PT; // muon transverse momentum
		Kinematic_t Eta; // muon pseudorapidity
		Kinematic_t Phi; // muon azimuthal angle

		Kinematic_t ImpactParameter; // radial distance of track close approach to beamline (signed by L_z)

		// x to various objects
		Kinematic_t xHardCore; // (Stored in Candidate::Tau[0])
		Kinematic_t xMinCore; //  (Stored in Candidate::Tau[1])
		Kinematic_t xTrue; //     (Stored in Candidate::Tau[2])

		// change in angle to matriarch when you add the muon a second time
		// to estimate the nuetrino
		Kinematic_t deltaTheta2MuHardCore; // (Stored in Candidate::FracPt[0])
		Kinematic_t deltaTheta2MuMinCore; //  (Stored in Candidate::FracPt[1])

		ID_t MotherID; // Mother is what emitted the muon (if tau, find tau's mother) (Stored in Candidate::BTag)
		ID_t MatriarchID; // Matriarch is primary hadron                              (Stored in Candidate::TauTag)

		Charge_t Charge; // muon charge

		const static Float_t Mass;

		static CompBase *fgCompare; //!
		const CompBase *GetCompare() const { return fgCompare; }

		TLorentzVector P4() const;

		ClassDef(TaggingEfficiencyMuon, 1)
};

#endif
