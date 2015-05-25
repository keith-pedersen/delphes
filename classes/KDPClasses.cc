#include "KDPClasses.h"
#include "TLorentzVector.h"
#include "DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "Pythia8/Pythia.h"
#include "DelphesFactory.h"
#include "TObjArray.h"
#include "TClonesArray.h"

const char* const PythiaParticle::PYTHIA_TREE_NAME = "Pythia8";
const char* const PythiaParticle::PYTHIA_EVENT_INFO_BRANCH_NAME = "EventInfo";
const char* const PythiaParticle::PYTHIA_EVENT_PARTICLE_BRANCH_NAME = "EventRecord";

//------------------------------------------------------------------------------

ExRootTreeWriter* PythiaParticle::NewPythiaWriter(TFile* const pythiaTreeFile)
{
	return new ExRootTreeWriter(pythiaTreeFile, PYTHIA_TREE_NAME);
}

//------------------------------------------------------------------------------

ExRootTreeBranch* PythiaParticle::CreateInfoBranch(ExRootTreeWriter* pythiaWriter)
{
	if(pythiaWriter)
		return pythiaWriter->NewBranch(PYTHIA_EVENT_INFO_BRANCH_NAME, HepMCEvent::Class());
	else return 0;
}

//------------------------------------------------------------------------------

ExRootTreeBranch* PythiaParticle::CreateParticleBranch(ExRootTreeWriter* pythiaWriter)
{
	if(pythiaWriter)
		return pythiaWriter->NewBranch(PYTHIA_EVENT_PARTICLE_BRANCH_NAME, PythiaParticle::Class());
	else return 0;
}

//------------------------------------------------------------------------------

ExRootTreeReader* PythiaParticle::NewPythiaReader(TFile* const pythiaTreeFile)
{
	if(pythiaTreeFile)
		return new ExRootTreeReader(static_cast<TTree*>(pythiaTreeFile->Get(PYTHIA_TREE_NAME)));
	else return 0;
}

//------------------------------------------------------------------------------

TClonesArray* PythiaParticle::GetInfoBranch(ExRootTreeReader* pythiaReader)
{
	if(pythiaReader)
		return pythiaReader->UseBranch(PYTHIA_EVENT_INFO_BRANCH_NAME);
	else return 0;
}

//------------------------------------------------------------------------------

TClonesArray* PythiaParticle::GetParticleBranch(ExRootTreeReader* pythiaReader)
{
	if(pythiaReader)
		return pythiaReader->UseBranch(PYTHIA_EVENT_PARTICLE_BRANCH_NAME);
	else return 0;
}

//------------------------------------------------------------------------------

// This routine fills from the Pythia event record to the output tree
void PythiaParticle::WritePythia(const Long64_t eventCounter, Pythia8::Pythia* const pythia,
	ExRootTreeBranch* const eventParticleBranch, ExRootTreeBranch* const eventInfoBranch, const Bool_t isPU)
{
	WriteInfo(eventCounter, pythia, eventInfoBranch);

	if(eventParticleBranch)
	{
		// Run through the Pythia event record, including the system line (index = 0)
		// This makes sure that 0 is a valid, but nonsensical, mother/daughter
		for(int i = 0; i < pythia->event.size(); ++i)
		{
			static_cast<PythiaParticle*>(eventParticleBranch->NewEntry())->Initialize(pythia->event[i], isPU);
		}
	}
}

//------------------------------------------------------------------------------

void PythiaParticle::WriteInfo(const Long64_t eventCounter,
	Pythia8::Pythia* const pythia, ExRootTreeBranch* const eventInfoBranch)
{
	if(eventInfoBranch)
	{
		HepMCEvent *event = static_cast<HepMCEvent *>(eventInfoBranch->NewEntry());

		event->Number = eventCounter;

		event->ProcessID = pythia->info.code();
		event->MPI = pythia->flag("PartonLevel:MPI"); // For some stupid reason, this function (Pythia::flag) is not const. This is the only reason we can't use a (Pythia const* const ) throughout
		event->Weight = pythia->info.weight();
		event->Scale = pythia->info.QRen();
		event->AlphaQED = pythia->info.alphaEM();
		event->AlphaQCD = pythia->info.alphaS();

		event->ID1 = pythia->info.id1();
		event->ID2 = pythia->info.id2();
		event->X1 = pythia->info.x1();
		event->X2 = pythia->info.x2();
		event->ScalePDF = pythia->info.QFac();
		event->PDF1 = pythia->info.pdf1();
		event->PDF2 = pythia->info.pdf2();

		// WARNING: Currently no stopwatch support
		event->ReadTime = 0.;
		event->ProcTime = 0.;
	}
}

//------------------------------------------------------------------------------

void PythiaParticle::Initialize(const Pythia8::Particle& original, const Bool_t isPU)
{
	PID = original.id();
	PythiaStatus = original.status();
	IsPU = isPU;
	Charge = original.charge();

	M1 = original.mother1();
	M2 = original.mother2();

	D1 = original.daughter1();
	D2 = original.daughter2();

	E  = original.e();
	Px = original.px();
	Py = original.py();
	Pz = original.pz();

	T = original.tProd();
	X = original.xProd();
	Y = original.yProd();
	Z = original.zProd();

	Mass = original.m();
	CTau = original.tau();
}

//------------------------------------------------------------------------------

/*
// This routine fills from the Pythia event record to the output tree
void PythiaParticle::FillDelphesFactoryFromPythiaTree(DelphesFactory* const factory,
	TObjArray* const allParticleOutputArray, TClonesArray const* const pythiaEventRecordArray,
	ExRootTreeBranch* const eventInfo, TClonesArray const* const pythiaEventInfoArray, const Float_t unitWeight)
{
	PythiaParticle::FillInfo(eventInfo, pythiaEventInfoArray, unitWeight);

	PythiaParticle const* pythiaParticle;
	TIterator* itEventRecord = pythiaEventRecordArray->MakeIterator();

	while((pythiaParticle = static_cast<PythiaParticle*>(itEventRecord->Next())))
	{
		Candidate* const candidate = factory->NewCandidate();
		allParticleOutputArray->Add(candidate);

		PythiaParticle::FillCandidate(candidate, pythiaParticle);
	}

	delete itEventRecord;
}
*/

//------------------------------------------------------------------------------

void PythiaParticle::CopyInfo(ExRootTreeBranch* const eventInfo, TClonesArray const* const pythiaEventInfoArray, const Float_t unitWeight)
{
	if(eventInfo && pythiaEventInfoArray)
	{
		HepMCEvent
			*copy = static_cast<HepMCEvent *>(eventInfo->NewEntry()),
			*original = static_cast<HepMCEvent *>(pythiaEventInfoArray->At(0));

		*copy = *original; // Default copy of event info
		if(unitWeight > 0.)
			copy->Weight = unitWeight; // Override weight with supplied weight
	}
}

//------------------------------------------------------------------------------

void PythiaParticle::FillCandidate(Candidate* const candidate, const Int_t indexShift) const
{
	candidate->PID = PID;
	candidate->Status = (PythiaStatus > 0) ? 1 : ((PythiaStatus < 0) ? 2 : 0);
	candidate->IsPU = (IsPU == kFALSE) ? 0 : 1;
	candidate->Charge = Charge;

	candidate->M1 = M1 + indexShift;
	candidate->M2 = M2 + indexShift;

	candidate->D1 = D1 + indexShift;
	candidate->D2 = D2 + indexShift;

	(candidate->Momentum).SetPxPyPzE(Px, Py, Pz, E);
	(candidate->Position).SetXYZT(X, Y, Z, T);

	candidate->Mass = Mass;
	candidate->CTau = CTau;
}

//------------------------------------------------------------------------------

CompBase* TaggingEfficiencyJet::fgCompare = CompPT<TaggingEfficiencyJet>::Instance();
CompBase* TaggingEfficiencyMuon::fgCompare = CompPT<TaggingEfficiencyMuon>::Instance();

//------------------------------------------------------------------------------

const Float_t TaggingEfficiencyMuon::Mass = 0.1056583715;

//------------------------------------------------------------------------------

TLorentzVector TaggingEfficiencyJet::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

//------------------------------------------------------------------------------

TLorentzVector TaggingEfficiencyMuon::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}
