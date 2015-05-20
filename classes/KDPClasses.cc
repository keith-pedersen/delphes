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
const char* const PythiaParticle::PYTHIA_EVENT_RECORD_BRANCH_NAME = "EventRecord";

ExRootTreeWriter* PythiaParticle::CreatePythiaOutputTree(TFile* const pythiaTreeFile)
{
	return new ExRootTreeWriter(pythiaTreeFile, PYTHIA_TREE_NAME);
}

ExRootTreeBranch* PythiaParticle::CreateInfoBranch(ExRootTreeWriter* pythiaWriter)
{
	return pythiaWriter->NewBranch(PYTHIA_EVENT_INFO_BRANCH_NAME, HepMCEvent::Class());
}

ExRootTreeBranch* PythiaParticle::CreateEventBranch(ExRootTreeWriter* pythiaWriter)
{
	return pythiaWriter->NewBranch(PYTHIA_EVENT_RECORD_BRANCH_NAME, PythiaParticle::Class());
}

ExRootTreeReader* PythiaParticle::LoadPythiaInputTree(TFile* const pythiaTreeFile)
{
	return new ExRootTreeReader(static_cast<TTree*>(pythiaTreeFile->Get(PYTHIA_TREE_NAME)));
}

TClonesArray* PythiaParticle::GetInfoBranch(ExRootTreeReader* pythiaReader)
{
	return pythiaReader->UseBranch(PYTHIA_EVENT_INFO_BRANCH_NAME);
}

TClonesArray* PythiaParticle::GetEventBranch(ExRootTreeReader* pythiaReader)
{
	return pythiaReader->UseBranch(PYTHIA_EVENT_RECORD_BRANCH_NAME);
}

// This routine fills from the Pythia event record to the output tree
void PythiaParticle::WritePythiaToTree(Long64_t eventCounter, Pythia8::Pythia* const pythia,
	ExRootTreeBranch* const eventInfoBranch, ExRootTreeBranch* const eventParticleBranch)
{
	HepMCEvent *event = static_cast<HepMCEvent *>(eventInfoBranch->NewEntry());

	event->Number = eventCounter;

	event->ProcessID = pythia->info.code();
	event->MPI = 1; //multi-particle interactions are on?
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

	// Run through the Pythia event record, including the system line
	// This makes sure that 0 is a completely invalid mother/daughter
	for(int i = 0; i < pythia->event.size(); ++i)
	{
		Pythia8::Particle &particleIn = pythia->event[i];
		PythiaParticle* const particleOut = static_cast<PythiaParticle*>(eventParticleBranch->NewEntry());

		particleOut->PID = particleIn.id();
		particleOut->PythiaStatus = particleIn.status();
		particleOut->IsPU = 0; // WARNING: Currently no pileup support

		particleOut->M1 = particleIn.mother1();
		particleOut->M2 = particleIn.mother2();

		particleOut->D1 = particleIn.daughter1();
		particleOut->D2 = particleIn.daughter2();

		particleOut->Charge = particleIn.charge();

		particleOut->Mass = particleIn.m();

		particleOut->E  = particleIn.e();
		particleOut->Px = particleIn.px();
		particleOut->Py = particleIn.py();
		particleOut->Pz = particleIn.pz();

		particleOut->T = particleIn.tProd();
		particleOut->X = particleIn.xProd();
		particleOut->Y = particleIn.yProd();
		particleOut->Z = particleIn.zProd();

		particleOut->CTau = particleIn.tau();
	}
}

// This routine fills from the Pythia event record to the output tree
void PythiaParticle::FillCandidatesFromPythiaTree(DelphesFactory* factory,
	const Float_t unitWeight, ExRootTreeBranch* eventInfo, TObjArray* allParticleOutputArray,
	TClonesArray* pythiaEventInfoArray, TClonesArray* pythiaEventRecordArray)
{
	HepMCEvent
		*copy = static_cast<HepMCEvent *>(eventInfo->NewEntry()),
		*original = static_cast<HepMCEvent *>(pythiaEventInfoArray->At(0));

	*copy = *original; // Default copy of event info
	if(unitWeight > 0.)
		copy->Weight = unitWeight; // Override weight with supplied weight

	PythiaParticle const* pythiaParticle;
	TIterator* itEventRecord = pythiaEventRecordArray->MakeIterator();

	// Run through the Pythia event record, including the system line
	// This makes sure that 0 is a completely invalid mother/daughter
	while((pythiaParticle = static_cast<PythiaParticle*>(itEventRecord->Next())))
	{
		Candidate* const candidate = factory->NewCandidate();
		allParticleOutputArray->Add(candidate);

		candidate->PID = pythiaParticle->PID;
		candidate->Status = (pythiaParticle->PythiaStatus > 0) ? 1 : ((pythiaParticle->PythiaStatus < 0) ? 2 : 0); // WARNING: Assume the PythiaStatus != 0

		candidate->M1 = pythiaParticle->M1;
		candidate->M2 = pythiaParticle->M2;

		candidate->D1 = pythiaParticle->D1;
		candidate->D2 = pythiaParticle->D2;

		candidate->Charge = pythiaParticle->Charge;
		candidate->Mass = pythiaParticle->Mass;
		candidate->CTau = pythiaParticle->CTau;

		candidate->IsPU = pythiaParticle->IsPU; // WARNING: Currently is always zero

		(candidate->Momentum).SetPxPyPzE(pythiaParticle->Px, pythiaParticle->Py, pythiaParticle->Pz, pythiaParticle->E);
		(candidate->Position).SetXYZT(pythiaParticle->X, pythiaParticle->Y, pythiaParticle->Z, pythiaParticle->T);
	}

	delete itEventRecord;
}

CompBase* TaggingEfficiencyJet::fgCompare = CompPT<TaggingEfficiencyJet>::Instance();
CompBase* TaggingEfficiencyMuon::fgCompare = CompPT<TaggingEfficiencyMuon>::Instance();

const Float_t TaggingEfficiencyMuon::Mass = 0.1056583715;

TLorentzVector TaggingEfficiencyJet::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}

TLorentzVector TaggingEfficiencyMuon::P4() const
{
  TLorentzVector vec;
  vec.SetPtEtaPhiM(PT, Eta, Phi, Mass);
  return vec;
}
