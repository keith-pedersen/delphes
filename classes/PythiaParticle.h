#ifndef PYTHIA_PARTICLE
#define PYTHIA_PARTICLE

#include "TObject.h"

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

		ClassDef(PythiaParticle, 1)
};

#endif
