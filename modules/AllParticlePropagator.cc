#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/KDPClasses.h"
#include "TLorentzVector.h"
#include "Pythia8/Pythia.h"
#include <cmath>
#include "AllParticlePropagator.h"

#include <cstdio>

//#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "TClonesArray.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"

const Double_t c_light = 299792458.; // [m/s]
const Double_t PI = acos(-1.);
const Double_t ABSURDLY_LARGE = pow(2., 71.);

const Double_t GeV_to_eV = 1.0E9;
const Double_t METERS_to_mm = 1E3;
//const Double_t mm_to_METERS = 1./METERS_to_mm;
//const Double_t CM_to_mm = 1E1;
//const Double_t microns_to_MM = 1E-3;
//const Double_t DEG_to_rad = 2*PI/360.;


// This function assumes that a > b > c > 0. If these conditions fail, then
// the result will not be as accurate as possible
Double_t KahanTriangleAreaPreSorted(const Double_t a, const Double_t b, const Double_t c)
{
	// EXTRA PARENTHESIS ARE DELIBERATE, DO NOT REMOVE
	return sqrt((a + (b + c))*(c - (a - b))*(c + (a - b))*(a + (b - c)))/4.;
	// See "Mathematics Written in Sand" by W. Kahan, bottom of page 10
	// http://www.cs.berkeley.edu/~wkahan/MathSand.pdf#page=10

	// My alternate form, with no double subtraction. Tests to be only a teensy
	// bit better (with limited data set), so not really worth the risk.
	//return sqrt((a + (b + c))*(a + (b - c))*(c + (b - a))*(a + (c - b)))/4.;
}

inline Double_t Squared(const Double_t arg) {return arg*arg;}

// This code deliberately eschews TMath

AllParticlePropagator::AllParticlePropagator():
	fRadius(0.), fRadius2(0.), fHalfLength(0.), 	fBz(0.), fMinTrackLength(0.), fMeanPileup(0.),
	fInputList(), currentInputArray(0),
	pileupStore(), pileupArraysToStore(), pileupFactory(0),
	fOutputArray(0), fChargedHadronOutputArray(0), fElectronOutputArray(0),	fMuonOutputArray(0)
{}
AllParticlePropagator::~AllParticlePropagator() {}

//------------------------------------------------------------------------------

void AllParticlePropagator::Init()
{
	// Find and import input arrays
	{
		ExRootConfParam inputArrayNames = GetParam("InputArray");
		const Long_t size = inputArrayNames.GetSize();

		for(int i = 0; i < size; ++i)
		{
			fInputList.push_back(ImportArray(inputArrayNames[i].GetString()));
		}
	}

	fRadius = GetDouble("Radius_m", 1.0) * METERS_to_mm; // Default = 1m
	fRadius2 = fRadius*fRadius;

	fHalfLength = GetDouble("HalfLength_m", 3.0) * METERS_to_mm; // Default = 3m

	fBz = GetDouble("Bz_Tesla", 0.0); // Default = 0 Tesla

	fMinTrackLength = GetDouble("MinTrackLengthRatio", .7) * fRadius; // Default = 70% of radius (~high purity tracks, minimal double-counting of energy)

	// Sanity tests
	// When these fail, the output arrays won't be created. So any subsequent module
	// which tries to link to this module during initialization will cause a run-time exception to be thrown.
	if(fRadius < 10.)
	{
		cout << "ERROR: magnetic field radius is too low\n";
		return;
	}
	if(fHalfLength < 10.)
	{
		cout << "ERROR: magnetic field length is too low\n";
		return;
	}

	// Create output arrays
	fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));
	fChargedHadronOutputArray = ExportArray(GetString("ChargedHadronOutputArray", "chargedHadrons"));
	fElectronOutputArray = ExportArray(GetString("ElectronOutputArray", "electrons"));
	fMuonOutputArray = ExportArray(GetString("MuonOutputArray", "muons"));

	// Deal with pileup
	// The fastest way to do this is to propagate all the pileup during initialization,
	// storing the resulting Candidates in sorted arrays. Then, during each event, the
	// stored pileup arrays can be copied to the master array for that event
	fMeanPileup = GetDouble("MeanPileup", 0.);

	PropagateAndStorePileup(std::string(GetString("PileupFile", "")));
}

//------------------------------------------------------------------------------

void AllParticlePropagator::Finish()
{
	// Delete all the Pileup Candidates
	delete pileupFactory;
}

//------------------------------------------------------------------------------

// Process() loops through all the Particles in the InputArrays, propagating
// any which have a non-zero lifetime.
void AllParticlePropagator::Process()
{
	// Loop through all input arrays
	for(vector<TObjArray*>::iterator itInputArray = fInputList.begin(); itInputArray != fInputList.end(); ++itInputArray)
	{
		currentInputArray = *itInputArray;
		TIterator* const itInput = currentInputArray->MakeIterator();

		// Loop through all particles sequentially
		while(Propagate(static_cast<Candidate*>(itInput->Next()), 0));
		// Once Propagate finds a particle that bends, it will recursively propagate
		// the entire decay chain, skipping ahead in the event record (but without
		// affecting itInput). Hence, once Propagate finally returns, the itInput loop
		// is guarenteed to iterate over particles which have already been propagated.
		// This is why Propagate always checks to see if TrackLength < 0 first (as
		// it will only be < 0 when the particle has yet to be Propagted).

		delete itInput;
	}

	// Add in a random number of pileup events
	FillPileup();
}

//------------------------------------------------------------------------------

// This function rotates the daughter's 4-momentum and sets its position
Candidate* AllParticlePropagator::Rotate(Candidate* const daughter, RotationXY const* const rotation, const TLorentzVector& mothersPosition)
{
	if(daughter)
	{
		daughter->Position = mothersPosition; // The daughter's creation vertex is the mother's decay vertex
		if(rotation)
			rotation->ForwardRotateDaughterMomentum(daughter); // Rotate the daughter's 4-momentum
	}

	return daughter;
}

//------------------------------------------------------------------------------

// Since AllParticlePropagator is intended to be used with longer decay chains than are
// required to get to the cylinder, once a particle makes it to the cylinder's edge,
// all of its decendents should be recursively flagged as null and "already propagated"
void AllParticlePropagator::NullifyDaughters(Candidate* const mother, const bool nullifyThisParticle)
{
	if(mother)
	{
		if(nullifyThisParticle) // We don't want to alter the mother that became final state
		{
			mother->Status = 0; // This particle is null; the decay never happened
			mother->TrackLength = 0.; // Make it appear that the particle was propagated.
		}

		// Now nullify any of it's daughters
		Int_t daughterIndex = mother->D1;

		if(daughterIndex > 0)
		{
			while(daughterIndex <= mother->D2)
			{
				NullifyDaughters(static_cast<Candidate*>( currentInputArray->At( daughterIndex++ ) ), true);
			}
		}
	}
	else
	{
		// NullifyDaughters is only called by Propagate or NullifyDaughters
		// To get here with an invalid mother pointer means that the calling instance of
		// NullifyDaughters looked up a daughter index that didn't correspond to a valid Candidate
		throw runtime_error("(AllParticlePropagator::NullifyDaughters): Daughter index does not correspond to valid Candidate*! Did you re-index?");
		// For more information, see #5 in the class description in the header file
	}
}

//------------------------------------------------------------------------------

// IMPORTANT conventions:
//
// 1) There will be less math/conversion factors if we store momentum in its dimensionless form
//         beta  ==  v / c  ==  p / energy
//    The complementary decision is to store time as a distance (ct)
//         dist  ==  v * t
//         dist  ==  beta * ct
//
// 2) Instead of converting stored distances from mm to meters, then back to mm once we're done,
//    we will keep all positions and times in mm
//
// 3) Vector notation: LaTeX typsetting is universal, but more difficult to read in its raw form
//    (and less dense). To keep comments shorter, a hybrid notation is used:
//
//       a) All vectors must use a lower-case letter for their main symbol (see f).
//
//       b) Unit vectors will be denoted with (^), with a space required after the carrot (so they are not confused with exponents):
//             x^   y^   z^
//       c) As such, exponents will be denoted in FORTRAN style
//             x**2 == pow(x, 2)
//       d) Non-unit vectors will be denoted with (>):
//             r>     p>
//       e) Subscripts will be denoted with (_), AFTER any vector symbols and without curly braces:
//             r>_helix     r0>hxPr
//       f) Magnitude is denoted by losing the vector symbol and switching to a captial letter:
//             |r>| == R, |r>_creation|**2 == R_creation**2
//       * Vector/cross product is lower-case ( x ), with a required space on either side):
//             l> == r> x p>    x^ x y^ = z^
//       * Scalar/dot product is period (.), no space required:
//             |r>|**2 = R**2 = r>.r>
//       * Normal multiplication is still (*):
//             x**2 = x*x


// This function propagates a particle in a homogenous magentic field parallel to z^
// It returns the VALIDITY of the candidate pointer (NOT whether the candidate was propagated)

// Whenever there is an active rotation (rotation != 0), all of this particle's daughters
// will be propagated recursively. If the particle rotates, then they will be transformed with
// a new rotation (which accumulates on top of the incoming rotation); otherwise, the will
// be transformed with the unaltered, incoming rotation.
bool AllParticlePropagator::Propagate(Candidate* const candidate,	RotationXY const* rotation)
{
	if(not candidate) // Test to make sure the candidate is dereferenceable
	{
		if(rotation)
		{
			// Propagate is only called by Propagate() or Process(). To get here with an
			// invalid candidate pointer, but a valid rotation, means that the calling instance of
			// Propagate() looked up a daughter index that didn't correspond to a valid Candidate
			throw runtime_error("(AllParticlePropagator::Propagate): Daughter index does not correspond to valid Candidate*! Did you re-index?");
			// For more information, see #5 in the class description in the header file
		}
		else return false;
	}

	// (candidate->TrackLength < 0) is a flag that indicates a candidate was not yet processed
	if(candidate->TrackLength >= 0.)
	{
		return true; // If the candidate was already processed, return immedietely
		// This is supposed to happen; see the long comment in Process().
	}

	// Whenever the particle is part of an active decay chain (i.e. fragmentation over).
	// all descendents of this particle will be propagated recursively. This can create
	// a lot of nested calls to Propagate. In order to keep memory space tidy, we'll keep
	// all temporary/unneccessary variables in a scope that will close before the recursive call.

	// If this particle rotates, then it will need to create a NEW rotation to use
	// on its daughters. Once all the daughters are propagated, this memory will
	// need to be freed. Hence, newRotation stores whether this particle
	// should "delete rotation" as its final task.
	bool newRotation = false;

	// Begin temporary variable scope
	{
		// Assume the particle's proper lifetime is stored in [mm]. If the particle is
		// final state (Status == 1), make cTau absurdly large. This will help later,
		// when we find the shortest of two times (decay or cylinder strike). It also
		// allows for the propagation of MASSLESS final state particles, whose cTau = 0.
		const Double_t cTau = ((candidate->Status == 1) ? ABSURDLY_LARGE : candidate->CTau);  //[mm]

		TLorentzVector& position = candidate->Position; // creation vertex
		// Store the particle's creation vertex in its Area field (see #8 in class descrption in header file)
		(candidate->Area) = position;
		const Double_t
			R02 = Squared(position.X()) + Squared(position.Y()),
			z0 = position.Z();

		{
			const Double_t creationRadius = sqrt(R02);

			// Check that the particle was created inside the cylinder
			if((creationRadius >= fRadius) or (abs(z0) >= fHalfLength))
			{
				stringstream message;
				message << "(AllParticlePropagator::Propagate): Particle (" << candidate->PID << ") created outside the cylinder ";
				message << "( " << creationRadius << " , " << z0 << " )!\n";
				message << "Did you include the entire decay chain? Did you re-index?\n";
				throw runtime_error(message.str());
				// For more information on this error, see #5 in the class description in the header file
			}

			// candidate::Creation radius is stored as a Float_t; check for sanity, then store, otherwise we could have precision conversion problems
			candidate->CreationRadius = creationRadius;
		}
                
                // Make sure the particle has energy.
                // This check was in response to the side-effectsof Pythia forcing a Particle with no energy to decay.
                // It's daughters subsequently have (nan) for their 4-momenta. This causes all sorts of havoc downstream.
                if(candidate->Momentum.E() == 0.)
                {
                        NullifyDaughters(candidate, true); // Nullify this particle and any daughters
                        return true; // Return to the calling instance immedietely
                }

		if(cTau <= 0.) // The particle decayed instantaneously (or has an invalid stored cTau)
		{
			candidate->TrackLength = 0.; // This indicates that the particle has been "processed
			// Nothing else to do
		}
		else
		{
			const TLorentzVector& momentum = candidate->Momentum; // 4-momentum

			// We'll be using cylindrical coordinates (r, z) for obvious reasons
			// Assume 4-position is stored in [mm]
			// Assume 4-momentum & mass is stored in [GeV]
			// Assume charge is stored in elementary charged [e]
			const VecXY r0(position.X(), position.Y());
			const Double_t
				energy = momentum.E(),
				q = candidate->Charge;

			// It's easier to work with beta than momentum
			const VecXY r0Beta( momentum.Px()/energy, momentum.Py()/energy );
			const Double_t R0Beta2 = r0Beta.Norm2(); // Transverse beta (used throughout instead of pT)
			const Double_t z0Beta = momentum.Pz()/energy;

			// The next 3 variables will be altered by the propagating routine, when the solutions are found.

			// "decays" lets us know whether or not the particle decays. Currently, we will assume it does.
			bool decays = true;

			// The total lab propagation time; currently, that's the time until decay (cTau * gamma).
			// For massless particles, this will be (inf). That's OK, because we will find the shortest time.
			// Use the absolute value JUST IN CASE the mass (or CTau) is accidently -0, since
			// we will be altering ctProp only if we find a shorter/smaller time
			Double_t ctProp = abs(cTau * (energy / candidate->Mass));

			// The particle's final position; default to initial position.
			VecXY rFinal(r0);

			// Figure out how long till the particle hits the endcap.
			{
				const Double_t ctEndcap = (copysign(fHalfLength, z0Beta) - z0) / z0Beta;
				// copysign() ensures that, if (z0Beta = (+/-)0.), we'll get (ctEndcap = inf)
				// (IFF the particle is inside the cylinder, which we already verified)

				if(ctEndcap <= ctProp)
				{
					ctProp = ctEndcap;
					decays = false; // Does not decay, exits the cylinder
				}
			}

			// Next, find any particle which CAN propagate helicly
			//    1. |q| != 0
			//    2. |R0Beta2| != 0
			if((q * R0Beta2) not_eq 0.)
			{
				// In order to model helical propagation, we'll need to know the gyration frequency:
				//
				//    omega = -q fBz /(gamma m)     [rad/s]
				//
				// Here, omega uses the sign convention of a RH xyz coordinate system:
				//
				//    x^ = cos(0)   and    y^ = sin(Pi/2)   with     z^ parallel to B>
				//
				// Thus, the LH rotation of a positively charged (q>0) particle is parameterized
				// by a steadily  *decreasing*  angle (a negative frequency).
				//
				// Since we use (ct) instead of (t), it is easier to work with omega/c:
				//
				//    deltaPhi == omega*t == (omega/c)*ct
				//
				//    omega / c == -q fBz /(gamma m c) == -q fBz / (energy / c)  [rad/m]
				//    omega / c == -q fBz c / energy / METERS_to_mm              [rad/mm]
				//
				const Double_t omegaOverC = (-q * fBz / energy) * (c_light / (GeV_to_eV * METERS_to_mm)); // [rad/mm]

				// PropagateHelicly() will propagate the particle. It has 4 returns. Its official
				// return indicates whether "rotation" (passed be reference) had an angle added
				// to it. Of  course all charged particles rotate, but rotation only needs to change when
				// the particle decays. If it strikes the cylinder, its daughters don't require
				// rotation (so don't bother adding its rotation).
				// The propagation solutions (ctProp and rFinal) are also passed by reference.
				newRotation = PropagateHelicly(candidate, decays,
					                            r0, z0,
					                            r0Beta, z0Beta,
					                            omegaOverC,
					                            rotation,
					                            ctProp, rFinal);
				decays = newRotation;
			}
			else
			{
				// We can use a transverse vector equation to solve for the straight-line propagation time.
				//
				//           rFinal> = r0> + ctBarrel * r0Beta>
				//
				// Square both sides, create the quadratic equation for ctBarrel:
				//
				//           a = R0Beta2       b = 2 * r0>.r0Beta       c = (R02 - fRadius2)
				//
				// In order to prevent catastrophic cancellation (for very small ctBarrel when b > 0),
				// we should use the two forms of the roots:
				//
				//         ctBarrel1 = -c / (b/2 + sign(b)*sqrt((b/2)**2 - ac)
				//         ctBarrel2 = -(b/2 + sign(b)*sqrt((b/2)**2 - ac))/a
				//
				// Because ct needs to be positive (and c < 0 and a > 0), we must choose:
				//
				//         ctBarrel1 when b > 0
				//         ctBarrel2 when b <= 0

				const Double_t halfB = r0.Dot(r0Beta);
				const bool bIsPositive = halfB > 0.; // Get the conditional jump in the pipeline
				const Double_t rootTerm = sqrt(halfB*halfB + (fRadius2 - R02)*R0Beta2);

				const Double_t ctBarrel =
					(bIsPositive ?
						(fRadius2 - R02)/(halfB + rootTerm) :
						(rootTerm - halfB)/R0Beta2);

				if(ctBarrel <= ctProp)
				{
					ctProp = ctBarrel;

					if(ctBarrel < 0.)
						throw runtime_error("(AllParticlePropagator::Propagate): Keith, you fat fuck, your straight-line math is all wrong!");
					decays = false;
				}

				rFinal += r0Beta*ctProp; // rFinal was initialized to initial position
			}

			// Set the final position
			position.SetXYZT(rFinal.x, rFinal.y,
				z0 + z0Beta * ctProp,
				position.T() + ctProp);

			// Check for propagation error (i.e. supposedly decays but actually outside of cylinder)
			// In reality, we should just declare that this particle does not decay
			// But, since I just found a different source of error in my code,
			// let's just keep this here for now until I'm confident that
			// we can safely make such a statement.
			if(decays and ((rFinal.Norm() >= fRadius) or (abs(position.Z()) >= fHalfLength)))
			{
				const Double_t omegaOverC = (-q * fBz / energy) * (c_light / (GeV_to_eV * METERS_to_mm)); // [rad/mm]
				printf("rFinal: %.16e\n", rFinal.Norm());
				printf("zFinal: %.16e\n", position.Z());
				printf("omega:  %.16e\n", omegaOverC);
				throw runtime_error("(AllParticleProagator::Propagate) Particle propagated to outside of cylinder!");
			}

			// Set the track length and use it to do a simple pre-sort of track candidates,
			// so that the charged hardron output array doesn't get filled with a bunch of invisible tracks
			{
				const Double_t trackLength = ctProp * sqrt(R0Beta2 + z0Beta*z0Beta);
				candidate->TrackLength = trackLength;

				if(q not_eq 0.) // Charged particles, add to tracked particle output arrays
				{
					Int_t absPID = abs(candidate->PID);

					if(absPID == 13) // Always add muons regardless of track length
					{
						// Even if they were created in last mm of tracker, they still might show up in the muon chamber
						fMuonOutputArray->Add(candidate);
					}
					else if(trackLength >= fMinTrackLength) // Add other particles only with sufficient track length
					{
						if(absPID == 11)
							fElectronOutputArray->Add(candidate);
						else
							fChargedHadronOutputArray->Add(candidate);
					}
				}
			}

			if(not decays)
			{
				// The candidate punctured the cylinder
				candidate->Status = 1; // Ensure the particle is now considered final state
				fOutputArray->Add(candidate); // Only candidates which puncture the cylinder should hit the calorimter
				NullifyDaughters(candidate); // Terminate its decay chain by nullifying its daughters
				return true; // Return to the calling instance immedietely
			}
		}//End propagation
	}//End temporary variable scope

	// Previously, we would only propagate recursively when there was an active rotation.
	// BUT, I started storing the pythia particles with binary32, while
	// propagating from the primary vertex with binary64. This means that
	// a totally nuetral decay lineage may no longer completely jive
	// (i.e. the creation vertex of a particle who's pregenetors are all
	// neutral is not neccessarily the same as their mother's decay position,
	// due to the differing precision of the calculations). This becomes
	// a problem near the edge of the cylinder (this problem
	// was discovered in 1 decay out of 10 billion). Thus, decay all
	// particles recursively when their final position is not at the primary vertex
	// (i.e. they have been propagated). This filters out colored particles.
	if(candidate->Position.Vect().Mag2() > 0.)
	{
		{
			Int_t daughterIndex = candidate->D1;

			if(daughterIndex > 0) // There are daughters
			{
				while(daughterIndex <= candidate->D2)
				{
					// The pointer to the daughter is sent to Rotate, which then returns it,
					// becoming the first argument of Propagate.
					// This prevents an unneccessary copy of each candidate from proliferating
					// at each level of recursion (see example below).
					Propagate(
						Rotate(static_cast<Candidate*>( currentInputArray->At( daughterIndex++ ) ),
							rotation, candidate->Position),
						rotation);

					// Example: The unneccessary pointer "daughter" sticks around until
					// Propagate returns, and is then promptly overwritten with the next daughter,
					// so there really was no point for it to wait in the stack, was there?
					/*
					Candidate* daughter = static_cast<Candidate*>( currentInputArray->At( daughterIndex++ );
					Rotate(daughter, rotation, candidate->Position);
					Propagte(daughter, rotation);
					*/
				}
			}
		}

		if(newRotation)
		{
			// The rotation pointed to by rotate was created for this particle
			// Now that all of its daughters have been propagated, we can delete it
			delete rotation;
		}
	}

	return true;
}

//------------------------------------------------------------------------------

// This function has 4 returns (3 via pass-by-reference)
// The actual return indicates whether (rotation) was altered (thus requiring deletion).
// Ideally, the body of this function could be contained inside Propagate, but
// keeping it as a seperate function keeps Propagate() much cleaner and easier to follow
bool AllParticlePropagator::PropagateHelicly(Candidate* const candidate, const bool decays,
	const VecXY& r0, const Double_t z0,
	const VecXY& r0Beta, const Double_t z0Beta,
	const Double_t omegaOverC,
	RotationXY const*& rotation,
	Double_t& ctProp, VecXY& rFinal)
{
	// In this function, the subscript (_) of a variable indicates its coordinate system

	// First we need to shift the origin to the center of the helix (the helix [hx] coordinate system).
	// If r0>_hx is the particle's initial position in the helix system, then we can start
	// with the initial magnetic angular momentum:
	//
	//     I * omega * z^ == r0>_helix x p>_T == r0>_helix x r0Beta> * (energy/c)
	//
	// Crossing in r0Beta> from the left, and rearranging terms (e.g. I = gamma*m*R_helix**2), we get:
	//
	//     r0>_helix == (r0Beta> x z^) / (omega / c)
	//
	const VecXY r0_hx = r0Beta.CrossZHat()/omegaOverC; // The vector from the center of the helix coordinate system to r0>
	const VecXY rBeam_hx = r0_hx - r0; // The vector from the helix origin to the beamline

	// We need the length of rBeam>_hx and r>_hx
	const Double_t
		R2_hx = r0_hx.Norm2(),
		RBeam2_hx = rBeam_hx.Norm2();
	const Double_t
		R_hx = sqrt(R2_hx),
		RBeam_hx = sqrt(RBeam2_hx);

	// If the particle isn't a looper, we'll need to calcualte it's barrel exit solution.
	// Let's get some calculations started now, that we'll need to use later.
	// Hopefully the compiler can get these in the pipe, to reduce latency when they are needed.
	const bool notLooper = (RBeam_hx + R_hx >= fRadius);
	std::vector<Double_t> smallToLargeTriangleSides;

	// For ease of working with angles, let's rotate our perspective to the helixPrime
	// coordinate system (hxPr). In this system, the helix's perigree (close approach) to
	// the beamline lies at {R_helix, 0}, with apogee at {-R_helix, 0}. We can parameterize these positions
	// using (phi), with (phi == 0) corresponding to perigree and (phi == sign(omega)*pi) corresponding
	// to agpogee.

	Double_t phi0; // The initial anglular position (in hxPr) will be set in a moment

	// Prepare for the edge intercept calculation
	if(notLooper)
	{
		// Sort the sides of the triangle (to be explained later)
		smallToLargeTriangleSides = {fRadius, R_hx, RBeam_hx};
		std::sort(smallToLargeTriangleSides.begin(), smallToLargeTriangleSides.end());
	}

	{// Find the coordinates of close approach between the track circle and the beamline.
		{// {x,y}
			// In the xy plane, the point of close approach lies on the line drawn through
			// the origin and helix center, either in front of behind the beamline. Thus,
			// the distance is determined by the difference in radii from the pertinent circles.
			Double_t RCloseApproach = R_hx - RBeam_hx;
			// The sign is currently backwards (negative if CloseAppraoch is between beam
			// and helix center), this is for good reason (we will correct it in just a moment).
			candidate->Dxy = RCloseApproach;

			// Project out the signed x and y components of the close approach position.
			// No trig required, rBeam_hx is already parallel to the line connecting the
			// beam line and helix center. But because rBeam_hx points from helix center to
			// the beam-line, we need a sign flip; this was already built into RCloseApproach.

			// Scale RCloseApproach to scale the components
			RCloseApproach /= RBeam_hx;
			candidate->Xd = rBeam_hx.x * RCloseApproach;
			candidate->Yd = rBeam_hx.y * RCloseApproach;

			// By convention, the impact parameter has the same sign as the scalar product between
			// the close approach vector and jet centroid. Since we haven't built any jets yet, we can
			// use the momentum of the particle.
			candidate->Dxy = copysign(candidate->Dxy, (r0Beta.x*candidate->Xd + r0Beta.y*candidate->Yd));
		}

		{
			// In the hxPr coordinate system, the angle of closest approach is at
			// rCloseApproach_hxPr> = {R_hx, 0}, which is parllel to rBeam_hx> in
			// the helix coordinate system.  Thus, the particle's initial angular position
			// (phi0) can be found from r0>_hx and rBeam>_hx.
			//
			// The most accurate way to find the (interior) angle between two vectors is NOT:
			//
			//     acos(a>.b>/sqrt(a>.a> * b>.b>))
			//
			// W. Kahan suggests a form in "How Futile are Mindless Assessments..." on the
			// top of p. 47, but (per my personal research) the best form is:
			//
			//     atan2(sqrt(|a> x b>|^2), a>.b>)
			//
			// This form may seem like slight overkill, but the acos form loses precision
			// for exactly the particles we need extra precision for: super
			// energetic particles that don't bend very much, with small angles between
			// r0_hx> and rBeam_hx>. Such particles really need to show up
			// in the correct Calorimeter cell. Plus, the acos form can sometimes
			// send cos(x) > 1 for very parallel vectors (normalization error),
			// which producess a NaN from acos.
			//
			// The sign of phi_0 depends on which side of perigree the particle
			// starts on (positive or negative deltaPhi, where positive is RH rotation)
			//
			//    sign(phi0) = sign(rBeam>_hx x r0>_hx)
			//
			// With atan2(y, x), the sign of x determines the magnitude of the angle,
			// while the sign of y determines the sign of the angle. Hence, we will
			// give the y argument the sign of the vector product. BUT, since we are working with
			// 2D vectors, the vector product is a scalar, so we can feed it in directly

			phi0 = atan2(rBeam_hx.Cross(r0_hx), rBeam_hx.Dot(r0_hx));

			// To find z at close approach, subtract the z-distance travelled
			// since (or until) close approach.

			candidate->Zd = z0 - z0Beta * (phi0 / omegaOverC);
		}
	}

	if(notLooper) // The particle crosses the barrel
	{
		// To find where the particle cross the barrel, we can use a vector equation:
		//
		//        (rBarrel>_hx - rBeam>_hx)**2 == rBarrel>**2
		//
		// Expanding to the dot product, we find the angle between rBarrel>_hx & rBeam>_hx:
		//
		//        cos(epsilon) == (RBeam_hx**2 + R_hx**2 - fRadisu2) / (2 * R_hx * RBeam_hx)
		//
		// Which, of course, maps to opposite angles.
		//
		//        phiBarrel == +/- epsilon
		//
		// To find the correct solution (the first angle encoutered by the particle),
		// we can use a very simple thought experiment. Imagine a positively charged
		// particle moving in its LH helix (omega < 0). Since it is constantly moving
		// to smaller phi. the only way that (+epsilon) is its first encounter
		// with the barrel is if phi0 > +epsilon. But, because the way we define
		// the hx_pr system, this would also mean that the particle started
		// OUTSIDE the barrel (which we already verified is false). Since the opposite
		// argument works for the negatively charged particle, we find:
		//
		//        phiBarrel == sign(Omega)*epsilon
		//
		// Now, using sin(epsilon) = sqrt(1-cos(epsilon)**2) loses a lot of
		// precision, so we're better off using the law of sines:
		//
		//        sin(epsilon) = sign(Omega)*2*AreaOfTriangle/(R_hx * RBeam_hx)
		//
		// Thus, to find epsilon, we can use atan2.

		const Double_t denom = R_hx * RBeam_hx;

		// The Kahan Triangle area requires pre-sorted sides (a > b > c), which we already did
		const Double_t sinEpsilon =
			copysign(2. * KahanTriangleAreaPreSorted(smallToLargeTriangleSides[2], smallToLargeTriangleSides[1], smallToLargeTriangleSides[0]) / denom,
				omegaOverC);

		// For the cos, we should always subtract the two terms which are closest in:
		//    const Double_t cosEpsilon = (RBeam2_hx + R2_hx - fRadius2)/(2.*denom);
		// This reduces rounding errors by about 3%, and we can accomplish with a simple inequality
		//
		// For (A + B - C), if |A-C| < |B-C|, then we want to use (B + (A - C)) (and vice versa)
		// Using (A-C)**2 < (B-C)**2, reduces to (A < B) ? (A + B > 2C) : (A + B < 2C)
		//
		// We should also express a**2 - b**2 as (a - b)*(a + b), because it is much more accurate
		// when a and b are very close.
		// This reduces rounding errors by about 10%

		/* Old form
		const Double_t cosEpsilon = ((RBeam2_hx < R2_hx) ?
			((RBeam2_hx + R2_hx > 2.*fRadius2) ? (R2_hx + (RBeam2_hx - fRadius2)) : (RBeam2_hx + (R2_hx - fRadius2))) :
			((RBeam2_hx + R2_hx < 2.*fRadius2) ? (R2_hx + (RBeam2_hx - fRadius2)) : (RBeam2_hx + (R2_hx - fRadius2))))/(2.*denom);
		*/

		const Double_t cosEpsilon = ((RBeam2_hx < R2_hx) ?
			((RBeam2_hx + R2_hx > 2.*fRadius2) ?
				(R2_hx + (RBeam_hx - fRadius)*(RBeam_hx + fRadius)) :
				(RBeam2_hx + (R_hx - fRadius)*(R_hx     + fRadius))) :
			((RBeam2_hx + R2_hx < 2.*fRadius2) ?
				(R2_hx + (RBeam_hx - fRadius)*(RBeam_hx + fRadius)) :
				(RBeam2_hx + (R_hx - fRadius)*(R_hx     + fRadius))))/(2.*denom);

		// Rounding errors will occasionally cause some non-perfect normalization. We can
		// adjust for this easily.
		const Double_t normalization = sqrt(sinEpsilon*sinEpsilon + cosEpsilon*cosEpsilon);

		// However, we should be seeing errors of O(epsilon) (i.e. 2e-16), so check for absurd errors
		if(abs(normalization-1.) > 1e-14)
		{
			printf("\nAbnormally high normalization error\n\n");
			printf("\n\nsin: %.16e\ncos: %.16e\n(norm-1): %.16e\n", sinEpsilon, cosEpsilon, normalization-1.);
			printf("r0 (m): %.2f\n", r0.Norm()/1000.);
			printf("gamma-1: %.3e\n\n", 1./sqrt(1.-r0Beta.Norm2())-1.);
		}

		// No need to correct for normalization before atan2, it is automatically done
		// (i.e. only normalize when the particle strikes the barrel)
		const Double_t ctBarrel = (atan2(sinEpsilon, cosEpsilon) - phi0) / omegaOverC;

		if(ctBarrel <= ctProp) // The barrel exit is the closest exit
		{
			ctProp = ctBarrel;

			// Check for a numerical error (negative propagation time to barrel)
			if(ctBarrel < 0.)
			{
				stringstream message;
				char messageChar[1024];
				sprintf(messageChar, "(AllParticlePropagator::PropagateHelicly): Keith, you fat fuck, your helix math is all wrong!\nctBarrel = %.16e\n", ctBarrel);
				message << messageChar;

				sprintf(messageChar, "r0: %.16e \nHelix radius: %.16e \nBeam distance: %.16e \nCylinder radius: %.16e \nOmegaOverC: %.16e \n",
					r0.Norm(), R_hx, RBeam_hx, fRadius, omegaOverC);
				message << messageChar;

				sprintf(messageChar, "cosEpsilon: %.16e \nsinEpsilon: %.16e\n", cosEpsilon, sinEpsilon);
				message << messageChar;

				sprintf(messageChar, "phi0: %.16e \nphiBarrel: %.16e\n", phi0, atan2(sinEpsilon, cosEpsilon));
				message << messageChar;

				throw runtime_error(message.str());
			}

			const Double_t normalizedR_hx = R_hx / normalization;

			// Because we know the cos and sin of epsilon, we know where the exit vector is in hxPr
			// Creating it in hxPr then rotating back to the helix system is the best way
			rFinal.x = normalizedR_hx*cosEpsilon;
			rFinal.y = normalizedR_hx*sinEpsilon;

			// In hxPr, rBeam>_hx is at phi = 0. Thus, to convert from hxPr to hx coordinate
			// systems, we can use rBeam>_hx as our rotation vector. To get back to the
			// beam frame, simply subtract rBeam>_hx
			(RotationXY(rBeam_hx.x/RBeam_hx, rBeam_hx.y/RBeam_hx).ForwardRotate(rFinal)) -= rBeam_hx;

			// The particle struck the barrel; no daughters, no new rotation
			return false;
		}
	}

	// The particle doesn't strike the barrel; to find it's exit position,
	// we simply rotate r0>_helix, then shift back to the beam frame
	RotationXY* newRotation = new RotationXY(ctProp * omegaOverC);
	rFinal = r0_hx;
	(newRotation->ForwardRotate(rFinal)) -= rBeam_hx;

	if(decays)
	{
		// Create a new cumulative rotation for the daughters by adding the rotation of the mother
		newRotation->Add(rotation);
		rotation = newRotation; // Return (by reference) the new rotation for this particle's daughters
		return true; // Pass ownership of the (rotation) memory to the calling instance of Propagate
	}
	else
	{
		delete newRotation;
		return false;
	}
}

//------------------------------------------------------------------------------

void AllParticlePropagator::PropagateAndStorePileup(const std::string& pileupFileName)
{
	if(pileupFileName.length() == 0)
	{
		if(fMeanPileup > 0.)
			throw runtime_error("No pileup file supplied, even though MeanPileup > 0");
	}
	else
	{
		cout << "=================> Propagating and Storing Pileup <=================\n";

		TFile* pileupFile = TFile::Open(pileupFileName.c_str(), "READ");
		if(not pileupFile)
			throw runtime_error("Can't find pileup file <" + pileupFileName + ">");

		ExRootTreeReader* pileupReader = PythiaParticle::NewPythiaReader(pileupFile);
		Long64_t pileupStoreSize = pileupReader->GetEntries(); // Long64_t is the type used by ExRootTreeReader::ReadEntry()

		if(pileupStoreSize == 0)
			throw runtime_error("Can't read Pythia EventRecord from <" + pileupFileName + "> or contains no events!\n");
		else if(Double_t(pileupStoreSize) < 10. * fMeanPileup)
			throw runtime_error("Pileup file <" + pileupFileName + "> does not contain enough events for supplied mean!\n Contains = " + to_string((long long int)pileupStoreSize) + ", Mean = " + to_string((long double)fMeanPileup) + "\n");

		// Get the TClonesArray for the pileup particles
		TClonesArray* pileupEventRecord = PythiaParticle::GetParticleBranch(pileupReader);

		// Create a TObjArray for the filled Candidates (initial capacity 1024)
		currentInputArray = new TObjArray(1024);

		// Now record the arrays that will be stored in pileupStore
		//  * stableParticles
		//  * chargedHadrons
		//  * electrons
		//  * muons
		pileupArraysToStore.push_back(fOutputArray);
		pileupArraysToStore.push_back(fChargedHadronOutputArray);
		pileupArraysToStore.push_back(fElectronOutputArray);
		pileupArraysToStore.push_back(fMuonOutputArray);

		// Ready the pileupStore by filling it with a bunch of empty vectors
		pileupStore = std::vector<std::vector<std::vector<Candidate*> > >(pileupStoreSize, std::vector<std::vector<Candidate*> >(pileupArraysToStore.size()));

		DelphesFactory* factory = GetFactory();

		// Since we will not have access to the full event record, there is no sense in keeping it arround.
		// Thus, to save space, we should only keep around the Candidates that are in the output vectors.
		// (This is partially a problem because the Candidate is so bloated with Jet stuff)

		// Step1. Create a temporary Candidate factory for pileup. This is where we will copy the full event record
		ExRootTreeBranch temporaryPileupFactory("TemporaryPileupFactory", Candidate::Class());

		// Also create a permanent factory for the pileup Candidate which stick around (the Delphes factory is cleared after each event,
		// whereas the pileup Candidate's need to survive the entire run)
		pileupFactory = new ExRootTreeBranch("PileupFactory", Candidate::Class());

		// Finally, create a map from the temporary to permanent pointers
		std::map<Candidate*, Candidate*> permanentPileup;

		for(Long64_t entry = 0; entry < pileupStoreSize; ++entry)
		{
			cout << "-"; // A progress bar

			// Copy the Pythia particles to Candidates
			{
				// Read the current Pythia event record
				pileupReader->ReadEntry(entry);
				// Clear the currentInputArray
				currentInputArray->Clear();

				temporaryPileupFactory.Clear();

				PythiaParticle const* pythiaParticle;
				Candidate* candidate;
				TIterator* itEventRecord = pileupEventRecord->MakeIterator();

				while((pythiaParticle = static_cast<PythiaParticle const*>(itEventRecord->Next())))
				{
					candidate = static_cast<Candidate*>(temporaryPileupFactory.NewEntry());
					candidate->Clear();
					currentInputArray->Add(candidate);

					pythiaParticle->FillCandidate(candidate);
				}

				delete itEventRecord;
			}

			// Propagate all Pileup particles
			{
				// Clear the output arrays
				for(std::vector<TObjArray*>::iterator itArray = pileupArraysToStore.begin(); itArray not_eq pileupArraysToStore.end(); ++itArray)
					(*itArray)->Clear();

				TIterator* itAllParticles = currentInputArray->MakeIterator();

				// Loop through all particles sequentially
				while(Propagate(static_cast<Candidate*>(itAllParticles->Next()), 0));
				// Once Propagate finds a particle that bends, it will recursively propagate
				// the entire decay chain, skipping ahead in the event record (but without
				// affecting itInput). Hence, once Propagate finally returns, the itInput loop
				// is guarenteed to iterate over particles which have already been propagated.
				// This is why Propagate always checks to see if TrackLength < 0 first (as
				// it will only be < 0 when the particle has yet to be Propagted).

				delete itAllParticles;
			}

			// Now loop through the output arrays and copy their objects to the storage vector
			// Step1: keep a running list of Candidates which exist in an output list
			{
				std::vector<std::vector<Candidate*> >& thisEntryStore = pileupStore[entry];

				Candidate* candidate = 0;

				for(size_t iArray = 0; iArray < pileupArraysToStore.size(); ++iArray)
				{
					std::vector<Candidate*>& thisStoreArray = thisEntryStore[iArray];
					TObjArray const* const thisAPProArray = pileupArraysToStore[iArray];

					thisStoreArray.reserve(thisAPProArray->GetEntriesFast());
					TIterator* itCandidate = thisAPProArray->MakeIterator();

					while((candidate = static_cast<Candidate*>(itCandidate->Next())))
					{
						// Map from the temporary to permanent Candidate
						Candidate*& permanentCandidate = permanentPileup[candidate];

						// Check to see if we've already instantiated a permanent version of this particle
						if(permanentCandidate == 0)
						{
							permanentCandidate = static_cast<Candidate*>(pileupFactory->NewEntry());
							candidate->Copy(*permanentCandidate);

							// Clear the Mother/Daughter indices, since they are intrinsically unreferenceable
							permanentCandidate->M1 = 0;
							permanentCandidate->M2 = 0;
							permanentCandidate->D1 = 0;
							permanentCandidate->D2 = 0;

							permanentCandidate->SetFactory(factory); // This is done so that the official Delphes factory is in charge of Candidate::Clone() (e.g. MomentumSmearing will Clone pileup)
							TProcessID::AssignID(permanentCandidate);
						}

						thisStoreArray.push_back(permanentCandidate);
					}

					delete itCandidate;
				}
			}
		}
		cout << "\n"; // End the progress bar

		delete currentInputArray;
		currentInputArray = 0;
		delete pileupReader;
		delete pileupFile;
	}

	fOutputArray->Clear();
	fChargedHadronOutputArray->Clear();
	fElectronOutputArray->Clear();
	fMuonOutputArray->Clear();
}

//------------------------------------------------------------------------------

void AllParticlePropagator::FillPileup()
{
	std::vector<UInt_t> pileupIndices;

	// Fill pileupIndices
	{
		UInt_t numPileup = gRandom->Poisson(fMeanPileup);

		const UInt_t sizePileupStore = pileupStore.size();
		while(numPileup > sizePileupStore)
		{
			// We already checked that the pileupStore was more than 10 times as large as the mean
			// In the rare case we get an insanely large number, simply re-draw.
			numPileup = gRandom->Poisson(fMeanPileup);
		}

		// Fill pileup WITHOUT replacement (by recording indices already selected)
		// std::vector<bool> is probably not faster in this case, since it's too space efficient.
		// To use std::vector<char>, 1 == free and 0 == used
		std::vector<char> freeIndex(sizePileupStore, 1);
		pileupIndices.reserve(numPileup);

		while(pileupIndices.size() < numPileup)
		{
			UInt_t possibleIndex = gRandom->Integer(sizePileupStore);

			char& isFree = freeIndex[possibleIndex];
			if(isFree == 1)
			{
				isFree = 0;
				pileupIndices.push_back(possibleIndex);
			}
		}
	}

	for(std::vector<UInt_t>::iterator itIndex = pileupIndices.begin(); itIndex not_eq pileupIndices.end(); ++itIndex)
	{
		std::vector<std::vector<Candidate*> >& thisEntryStore = pileupStore[*itIndex];

		// Currently, we do not include the entire event record -- only the particles are tracked, or hit the edge
		// Thus, it is not meaningful (and potentially dangerous) to use the Mother/Daughter indices of pileup to trace the lineage
		// If we were to include the entire event record (at a later date) we would need to re-index the pileup EACH time it is added
		// to "Delphes/allParticleArray" (reading the mother/daughter index of the system line (At(0)) to figure out what the LAST shift was,
		// then finding the difference between the new shift to re-shift).
		for(size_t iArray = 0; iArray < pileupArraysToStore.size(); ++iArray)
		{
			std::vector<Candidate*>& thisStoreArray = thisEntryStore[iArray];
			TObjArray* const thisAPProArray = pileupArraysToStore[iArray];

			for(std::vector<Candidate*>::iterator itCandidate = thisStoreArray.begin(); itCandidate not_eq thisStoreArray.end(); ++itCandidate)
			{
				//cout << (*itCandidate)->IsPU << "\t";
				thisAPProArray->Add(*itCandidate);
			}
		}
	}
	//cout << endl << "exiting loop" << endl;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// VecXY

VecXY::VecXY():
	x(0.), y(0.) {}

//------------------------------------------------------------------------------

VecXY::VecXY(int) {/*Don't initialize my vector, I know what I'm doing.*/}

//------------------------------------------------------------------------------

VecXY::VecXY(const Double_t x_in, const Double_t y_in):
	 x(x_in), y(y_in) {}

//------------------------------------------------------------------------------

VecXY& VecXY::operator += (const VecXY& otherVec)
{
	x += otherVec.x;
	y += otherVec.y;
	return *this;
}

//------------------------------------------------------------------------------

VecXY& VecXY::operator -= (const VecXY& otherVec)
{
	x -= otherVec.x;
	y -= otherVec.y;
	return *this;
}

//------------------------------------------------------------------------------

VecXY& VecXY::operator *= (const Double_t scalar)
{
	x *= scalar;
	y *= scalar;
	return *this;
}

//------------------------------------------------------------------------------

VecXY& VecXY::operator /= (const Double_t scalar)
{
	x /= scalar;
	y /= scalar;
	return *this;
}

//------------------------------------------------------------------------------

VecXY& VecXY::operator ~ ()
{
	x = -x;
	y = -y;
	return *this;
}

//------------------------------------------------------------------------------

VecXY VecXY::operator + (const VecXY& otherVec) const
{
	VecXY newVec(*this);
	newVec += otherVec;
	return newVec;
}

//------------------------------------------------------------------------------

/*
VecXY VecXY::operator + (VecXY&& otherVec) const
{
	otherVec += *this;
	return otherVec;
}
*/

//------------------------------------------------------------------------------

VecXY VecXY::operator - (const VecXY& otherVec) const
{
	VecXY newVec(*this);
	newVec -= otherVec;
	return newVec;
}

VecXY VecXY::operator - () const
{
	return VecXY(-x, -y);
}

//------------------------------------------------------------------------------

/*
VecXY VecXY::operator - (VecXY&& otherVec) const
{
	(~otherVec) += *this;
	return otherVec;
}
*/

//------------------------------------------------------------------------------

VecXY VecXY::operator * (const Double_t scalar) const
{
	VecXY newVec(*this);
	newVec *= scalar;
	return newVec;
}

//------------------------------------------------------------------------------

VecXY VecXY::operator / (const Double_t scalar) const
{
	VecXY newVec(*this);
	newVec /= scalar;
	return newVec;
}

//------------------------------------------------------------------------------

Double_t VecXY::Norm() const
{
	// The naive form is more accurate unless we're dealing with x or y so
	// large that the result will under/overflow. In this application ... unlikely.
	return sqrt(x*x + y*y);
	//return hypot(x, y);
}

//------------------------------------------------------------------------------

Double_t VecXY::Norm2() const
{
	return x*x + y*y;
}

//------------------------------------------------------------------------------

Double_t VecXY::Dot(const VecXY& otherVec) const
{
	return (x*otherVec.x + y*otherVec.y);
}

//------------------------------------------------------------------------------

Double_t VecXY::Cross(const VecXY& otherVec) const
{
	return (x*otherVec.y - y*otherVec.x);
}

//------------------------------------------------------------------------------

// thisVec> x z^
VecXY VecXY::CrossZHat() const
{
	return VecXY(y, -x);
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// RotationXY

RotationXY::RotationXY(const Double_t angle):
	cosine(cos(angle)), sine(sin(angle)) {}

//------------------------------------------------------------------------------

RotationXY::RotationXY(const Double_t cos_in, const Double_t sin_in):
	cosine(cos_in), sine(sin_in) {} // No check for unitarity, trust the user

//------------------------------------------------------------------------------

RotationXY::RotationXY(const VecXY& vec)
{
	const Double_t norm = vec.Norm();
	cosine = vec.x/norm;
	sine = vec.y/norm;
}

//------------------------------------------------------------------------------

RotationXY& RotationXY::Add(RotationXY const* const otherRotation)
{
	if(otherRotation)
	{
		{
			const Double_t
				cosine0 = cosine,
				sine0 = sine;

			// Add Angles
			// The naive form of angular addition for cosine is ripe for catastrophic cancellation
			// when the new cosine term is approximately zero.
			/// cosine = cosine0 * otherRotation->cosine  - sine0 * otherRotation->sine;

			// This alternate form gives the same result with only benign cancellation
			// (and only 4 additional additions and a final exponent shift)
			cosine = 0.5*((cosine0 + sine0)*(otherRotation->cosine - otherRotation->sine) +
							  (cosine0 - sine0)*(otherRotation->cosine + otherRotation->sine));
			sine = cosine0 * otherRotation->sine + sine0 * otherRotation->cosine;
		}

		// Ensure unitarity (multiplying by reciprocal would be faster, but less accurate).
		const Double_t normalization = sqrt(cosine*cosine + sine*sine);

		// However, we should be seeing errors of O(epsilon) (i.e. 2e-16), so check for absurd errors
		if(abs(normalization-1.) > 1e-14)
		{
			printf("\n\nsin: %.16e\ncos: %.16e\n(norm-1): %.16e\n\n", sine, cosine, normalization-1.);
			throw runtime_error("(RotationXY::Add): normalization too high, numerical error!");
		}

		cosine /= normalization;
		sine /= normalization;
	}

	return *this;
}

//------------------------------------------------------------------------------

// Rotate the candidate's 4-momentum in the XY plane
void RotationXY::ForwardRotateDaughterMomentum(Candidate* const daughter) const
{
	TLorentzVector& pMu = daughter->Momentum;
	// Copies of the original, since we're changing it
	Double_t
		px0 = pMu.Px(),
		py0 = pMu.Py();

	// Apply an ACTIVE rotation to the 4-momentum in the XY plane
	pMu.SetPx(  cosine * px0  - sine   * py0 );
	pMu.SetPy(  sine   * px0  + cosine * py0 );
}

//------------------------------------------------------------------------------

VecXY& RotationXY::ForwardRotate(VecXY& vec) const
{
	// Copies of the original, since we're changing it
	const Double_t
		x0 = vec.x,
		y0 = vec.y;

	// Apply an ACTIVE rotation to the vector
	vec.x =   cosine * x0  - sine   * y0;
	vec.y =   sine   * x0  + cosine * y0;

	return vec;
}

//------------------------------------------------------------------------------

VecXY RotationXY::ForwardRotate(VecXY&& vec) const
{
	// Copies of the original, since we're changing it
	const Double_t
		x0 = vec.x,
		y0 = vec.y;

	// Apply an ACTIVE rotation to the vector
	vec.x =   cosine * x0  - sine   * y0;
	vec.y =   sine   * x0  + cosine * y0;

	return vec;
}

//------------------------------------------------------------------------------

VecXY& RotationXY::ReverseRotate(VecXY& vec) const
{
	// Copies of the original, since we're changing it
	const Double_t
		x0 = vec.x,
		y0 = vec.y;

	// Reverse an ACTIVE rotation to the vector
	vec.x =   cosine * x0  + sine   * y0;
	vec.y = - sine   * x0  + cosine * y0;

	return vec;
}

//------------------------------------------------------------------------------

VecXY RotationXY::ReverseRotate(VecXY&& vec) const
{
	// Copies of the original, since we're changing it
	const Double_t
		x0 = vec.x,
		y0 = vec.y;

	// Reverse an ACTIVE rotation to the vector
	vec.x =   cosine * x0  + sine   * y0;
	vec.y = - sine   * x0  + cosine * y0;

	return vec;
}

//------------------------------------------------------------------------------

VecXY RotationXY::ForwardRotateCopy(const VecXY& vec) const
{
	// Apply an ACTIVE rotation to the vector
	return VecXY(
		 cosine * vec.x  - sine   * vec.y,
		 sine   * vec.x  + cosine * vec.y);
}

//------------------------------------------------------------------------------

VecXY RotationXY::ReverseRotateCopy(const VecXY& vec) const
{
	// Reverse an ACTIVE rotation to the vector
	return VecXY(
		  cosine * vec.x + sine   * vec.y,
		- sine   * vec.x + cosine * vec.y);
}

//------------------------------------------------------------------------------

Double_t RotationXY::CalculateAngle() const
{
	return atan2(sine, cosine);
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// DecayChainExtractor, accessible through its member function AddFullDecayChain(...),
// takes a PYTHIA event record and sucks out the full decay chains. So as not to be
// dependent on a specific implementation of the Candidate class, the function that
// initializes each Candidate is a callback function.
//
// Example usage:

/*
void MyPythiaParticleToDelphesCandidateFunction(Pythia8::Particle const&, Candidate* candidate)
{
	// Your initialization code
	// Setting the mother/daughter indices is pointless, as they will be subsequently overwritten
}

void SomeOtherFunction()
{
	// Some code
	DecayChainExtractor::AddFullDecayChain(pythia->event, factory, outputArray, &MyPythiaParticleToDelphesCandidateFunction);
}
*/

// Must declare static members
const Int_t DecayChainExtractor::PROCESSED_STATUS = 223;
Pythia8::Event* DecayChainExtractor::event = 0;
DelphesFactory* DecayChainExtractor::factory = 0;
TObjArray* DecayChainExtractor::outputArray = 0;
Int_t DecayChainExtractor::outputArrayCurrentSize = 0;
void (*DecayChainExtractor::PythiaParticleToDelphesCandidate)(Pythia8::Particle const&, Candidate* const) = 0;

void DecayChainExtractor::AddFullDecayChain(Pythia8::Event& event_in, DelphesFactory* const factory_in,
	TObjArray* const outputArray_in,
	void (* const PythiaParticleToDelphesCandidate_in)(Pythia8::Particle const&, Candidate* const))
{
	event = &event_in;
	factory = factory_in;
	outputArray = outputArray_in;
	PythiaParticleToDelphesCandidate = PythiaParticleToDelphesCandidate_in;

	outputArrayCurrentSize = outputArray->GetEntriesFast();

	// Loop through particles, skipping the system line (0)
	for(Int_t pythiaIndex = 1; pythiaIndex <= event->size(); ++pythiaIndex)
	{
		Pythia8::Particle& particle = event->operator[](pythiaIndex);

		if((particle.status() not_eq PROCESSED_STATUS) and (particle.tau() > 0.))
		{
			// A particle with a non-zero lifetime which hasn't been processed; the
			// start of a decay chain. Since this is the FIRST particle in the decay
			// chain that will be stored, it's mother won't be dereferencible.

			// Add a Candidate* to store this particle, which will therefore exist
			// at index (outputArrayCurrentSize)
			outputArray->Add(factory->NewCandidate());
			// Calculate the difference between Pythia's and Delphes' indexing for this particle
			const Int_t pythiaIndexShift = pythiaIndex - outputArrayCurrentSize;
			// Now add all its daughters (passing 0 as the last argument to
			// indicate that this particle's mother is not dereferencable).
			// Also increment the size of the outputArray b/c it just grew
			AddDaughters(outputArrayCurrentSize++, pythiaIndexShift, 0);
		}
	}
}

//------------------------------------------------------------------------------

void DecayChainExtractor::AddDaughters(const Int_t particleIndex, Int_t pythiaIndexShift, Int_t motherThenLastDaughterIndex)
{
	Int_t currentDaughterIndex;
	{
		Pythia8::Particle& particle = event->operator[](particleIndex + pythiaIndexShift);
		Candidate* const candidate = static_cast<Candidate*>(outputArray->At(particleIndex));

		PythiaParticleToDelphesCandidate(particle, candidate); // Initialize the candidate
		candidate->M1 = motherThenLastDaughterIndex; // Set the index of the mother
		candidate->M2 = 0; // A particle in a decay chain can't have two mothers

		particle.status(PROCESSED_STATUS); // Indicate the particle has been processed

		if(particle.daughter1() == 0)
		{
			// For absolute safety, because we can't control PythiaParticleToDelphesCandidate(...)
			candidate->D1 = 0;
			candidate->D2 = 0;
			return; // Particle has no daughters, nothing more to do
		}

		// We need to re-index this particle's daughter indices. It's first daughter
		// will be the next particle added to outputArray
		candidate->D1 = outputArrayCurrentSize;
		// The incoming (pythiaShift) applies to this particle and her siblings.
		// However, this particle's daughters can have a different pythiaShift (e.g
		// if her daughters aren't positioned directly after her in the event record).
		pythiaIndexShift = particle.daughter1() - candidate->D1;
		candidate->D2 = particle.daughter2() - pythiaIndexShift;
		// We already used (motherThenLastDaughterIndex) to set the mother index. Now
		// it becomes the index of the last daughter. This re-use keeps the recursive
		// memory space as compact as possible.
		motherThenLastDaughterIndex = candidate->D2;

		//	We will now add a Candidate* to the outputArray for every daughter this particle has.
		// It is important to do this now because, if the daughter itself has daughters, the daughter
		// will add them all recursively, and we don't want her putting them in the range [D1, D2].
		while(outputArrayCurrentSize <= motherThenLastDaughterIndex)
		{
			outputArray->Add(factory->NewCandidate());
			++outputArrayCurrentSize;
		}

		// Now set (currentDaughterIndex) to the first daughter for initializing loop
		currentDaughterIndex = candidate->D1;
	}

	// Add all the daughters (and they're daughters)
	for(; currentDaughterIndex <= motherThenLastDaughterIndex; currentDaughterIndex++)
	{
		AddDaughters(currentDaughterIndex, pythiaIndexShift, particleIndex);
	}
}

