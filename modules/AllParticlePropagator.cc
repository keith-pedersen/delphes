#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "TLorentzVector.h"
#include "Pythia8/Pythia.h"
#include <cmath>
#include "AllParticlePropagator.h"

#include <cstdio>

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
// so will the result.
Double_t KahanTriangleAreaPreSorted(const Double_t a, const Double_t b, const Double_t c)
{
	// EXTRA PARENTHESIS ARE DELIBERATE, DO NOT REMOVE
	return sqrt((a + (b + c))*(c - (a - b))*(c + (a - b))*(a + (b - c)))/4.;
	// See "Mathematics Written in Sand" by W. Kahan, bottom of page 10
	// http://www.cs.berkeley.edu/~wkahan/MathSand.pdf#page=10
}


// This code deliberately eschews TMath

AllParticlePropagator::AllParticlePropagator():
	fRadius(0.), fRadius2(0.), fHalfLength(0.), 	fBz(0.), fMinTrackLength(0.),
	fInputList(), currentInputArray(0), fOutputArray(0), fChargedHadronOutputArray(0), fElectronOutputArray(0),
	fMuonOutputArray(0) {}
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

	fBz = GetDouble("Bz_Tesla", 0.0); // Default = 0 Tesla (will cause fatal exception)

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
}

//------------------------------------------------------------------------------

void AllParticlePropagator::Finish() {/* No memory to clean up */}

//------------------------------------------------------------------------------

// This function rotates the daughter's 4-momentum and sets its position
Candidate* AllParticlePropagator::Rotate(Candidate* const daughter, RotationXY const* const rotation, const TLorentzVector& mothersPosition)
{
	if(daughter)
	{
		daughter->Position = mothersPosition; // The daughter's creation vertex is the mother's decay vertex
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
		if(nullifyThisParticle) // We don't want to alter the mother that initiates the recursion
		{
			mother->Status = 0; // This particle is null; the decay never happened
			mother->TrackLength = ABSURDLY_LARGE; // Make it appear that the particle was propagated.
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

	// Whenever there is an active rotation, all descendents of this particle will be
	// propagated recursively. This can create a lot of nested calls to Propagate.
	// In order to keep memory space tidy, we'll keep all temporary/unneccessary
	// variables in a scope that will close before the recursive call.

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

		if(cTau <= 0.) // The particle decayed instantaneously (or has an invalid cTau)
		{
			candidate->TrackLength = 0.; // This indicates that the particle has been "processed

			const TLorentzVector& position = candidate->Position;
			candidate->Area = position; // It didn't go anywhere, but give it a valid initial position

			const Double_t creationRadius = sqrt(position.X()*position.X() +  position.Y()*position.Y());

			// Check that the particle was created inside the cylinder (sanity check)
			if((creationRadius >= fRadius)	or (abs(position.Z()) >= fHalfLength))
			{
				printf("r0: %.16e\n", creationRadius);

				stringstream message;
				message << "(AllParticlePropagator::Propagate1): Particle created outside the cylinder ";
				message << "( " << creationRadius << " , " << position.Z() << " )!\n";
				message << "Did you include the entire decay chain? Did you re-index?\n";
				throw runtime_error(message.str());
				// For more information on this error, see #5 in the class description in the header file
			}

			candidate->CreationRadius = creationRadius;
		}
		else
		{
			const TLorentzVector& momentum = candidate->Momentum; // 4-momentum
			TLorentzVector& position = candidate->Position; // creation vertex

			// Store the particle's creation vertex in its Area field (see #8 in class descrption in header file)
			(candidate->Area) = position;

			// We'll be using cylindrical coordinates (r, z) for obvious reasons
			// Assume 4-position is stored in [mm]
			// Assume 4-momentum & mass is stored in [GeV]
			// Assume charge is stored in elementary charged [e]
			const VecXY r0(position.X(), position.Y());
			const Double_t
				energy = momentum.E(),
				R02 = r0.Norm2(),
				z0 = position.Z(),
				q = candidate->Charge;

			{
				const Double_t creationRadius = sqrt(R02);

				// Check that the particle was created inside the cylinder
				if((creationRadius >= fRadius) or (abs(z0) >= fHalfLength))
				{
					cout << "rotation: " << rotation << endl;
					printf("r0: %.16e\n", creationRadius);
                                        printf("r0 (no intermediate): %.16e\n", sqrt(position.X()*position.X() + position.Y()*position.Y()));
					cout << "M1: " << candidate->M1 << endl;
					cout << "M2: " << candidate->M2 << endl;
					
					Candidate* const mother = static_cast<Candidate*>(currentInputArray->At(candidate->M1));
					cout << "Mother PID: " << mother->PID << endl;
					printf("Mother CTau: %.16e\n" ,mother->CTau);
					cout << "Mother status: " << mother->Status << endl;
					
					const TLorentzVector& motherP4 = mother->Momentum;
					cout << "Mother pt: " << motherP4.Pt() << endl;

					const TLorentzVector& motherPos4 = mother->Position;
					
					printf("M r0: %.16e\n", sqrt(motherPos4.X()*motherPos4.X() + motherPos4.Y()*motherPos4.Y()));

					stringstream message;
					message << "(AllParticlePropagator::Propagate2): Particle (" << candidate->PID << ") created outside the cylinder ";
					message << "( " << creationRadius << " , " << z0 << " )!\n";
					message << "Did you include the entire decay chain? Did you re-index?\n";
					throw runtime_error(message.str());
					// For more information on this error, see #5 in the class description in the header file
				}

				candidate->CreationRadius = creationRadius;
			}

			// This maximum representable non-infinite gamma is gammaMax = 1/sqrt(MachineEpsilon).
			// For double precision, MachineEpsilon ~= 1.1E-16, so gammaMax = 9.5e7.
			// Thus, until we start dealing with gamma > 10^7, we don't have to
			// worry about beta losing precision (and also not unless the particle
			// is completely parallel to one of the axes). This works out to
			// a {5 TeV e+} or a {1.4 PeV pi+})
			const VecXY r0Beta( momentum.Px()/energy, momentum.Py()/energy );
			const Double_t R0Beta2 = r0Beta.Norm2(); // Transverse beta (used throughout instead of pT)
			const Double_t z0Beta = momentum.Pz()/energy;

			// The next 3 variables will be set by the propagating routine, when the solutions are found.

			// This bool lets us know whether or not the particle decays. Currently, we will assume it does.
			bool decays = true;

			// The total lab propagation time; currently, that's the time until decay (cTau * gamma).
			// Use the absolute value JUST IN CASE the mass (or CTau) is accidently -0.
			Double_t ctProp = abs(cTau * (energy / candidate->Mass));
			// For massless particles, this will be (inf). That's OK, because we will find the shortest time.

			// The particle's final position; default to initial position.
			VecXY rFinal(r0);

			{// Figure out how long till the particle hits the endcap:
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
				//    omega / c == -q fBz /(gamma m c) == -q fBz c / energy    [rad/m]
				//    omega / c == -q fBz c / energy / METERS_to_mm           [rad/mm]
				//
				const Double_t omegaOverC = (-q * fBz / energy) * (c_light / (GeV_to_eV * METERS_to_mm)); // [rad/mm]

				// PropagateHelicly() will propagate the particle. It has 4 returns. Its official
				// return indicates whether "rotation" (passed be reference) was changed. Of
				// course all charged particles rotate, but rotation only needs to change when
				// the particle decays. If it strikes the cylinder, its daughters doesn't require
				// rotation. The propagation solutions (ctProp and rFinal) are also passed by reference.
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

	if(rotation) // If we have an active rotation, Propagate all daughters immedietely
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
	// Hopefully the compiler can get these in the pipe, to reduce latency if they are needed.
	const bool notLooper = (RBeam_hx + R_hx >= fRadius);
	std::vector<Double_t> smallToLargeTriangleSides;

	// For ease of working with angles, let's rotate our perspective to the helixPrime
	// coordinate system (hxPr). In this system, the helix's perigree (close approach) to
	// the beamline lies at {R_helix, 0}, with apogee at {-R_helix, 0}. We can parameterize these positions
	// using (phi), with (phi == 0) corresponding to perigree and (phi == sign(omega)*pi) corresponding
	// to agpogee. This coordinate system (a 180 degree flip) keeps phi0 and
	// deltaPhiBarrel small for very enegetic particles, which need the most precision
	// (hence we don't need to do any angular subtraction, and subject our angle to epsilon error).

	Double_t phi0; // The initial anglular position (in hxPr) will be set in a moment
	const VecXY beta = rBeam_hx * R_hx; // This vector will be explained just before phi0 is set

	// Finish the notLooper calculations
	if(notLooper)
	{
		// Sort the sides of the triangle (to be explained later)
		smallToLargeTriangleSides = {fRadius, R_hx, RBeam_hx};
		std::sort(smallToLargeTriangleSides.begin(), smallToLargeTriangleSides.end());
	}

	{// Find the coordinates of close approach between the track circle and the beamline.
		{
			// In the xy plane, the point of close approach lies on the line drawn through
			// the origin and helix center, either in front of behind the beamline.
			Double_t RCloseApproach = R_hx - RBeam_hx;
			candidate->Dxy = copysign(RCloseApproach, r0.Cross(r0Beta));
			// By convention, the impact parameter has the same sign as the angular momentum

			// Project out the signed x and y components of the close approach position
			// No trig required, use the information we already have
			RCloseApproach /= RBeam_hx;
			candidate->Xd = rBeam_hx.x * RCloseApproach;
			candidate->Yd = rBeam_hx.y * RCloseApproach;
		}

		{
			// In the hxPr coordinate system, the angle of closest approach is at
			// rCloseApproach_hxPr> = {R_hx, 0}, which is parllel to rBeam_hx> in
			// the helix coordinate system.  Thus, the particle's initial angular position
			// (phi0) can be found from r0>_hx and rBeam>_hx.
			//
			// The most accurate way to find the angle between two vectors is NOT:
			//     acos(a>.b>/sqrt(a>.a> * b>.b>))
			// but (per W. Kahan):
			//
			//     2*atan(sqrt((alpha> - beta>)**2/(alpha> + beta>)**2)
			//
			// where (alpha> = Beta*alpha>) and (beta> = Alpha*beta>). This
			// form may seem like overkill, but the acos form loses precision
			// for exactly the particles we need extra precision for: super
			// energetic particles that don't bend very much, and thus start
			// very close to apogee. Such particles really need to show up
			// in the correct Calorimeter cell. Plus, the acos can sometimes
			// send a cos(x) > 1 for very parallel vectors (normalization error).
			// This producess a NaN from acos.
			//
			// We'll still need to find the sign of phi0 (which doesn't depend
			// on the charge ... some sketches may be required to understand why):
			//
			//    sign(phi0) = sign(rBeam>_hx x r0>_hx)

			const VecXY alpha = r0_hx * RBeam_hx;

			phi0 = copysign(2*atan(sqrt((alpha - beta).Norm2()/(alpha + beta).Norm2())),
					rBeam_hx.Cross(r0_hx));
			// Here copysign is safe because if rBeam_hx.Cross(r0_hx) = +/- 0,
			// then either the particle is at perigree or apogee. If it's at
			// perigree, then phi = +/-0 (and the sign of zero shouldn't matter).
			// If it's at apogee, then it's a looper, so phi_0 is only used for
			// Zd (which now is a meaningless number, since the particle is at
			// apogee, so close approach is equally in front and behind).

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
		// with the barrel is if phi0 starts BETWEEN the two solutions. But since
		// (RCenter + R_helix >= fRadius), this would also mean that the particle started
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
		const Double_t sinEpsilon =
			copysign(2. * KahanTriangleAreaPreSorted(smallToLargeTriangleSides[2], smallToLargeTriangleSides[1], smallToLargeTriangleSides[0]) / denom,
				omegaOverC);
		// Assume RBeam2_hx > R2_hx > fRadius2
		const Double_t cosEpsilon = (RBeam2_hx + (R2_hx - fRadius2))/(2.*denom);

		const Double_t ctBarrel = (atan2(sinEpsilon, cosEpsilon) - phi0) / omegaOverC;

		if(abs(cosEpsilon) > 1.)
			throw runtime_error("(AllParticlePropagator::PropagateHelicly): |cosEpsilon| > 1, numerical error!");
		if(abs(sinEpsilon) > 1.)
			throw runtime_error("(AllParticlePropagator::PropagateHelicly): |sinEpsilon| > 1, numerical error!");

		if(ctBarrel <= ctProp) // The barrel exit is the closest exit
		{
			ctProp = ctBarrel;

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

			// Because we know the cos and sin of epsilon, we know where the exit vector is in hxPr
			// Creating it in hxPr then rotating back to helix is much faster than using additional trig functions
			rFinal.x = R_hx*cosEpsilon;
			rFinal.y = R_hx*sinEpsilon;

			// In hxPr, rBeam>_hxPr is at phi = 0. Thus, to rotate our exit vector to
			// the helix coordinate system, we can use the coordinates of rBeam>_hx.
			// Then we simply have to subtract rBeam>_hx to get back to the beam frame
			(RotationXY(rBeam_hx.x/RBeam_hx, rBeam_hx.y/RBeam_hx).ForwardRotate(rFinal)) -= rBeam_hx;

			// The particle struck the barrel; no daughters, no new rotation
			return false;
		}
	}

	// The particle doesn't strike the barrel; to find it's exit position,
	// we simply rotate r0>_helix
	RotationXY* thisRotation = new RotationXY(ctProp * omegaOverC);
	rFinal = r0_hx;
	(thisRotation->ForwardRotate(rFinal)) -= rBeam_hx;

	if(decays)
	{
		// Create a new cumulative rotation for the daughters by adding the rotation of the mother
		thisRotation->Add(rotation);
		rotation = thisRotation; // Return (by reference) the new rotation for this particle's daughters
		return true; // Pass ownership of the (rotation) memory to the calling instance of Propagate
	}
	else
	{
		delete thisRotation;
		return false;
	}
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
	// The original is more accurate unless we're dealing with x or y so
	// large that the result will under/overflow. In this application, unlikely.
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
	cosine(cos_in), sine(sin_in) {} // No check for unitarity

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
		const Double_t
			cosine0 = cosine,
			sine0 = sine;

		// Add Angles
		cosine = cosine0 * otherRotation->cosine  - sine0 * otherRotation->sine;
		sine =   cosine0 * otherRotation->sine    + sine0 * otherRotation->cosine;
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

