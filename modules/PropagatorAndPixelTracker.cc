#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "TFolder.h"
//#include "TLorentzVector.h"
#include "Pythia8/Pythia.h"
#include <cmath>
#include "PropagatorAndPixelTracker.h"

const Double_t c_light = 299792458.; // [m/s]
const Double_t PI = acos(-1.);
const Double_t ABSURDLY_LARGE = pow(2., 71.);

const Double_t GeV_to_eV = 1.0E9;
const Double_t METERS_to_mm = 1E3;
//const Double_t mm_to_METERS = 1./METERS_to_mm;
const Double_t CM_to_mm = 1E1;
const Double_t microns_to_MM = 1E-3;
const Double_t DEG_to_rad = 2*PI/360.;

// This code eschews TMath, as its safety checks are not needed

PropagatorAndPixelTracker::PropagatorAndPixelTracker():
	fRadius(0.), fRadius2(0.), fHalfLength(0.), 	fBz(0.), fMinTrackLength(0.), 
	fPixelFolder(0), fPixelArrays(0), fInputList(), currentInputArray(0), 
	fOutputArray(0), fChargedHadronOutputArray(0), fElectronOutputArray(0), 
	fMuonOutputArray(0), fBarrelLayers(0), fItBarrel(0), fEndcapLayers(0), fItEndcap(0) {}
PropagatorAndPixelTracker::~PropagatorAndPixelTracker() {}

void PropagatorAndPixelTracker::Init()
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
	
	fMinTrackLength = GetDouble("MinTrackLengthRatio", .7) * fRadius; // Default = 70% of radius (~high purity tracks, minimal double-counting)
	
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
	
	// Pixel arrays
	{
		// Get pixel dimensions
		const Double_t 
			barrelRPhiWidth = GetDouble("BarrelRPhiWidth_microns", 100.) * microns_to_MM,
			barrelZWidth = GetDouble("BarrelZWidth_microns", 150.) * microns_to_MM,
			barrelThickness = GetDouble("BarrelThickness_microns", 274.) * microns_to_MM,
			lorentzAngle = GetDouble("LorentzAngle_deg", 22.) * DEG_to_rad,
		
			endcapRPhiWidth = GetDouble("EndcapRPhiWidth_microns", 100.) * microns_to_MM,
			endcapRWidth = GetDouble("EndcapRWidth_microns", 100.) * microns_to_MM,
			endcapThickness = GetDouble("EndcapThickness_microns", 274.) * microns_to_MM,
			endcapBladeAngle = GetDouble("EndcapBladeAngle_deg", 20.) * DEG_to_rad;
	
		// The pixel TObjArray(s) keep the barrel and endcap objects; these objects
		// internally store arrays of pixel hits. 
		// We want the hits to be cleared before each event, but not the TObjArray(s)
		// storing the barrel/endcap objects themselves.
		// This means that we'll have to create a new, special folder and branch
		// for storing the Pixel arrays, because the DelphesModule::ExportArray framework is 
		// designed for Candidate(s), and it clears all arrays before each event.
		fPixelFolder = NewFolder("Pixel");
		fPixelArrays = new ExRootTreeBranch("PixelArrays", TObjArray::Class(), 0);
		
		fBarrelLayers = NewPixelArray(GetString("PixelBarrelOutputArray", "barrelLayers"));
		fItBarrel = fBarrelLayers->MakeIterator();
		
		fEndcapLayers = NewPixelArray(GetString("PixelEndcapOutputArray", "endcapLayers"));
		fItEndcap = fEndcapLayers->MakeIterator();
		
		{
			ExRootConfParam pixelBarrelRadii = GetParam("PixelBarrelRadii_cm");
			ExRootConfParam pixelBarrelHalfLengths = GetParam("PixelBarrelHalfLengths_cm");
			
			const int numBarrelLayers = std::min(pixelBarrelRadii.GetSize(), pixelBarrelHalfLengths.GetSize());
	
			for(int layer = 0; layer < numBarrelLayers; ++layer)
			{
				const Double_t 
					barrelLayerRadius = pixelBarrelRadii[layer].GetDouble() * CM_to_mm,
					barrelLayerHalfLength = pixelBarrelHalfLengths[layer].GetDouble() * CM_to_mm;
					
				fBarrelLayers->Add(
					new PixelBarrel(barrelLayerRadius, barrelLayerHalfLength,
						             barrelRPhiWidth, barrelZWidth, barrelThickness,
						             lorentzAngle));
			}
		}
		
		{
			ExRootConfParam pixelEndcapInnerRadii = GetParam("PixelEndcapInnerRadii_cm");
			ExRootConfParam pixelEndcapOuterRadii = GetParam("PixelEndcapOuterRadii_cm");
			ExRootConfParam pixelEndcapZPositions = GetParam("PixelEndcapZPositions_cm");
			
			const int numEndcapLayers = std::min( std::min(pixelEndcapInnerRadii.GetSize(), pixelEndcapOuterRadii.GetSize()), 
			                                     pixelEndcapZPositions.GetSize());
	
			for(int layer = 0; layer < numEndcapLayers; ++layer)
			{
				const Double_t 
					endcapLayerInnerR = pixelEndcapInnerRadii[layer].GetDouble() * CM_to_mm,
				   endcapLayerOuterR = pixelEndcapOuterRadii[layer].GetDouble() * CM_to_mm,
				   endcapLayerZPosition = pixelEndcapZPositions[layer].GetDouble() * CM_to_mm;
					
				fEndcapLayers->Add( 
					new PixelEndcap(endcapLayerInnerR, endcapLayerOuterR, endcapLayerZPosition,
						             endcapRPhiWidth, endcapRWidth, endcapThickness,
						             lorentzAngle, endcapBladeAngle));
			}
		}
	}
}

TObjArray* PropagatorAndPixelTracker::NewPixelArray(const char* name)
{
	TObjArray* newArray = static_cast<TObjArray *>(fPixelArrays->NewEntry());
	newArray->Clear();
	newArray->SetName(name);
	fPixelFolder->Add(newArray);
	return newArray;
}

// Process() loops through all the Particles in the InputArrays, propagating
// any which have a non-zero lifetime. 
void PropagatorAndPixelTracker::Process()
{
	// Clear all pixel hits
	{
		{//Barrel
			fItBarrel->Reset();
			PixelBarrel* barrelLayer;
	
			while((barrelLayer = static_cast<PixelBarrel*>(fItBarrel->Next())))
				barrelLayer->Clear();
		}
		
		{//Endcap
			fItEndcap->Reset();
			PixelEndcap* endcapLayer;
	
			while((endcapLayer = static_cast<PixelEndcap*>(fItEndcap->Next())))
				endcapLayer->Clear();
		}
	}

	// Loop through all input arrays and propagate all particles
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

void PropagatorAndPixelTracker::Finish()
{
	// Delete all pixel objects
	{
		{//Barrel
			fItBarrel->Reset();
			PixelBarrel* barrelLayer;
	
			while((barrelLayer = static_cast<PixelBarrel*>(fItBarrel->Next())))
				delete barrelLayer;
		}
		
		{//Endcap
			fItEndcap->Reset();
			PixelEndcap* endcapLayer;
	
			while((endcapLayer = static_cast<PixelEndcap*>(fItEndcap->Next())))
				delete endcapLayer;
		}
	}
		
	delete fPixelArrays; // The branch will delete its objects
	delete fItBarrel;
	delete fItEndcap;
}

// This function rotates the daughter's 4-momentum and sets its position
Candidate* PropagatorAndPixelTracker::Rotate(Candidate* const daughter, RotationXY const* const rotation, const TLorentzVector& mothersPosition)
{
	if(daughter)
	{
		daughter->Position = mothersPosition; // The daughter's creation vertex is the mother's decay vertex
		rotation->ForwardRotateDaughterMomentum(daughter); // Rotate the daughter's 4-momentum
	}
	
	return daughter;
}

// Since PropagatorAndPixelTracker is intended to be used with longer decay chains than are 
// required to get to the cylinder, once a particle makes it to the cylinder's edge, 
// all of its decendents should be recursively flagged as null and propagated
void PropagatorAndPixelTracker::NullifyDaughters(Candidate* const mother, const bool nullifyThisParticle)
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
		throw runtime_error("(PropagatorAndPixelTracker::NullifyDaughters): Daughter index does not correspond to valid Candidate*! Did you re-index?");
		// For more information, see #6 in the class description in the header file
	}
}

// This function propagates a particle in the magentic field
// It returns the VALIDITY of the candidate pointer (NOT whether the candidate was propagated)
//
// IMPORTANT conventions:
//
// 1) There will be less math/conversion factors if we store momentum in its dimensionless form
//         beta  ==  v / c  ==  p(GeV) / energy(GeV)
//    with time as a distance (ct)
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
//             r>_creation     p>_T
//       f) Magnitude is denoted by losing the vector symbol and switching to a captial letter:
//             |r>| == R, |r>_creation|**2 == R_creation**2 
//       * Vector/cross product is lower-case ( x ), with a required space on either side):
//             L> == r> x p>    x^ x y^ = z^
//       * Scalar/dot product is period (.), no space required: 
//             |r>|**2 = R**2 = r>.r>
//       * Normal multiplication is still (*): 
//             x**2 = x*x
//       * The substrcipt T is used throughout to denote transverse (xy)
//           
bool PropagatorAndPixelTracker::Propagate(Candidate* const candidate,	RotationXY const* rotation) 
{
	if(not candidate) // test to make sure we can dereference the candidate
	{
		if(rotation)
		{
			// Propagate is only called by Propagate() or Process(). To get here with an 
			// invalid candidate pointer, but a valid rotation, means that the calling instance of 
			// Propagate() looked up a daughter index that didn't correspond to a valid Candidate
			throw runtime_error("(PropagatorAndPixelTracker::Propagate): Daughter index does not correspond to valid Candidate*! Did you re-index?");
			// For more information, see #6 in the class description in the header file
		}
		else return false;
	}
	
	// candidate->TrackLength < 0 is a flag that indicates a candidate was not yet processed
	if(candidate->TrackLength >= 0.)
	{
		return true; // If candidate was already processed, return immedietely
		// This is supposed to happen; see the long comment in Process().
	}
	
	// newRotation will be used as a flag to indicate that this function should "delete rotation" at the end
	bool newRotation = false;
	
	// Whenever there is an active rotation, all descendents of this particle will be
	// propagated recursively. This can create a lot of nested calls to Propagate. 
	// In order to keep memory space tidy, we'll keep all temporary/unneccessary 
	// variables in a scope that will close before the recursive call.
	
	// Begin temporary variable scope
	{
		// Assume the particle's proper lifetime is stored in [mm]. If the particle is
		// final state (Status == 1), make cTau absurdly large. This will help later,
		// when we find the shortest of two distances (decay or cylinder strike). It also
		// allows for the propagation of MASSLESS final state particles, whose cTau = 0.
		const Double_t cTau = ((candidate->Status == 1) ? ABSURDLY_LARGE : candidate->CTau);  //[mm]
		
		if(cTau <= 0.)
		{
			// The particle decayed instantaneously (or has an invalid cTau)
			candidate->TrackLength = 0.; // This is also a flag that says this particle has been processed
			
			const TLorentzVector& position = candidate->Position;
			candidate->Area = position; // It didn't go anywhere, but give it a valid initial position
			
			// Set CreationRadius, then check that the particle was created inside the cylinder
			if(((candidate->CreationRadius = hypot(position.X(), position.Y())) > fRadius)
				or (abs(position.Z()) > fHalfLength))
			{
				throw runtime_error("(PropagatorAndPixelTracker::Propagate): Particle created outside the cylinder! Did you include the entire decay chain? Did you re-index?");
				// For more information on this error, see #6 in the class description in the header file
			}
		}
		else
		{
			const TLorentzVector& momentum = candidate->Momentum; // 4-momentum
			TLorentzVector& position = candidate->Position; // creation vertex
	
			// Store the particle's creation vertex in its Area field (see #9 in class descrption in header file)
			(candidate->Area) = position;
		
			// We'll be using cylindrical coordinates (r, z) for obvious reasons
			// Assume 4-position is stored in [mm]
			// Assume 4-momentum & mass is stored in [GeV]
			// Assume charge is stored in elementary charged [-e]
			const VecXY r0(position.X(), position.Y()); 
			const Double_t
				energy = momentum.E(),
				z0 = position.Z(),
				q = candidate->Charge; 
				
			// Set CreationRadius, then check that the particle was created inside the cylinder
			if(((candidate->CreationRadius = r0.Norm()) > fRadius) or (abs(z0) > fHalfLength))
			{
				throw runtime_error("(PropagatorAndPixelTracker::Propagate): Particle created outside the cylinder! Did you include the entire decay chain? Did you re-index?");
				// For more information on this error, see #6 in the class description in the header file
			}
			
			// For double precision, machine epsilon ~= 2.2E-16. Thus, we don't
			// have to worry about these betas losing resolution until we're 
			// dealing with gamma > 10^7  (e.g. a {5 TeV e+} or a {1.4 PeV pi+})
			const VecXY r0Beta( momentum.Px()/energy, momentum.Py()/energy );
			const Double_t R0Beta = r0Beta.Norm(); // Transverse beta (used throughout instead of pT)
			const Double_t z0Beta = momentum.Pz()/energy;
			
			// The next 3 variables will be set by the propagating routine, when the solutions are found.
			
			// This bool lets us know Whether or not the particle decays. Currently, we will assume it does.
			bool decays = true; 
			
			// The total lab propagation time; currently, that's the time until decay (c * Tau * gamma).
			Double_t ctProp = cTau * (energy / candidate->Mass);
			// For massless particles, this will be (inf). That's OK, because we will repeatedly find
			// the shortest time. 
			
			// The particle's final position. Setting it to the original position allows the
			// final position to be calculated by += the displacement vector.
			VecXY rFinal(r0); 
			
			{// Figure out how long till the particle hits the endcap:
				const Double_t ctEndcap = (copysign(fHalfLength, z0Beta) - z0) / z0Beta; 
				// copysign() ensures that, if (z0Beta = (+/-)0.), we'll get (ctEndcap = inf)
				// (assuming the particle is inside the cylinder, which we already checked for)
								
				if(ctEndcap < ctProp)
				{
					ctProp = ctEndcap;
					decays = false;
				}
			}
			
			// Next, find any particle which CAN propagate helicly
			//    1. |q| != 0
			//    2. |R0Beta| != 0
			if((q * R0Beta) not_eq 0.) 
			{
				const Double_t omegaOverC = (-q * fBz / energy) * (c_light / (GeV_to_eV * METERS_to_mm)); // [rad/mm]
				
				newRotation =
				PropagateHelicly(candidate, decays,
					r0, z0, 
					r0Beta, R0Beta, z0Beta,
					omegaOverC,
					rotation, 
					ctProp, rFinal);
									
				// The only time we need a new rotation is when the particle decays 
				// (i.e. we'll have to propagate its daughters).
				decays = newRotation;
			}
			else // 6. Propagate particle in a straight line 
			{
				const VecXY r0BetaHat = r0Beta/R0Beta; //Normalize the transverse momentum
				const Double_t LzOverPt = r0.Cross(r0BetaHat);
				const Double_t ctBarrel = (-r0BetaHat.Dot(r0) + sqrt(fRadius2 - LzOverPt*LzOverPt))/R0Beta;
				
				if(ctBarrel < ctProp)
				{
					ctProp = ctBarrel;
					decays = false;
				}
				
				rFinal += r0Beta*ctProp;
			}
			
			// Set the final position
			position.SetXYZT(rFinal.x, rFinal.y,
				z0 + z0Beta * ctProp,
				position.T() + ctProp);
			
			// Set the track length and use it to do a simple pre-sort of track candidates,
			// so that the charged hardron output array doesn't get filled with a bunch of invisible tracks
			{
				const Double_t trackLength = ctProp * sqrt(R0Beta*R0Beta + z0Beta*z0Beta);
				
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
				candidate->Status = 1; // Ensure it is now considered final state
				fOutputArray->Add(candidate); // Only candidates which puncture the cylinder should hit the calorimter
				NullifyDaughters(candidate); // Terminate its decay chain by nullifying its daughters
				return true;
			}
		}//End cTau > 0.
	}//End temporary variable scope

	if(rotation) // If we have an active rotation, Propagate all daughters immedietely 
	{
		Int_t daughterIndex = candidate->D1;
	
		if(daughterIndex > 0) // There are daughters
		{
			while(daughterIndex <= candidate->D2)
			{
				// The pointer to the daughter is sent to Rotate, which then returns it, 
				// so it becomes the first argument of Propagate. 
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
	
		if(newRotation)
		{
			// The rotation pointed to by rotate was created for this particle
			// Now that all of its daughters have been propagated, we can delete it
			delete rotation;
		}
	}
	
	return true;
}

// This function has 4 returns (3 via pass-by-reference)
// The actual return indicates whether (rotation) was created anew.
// Ideally, the body of this function could be contained inside Propagate, but 
// keeping it as a seperate function keeps Propagate() much cleaner
bool PropagatorAndPixelTracker::PropagateHelicly(Candidate* const candidate, bool decays,
	const VecXY& r0, const Double_t z0,
	const VecXY& r0Beta, const Double_t R0Beta, const Double_t z0Beta,
	const Double_t omegaOverC, 
	RotationXY const*& rotation,
	Double_t& ctProp, VecXY& rFinal)
{
	const VecXY r0_helix = r0Beta.CrossZHat()/omegaOverC; // The vector from the center of the helix coordinate system to r0>
	const VecXY rCenter = (r0 - r0_helix); // The vector from the beam-origin to the center of the helix coordinate system
	
	const Double_t 
		R_helix = R0Beta/abs(omegaOverC),
		RCenter = rCenter.Norm();
		
	// For ease of working with angles, we will rotate to the normal, helixPrime coordinate system (hxPr).
	// The helix's close approach to the beamline lies at [-1, 0] and farthest withdrawl at [1, 0]. 
	// We don't know (chi), the angle required to rotate from the helix system to helixPrime, 
	// but we do know the cosine/sine of chi, so we don't really need the angle. 
	const RotationXY chi(rCenter.x/RCenter, rCenter.y/RCenter);
		
	// Now PASSIVELY rotate a normalized r0>_helix to the helixPrime frame (rho indicates unit vector)
	const VecXY rho0_hxPr = chi.ReverseRotate(r0_helix/R_helix); // Passive rotation is the same as reversing an active rotation
	
	const Double_t signQ = - copysign(1., omegaOverC);
	Double_t phi0; // We will set this based on information used to find the close approach
	
	// 4. Find the coordinates of close approach between the track circle and the beamline.
	{
		{
			Double_t RCloseApproach = RCenter - R_helix;
			candidate->Dxy = copysign(RCloseApproach, r0.Cross(r0Beta));
			// By convention, the impact parameter has the same sign as the angular momentum 
	
			RCloseApproach /= RCenter;
			candidate->Xd = rCenter.x * RCloseApproach;
			candidate->Yd = rCenter.y * RCloseApproach;
		}
		
		{
			// In the hxPr coordinate system, the angle of closest approach is at [-1, 0]
			// Thus, we can find the particle's initial displacement from close approach 
			// (deltaPhiCloseApproach = (phiCloseApproach - phi0)) by dotting rho0>_hxPr with 
			// [-1, 0], then giving it the same sign as y
			const Double_t deltaPhiCloseApproach = copysign(acos(-rho0_hxPr.x), rho0_hxPr.y);
			// Once again, copysign is safe because if y == (+/-)0, acos(-x) == 0.
	
			candidate->Zd = z0 + z0Beta * (deltaPhiCloseApproach / omegaOverC);
		
			phi0 = signQ*PI - deltaPhiCloseApproach; // Close approach is at sign(q) * Pi
		}
	}
						
	// cos(epsilon) = (R_barrel**2 - RCenter**2 - R_helix**2)/(2*R_helix*RCenter)
	const Double_t cosEpsilonDenominator = 2*R_helix*RCenter;
	const Double_t cosEpsilonRadiiTerm = -(RCenter*RCenter + R_helix*R_helix)/cosEpsilonDenominator;
	
	bool oneOutwardHit = false;
	RotationXY* thisRotation = 0;
							
	if(RCenter + R_helix >= fRadius) // The particle crosses the barrel
	{
		const Double_t cosEpsilonExit = fRadius2/cosEpsilonDenominator + cosEpsilonRadiiTerm;
		const Double_t phiBarrel = signQ*acos(cosEpsilonExit); // The particle is inside the barrel, only one solution
					
		const Double_t ctBarrel = (phiBarrel - phi0) / omegaOverC;
		
		if(abs(phi0) < PI)
			oneOutwardHit = true; // particle is past close approach
		
		if(ctBarrel < ctProp) // The barrel exit is the closest exit
		{
			if(ctBarrel < 0.)
				throw runtime_error("(PropagatorAndPixelTracker::Propagate): Keith, you fat fuck, your math is all wrong!");
			
			ctProp = ctBarrel;
			decays = false;
						
			// Set fFinal in the hxPr frame, scaled by R_helix
			rFinal.x = R_helix*cosEpsilonExit;
			rFinal.y = signQ*R_helix*sqrt(1.0 - cosEpsilonExit*cosEpsilonExit);
			
			// Move rFinal from the hxPr frame back to the helix frame, then add rCenter to take it back to the beam frame
			(chi.ForwardRotate(rFinal)) += rCenter;
			
			goto GeneratePixelHits;
		}
	}
	
	// Figure out the final rotation/position based on the propagation time
	thisRotation = new RotationXY(ctProp * omegaOverC);
	
	// Generate the rFinal
	rFinal = r0_helix;
	(thisRotation->ForwardRotate(rFinal)) += rCenter;
	
	GeneratePixelHits:
	// This was originally a seperate function, to keep the code sequestered.
	// However, it needed SO many arguments of common type that it was too risky (argument confusion)
	// A label, though less aesthetically pleasing, is the better solution. 
	{
		// Loop over barrel layers
		{
			fItBarrel->Reset();
			PixelBarrel* barrelLayer;
	
			while((barrelLayer = static_cast<PixelBarrel*>(fItBarrel->Next())))
			{
				// Check to make sure that particles with only one outward hit started inside the layer
				if(oneOutwardHit && (candidate->CreationRadius > barrelLayer->radius))
					continue;
				
				const Double_t zLastEdge = copysign(barrelLayer->halfLength, z0Beta);
				const Double_t ctLeaveBarrel = (zLastEdge - z0)/z0Beta;
				if(ctLeaveBarrel < 0.) // Check to make sure that particle isn't already past the edge of the barrel layer
					continue;
				
				const Double_t cosEpsilon = (barrelLayer->radius2)/cosEpsilonDenominator + cosEpsilonRadiiTerm;
				Double_t sinEpsilon = 1. - cosEpsilon*cosEpsilon;
				if(sinEpsilon < 0.) // Check to make sure that helix actually intersects the layer (i.e. |cosEpsilon| <= 1)
					continue; // This helps filter out the low-pt particles
							
				// We've passed important sanity checks, calculate epsilon and sin(epsilon)
				sinEpsilon = sqrt(sinEpsilon);
				const Double_t 
					epsilon = acos(cosEpsilon),
					ctMax = std::min(ctProp, ctLeaveBarrel); // How long till the last hit
					
				if(oneOutwardHit)
				{
					// Most particles only pass through one pixel layer; a quick treatment is faster
					const Double_t ctHit = (signQ*epsilon - phi0)/omegaOverC;
					
					if(ctHit > ctMax) // ctHit > 0 ensured with CreationRadius check
						continue;
					else if(ctHit < 0.)
						throw runtime_error("(PropagatorAndPixelTracker::Propagate): How'd we get here!");
					
					// Construct the hit in the helixPrime system, then rotate it back to helix system
					VecXY rHit = chi.ForwardRotate(VecXY(cosEpsilon * R_helix, signQ * sinEpsilon * R_helix));
					// This allows us to calculate the velocity vector
					const VecXY rBetaHit = rHit.CrossZHat() * omegaOverC;
					// Now we put the hit vector into the beam system
					rHit += rCenter;
						
					barrelLayer->AddHit(candidate, (candidate->Position).T() + ctHit, rHit, z0 + z0Beta*ctHit, rBetaHit, z0Beta);
				}
				else // Full treatment (loopers)
				{
					const Double_t ctPeriod = 2*PI*(-signQ/omegaOverC); // How long to make one loop
					Double_t	
						ctMin = (-zLastEdge - z0)/z0Beta, // calculate distance from other edge
						ctShift;
					
					if(ctMin > 0.) // See if the particle is outside the barrel, moving back into it
					{
						if(ctMin > ctProp)
							continue;
						// Propagate the particle an INTEGER number of periods, while still outside the barrel
						ctShift = floor(ctMin / ctPeriod) * ctPeriod;
					}
					else
					{
						ctMin = 0.;
						ctShift = 0.;
					}
				
					// The [0] hit is when the particle is moving towards close approach, [1] is when its moving away				
					const Double_t deltaCtHit[2] = {(signQ*(2*PI - epsilon) - phi0)/omegaOverC, (signQ*epsilon - phi0)/omegaOverC};
					const Double_t ySignHit[2] = {-signQ, signQ};

					for(unsigned int i = 0; i < 2; ++i) // Loop should unroll with optimization
					{
						Double_t ctHit = ctShift + deltaCtHit[i];
					
						// Ensure the first hits are on the barrel
						while(ctHit < ctMin) ctHit += ctPeriod;	
					
						// Check to make sure there are actually hits before constructing the hit vector
						if(ctHit <= ctMax)
						{
							// Construct the hit in the helixPrime system, then rotate it back to helix system
							VecXY rHit = chi.ForwardRotate(VecXY(cosEpsilon * R_helix, ySignHit[i] * sinEpsilon * R_helix));
							// This allows us to calculate the velocity vector
							const VecXY rBetaHit = rHit.CrossZHat() * omegaOverC;
							// Now we put the hit vector into the beam system
							rHit += rCenter;
						
							const Double_t ct0 = (candidate->Position).T();						
							for(; ctHit <= ctMax; ctHit += ctPeriod)
							{
								barrelLayer->AddHit(candidate, ct0 + ctHit, rHit, z0 + z0Beta*ctHit, rBetaHit, z0Beta);
							}
						}
					} 
				}// End full treatment
			}
		}
	
		// Loop over endcap layers
		{
			fItEndcap->Reset();
			PixelEndcap* endcapLayer;
		
			while((endcapLayer = static_cast<PixelEndcap*>(fItEndcap->Next())))
			{
				{
					// Check to make sure that the helix actually intersects the inner radius
					// of the layer. This helps filter out the low pt particles.
					const Double_t cosEpsilon = (endcapLayer->innerRadius2)/cosEpsilonDenominator + cosEpsilonRadiiTerm;
					if(abs(cosEpsilon) > 1.)
						continue;
				}
				
				// [0] is the forward hit, [1] is the backward hit
				const Double_t ctEndcap[2] = {(endcapLayer->zPosition - z0)/z0Beta, (-endcapLayer->zPosition - z0)/z0Beta};
				const bool forwardHit[2] = {true, false};
				
				// This loop checks for both the forward and backward hit.
				// It may seem unneccessary to check for both hits from the same particle, but it does happen occasionaly
				for(unsigned int i = 0; i < 2; ++i)
				{
					if(ctEndcap[i] < 0. or ctEndcap[i] > ctProp)
						continue; //The particle doesn't make it to the endcap
							
					// Propagate the particle to its strike position and calculate velocity
					VecXY rHit = RotationXY(ctEndcap[i] * omegaOverC).ForwardRotateCopy(r0_helix);
					const VecXY rBetaHit = rHit.CrossZHat() * omegaOverC;
					rHit += rCenter;
					
					const Double_t RHit2 = rHit.Norm2();
		
					// Make sure the hit is on the disk
					if((RHit2 < endcapLayer->innerRadius2) or (RHit2 > endcapLayer->outerRadius2))
						continue;
		
					if(forwardHit[i])
						endcapLayer->AddForwardHit(candidate, (candidate->Position).T() + ctEndcap[i], rHit, rBetaHit, z0Beta);
					else
						endcapLayer->AddBackwardHit(candidate, (candidate->Position).T() + ctEndcap[i], rHit, rBetaHit, z0Beta);
				}
			}
		}
	}
	
	if(decays)
	{
		thisRotation->Add(rotation); // Add this particle's rotation to that of its mother
		rotation = thisRotation; // Return (by reference) the new rotation for this particle's daughters
		return true; // Pass ownership of the rotation pointer to the calling instance of Propagate
	}
	else
	{
		delete thisRotation;
		return false;
	}
}

//PixelHit

PixelHit::PixelHit(Candidate* const particle_in, 
	                const Double_t ct_in, const VecXY& r_in, const Double_t z_in, 
	                const VecXY& rBeta_in, const Double_t zBeta_in):
particle(particle_in), r(r_in), rBeta(rBeta_in), ct(ct_in), z(z_in), zBeta(zBeta_in) {}

//PixelBarrel

PixelBarrel::PixelBarrel(const Double_t radius_in, const Double_t halfLength_in,
	                      const Double_t rPhiWidth_in, const Double_t zWidth_in, const Double_t thickness_in, 
	                      const Double_t lorentzAngle_in):
	TObject(), 
	hits(),
	radius(radius_in), radius2(radius*radius), halfLength(halfLength_in),
	rPhiWidth(rPhiWidth_in), zWidth(zWidth_in), thickness(thickness_in), lorentzAngle(lorentzAngle_in)
{}
		
void PixelBarrel::AddHit(Candidate* const particle,
		                   const Double_t ct, const VecXY& r, const Double_t z, 
	                      const VecXY& rBeta, const Double_t zBeta)
{
	hits.push_back(PixelHit(particle, ct, r, z, rBeta, zBeta));
	//hits.emplace_back(particle, ct, r, z, rBeta, zBeta); // C++11 
}

void PixelBarrel::Clear()
{
	hits.clear();
}

const std::vector<PixelHit>& PixelBarrel::GetHits() const
{
	return hits;
}




// PixelEndcap


PixelEndcap::PixelEndcap(const Double_t innerRadius_in, const Double_t outerRadius_in, const Double_t zPosition_in,
	                      const Double_t rPhiWidth_in, const Double_t rWidth_in, const Double_t thickness_in, 
	                      const Double_t lorentzAngle_in, const Double_t bladeAngle_in):
TObject(),
forward(), backward(),
innerRadius2(innerRadius_in*innerRadius_in), outerRadius2(outerRadius_in*outerRadius_in),
zPosition(zPosition_in),
rPhiWidth(rPhiWidth_in), rWidth(rWidth_in), thickness(thickness_in), 
lorentzAngle(lorentzAngle_in), bladeAngle(bladeAngle_in)
{}
		
void PixelEndcap::AddForwardHit(Candidate* const particle,
		             const Double_t ct, const VecXY& r,
	                const VecXY& rBeta, const Double_t zBeta)
{
	forward.push_back(PixelHit(particle, ct, r, zPosition, rBeta, zBeta));
	//forward.emplace_back(particle, ct, r, zPosition, rBeta, zBeta); // C++11 
}

void PixelEndcap::AddBackwardHit(Candidate* const particle,
		              const Double_t ct, const VecXY& r,
	                 const VecXY& rBeta, const Double_t zBeta)
{
	backward.push_back(PixelHit(particle, ct, r, -zPosition, rBeta, zBeta));
	//backward.emplace_back(particle, ct, r, -zPosition, rBeta, zBeta); // C++11 
}
	
const std::vector<PixelHit>&  PixelEndcap::GetForwardHits() const {return forward;}
const std::vector<PixelHit>&  PixelEndcap::GetBackwardHits() const {return backward;}




// VecXY

VecXY::VecXY(): 
	x(0.), y(0.) {}
	
VecXY::VecXY(int) {/*Don't initialize my vector, I know what I'm doing.*/}
	
VecXY::VecXY(const Double_t x_in, const Double_t y_in):
	 x(x_in), y(y_in) {}
		
VecXY& VecXY::operator += (const VecXY& otherVec)
{
	x += otherVec.x;
	y += otherVec.y;
	return *this;
}
		
VecXY& VecXY::operator -= (const VecXY& otherVec)
{
	x -= otherVec.x;
	y -= otherVec.y;
	return *this;
}
		
VecXY& VecXY::operator *= (const Double_t scalar)
{
	x *= scalar;
	y *= scalar;
	return *this;
}
		
VecXY& VecXY::operator /= (const Double_t scalar)
{
	x /= scalar;
	y /= scalar;
	return *this;
}

VecXY& VecXY::operator ~ ()
{
	x = -x;
	y = -y;
	return *this;
}	

VecXY VecXY::operator + (const VecXY& otherVec) const
{
	VecXY newVec(*this);
	newVec += otherVec;
	return newVec;
}

VecXY VecXY::operator + (VecXY&& otherVec) const
{
	otherVec += *this;
	return otherVec;
}

VecXY VecXY::operator - (const VecXY& otherVec) const
{
	VecXY newVec(*this);
	newVec -= otherVec;
	return newVec;
}

VecXY VecXY::operator - (VecXY&& otherVec) const
{
	(~otherVec) += *this;
	return otherVec;
}

VecXY VecXY::operator * (const Double_t scalar) const
{
	VecXY newVec(*this);
	newVec *= scalar;
	return newVec;
}

VecXY VecXY::operator / (const Double_t scalar) const
{
	VecXY newVec(*this);
	newVec /= scalar;
	return newVec;
}

Double_t VecXY::Norm() const
{
	return sqrt(x*x + y*y);
}

Double_t VecXY::Norm2() const
{
	return x*x + y*y;
}

Double_t VecXY::Dot(const VecXY& otherVec) const
{
	return (x*otherVec.x + y*otherVec.y);
}

Double_t VecXY::Cross(const VecXY& otherVec) const
{
	return (x*otherVec.y - y*otherVec.x);
}

// thisVec> x z^
VecXY VecXY::CrossZHat() const
{
	return VecXY(y, -x);
}


// RotationXY
RotationXY::RotationXY(const Double_t angle):
	cosine(cos(angle)), sine(sin(angle)) {}
	
RotationXY::RotationXY(const Double_t cos_in, const Double_t sin_in):
	cosine(cos_in), sine(sin_in) {} // No check for unitarity
	
RotationXY::RotationXY(const VecXY& vec)
{
	const Double_t norm = vec.Norm();
	cosine = vec.x/norm;
	sine = vec.y/norm;	
}
	
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
		
VecXY RotationXY::ForwardRotateCopy(const VecXY& vec) const
{
	// Apply an ACTIVE rotation to the vector
	return VecXY(
		 cosine * vec.x  - sine   * vec.y,
		 sine   * vec.x  + cosine * vec.y);		
}

VecXY RotationXY::ReverseRotateCopy(const VecXY& vec) const
{
	// Reverse an ACTIVE rotation to the vector
	return VecXY(
		  cosine * vec.x + sine   * vec.y,
		- sine   * vec.x + cosine * vec.y);
}

Double_t RotationXY::CalculateAngle() const
{
	return atan2(sine, cosine);
}
