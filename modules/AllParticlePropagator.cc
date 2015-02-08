#include "classes/DelphesClasses.h"
#include "TLorentzVector.h"
#include "Pythia8/Pythia.h"
#include <cmath>
#include "AllParticlePropagator.h"

// Create a brand new DaughterRotation
DaughterRotation::DaughterRotation(const Double_t theta_in):
	theta(theta_in), cosTheta(cos(theta)), sinTheta(sin(theta)) {}

// Create a DaughterRotation based on a previous one, with an additional rotation.
DaughterRotation::DaughterRotation(DaughterRotation const* const grandmother, const Double_t deltaTheta):
	theta( (grandmother ? (grandmother->theta + deltaTheta) : deltaTheta) ), cosTheta(cos(theta)), sinTheta(sin(theta)) {}

// Rotate the candidate's 4-momentum in the XY plane
void DaughterRotation::Rotate(Candidate* const candidate) const
{
	TLorentzVector& pMu = candidate->Momentum;
	Double_t pxOrig = pMu.Px();
	Double_t pyOrig = pMu.Py();
	
	// Apply an ACTIVE rotation to the 4-momentum in the XY plane
	pMu.SetPx(cosTheta * pxOrig - sinTheta * pyOrig);
	pMu.SetPy(sinTheta * pxOrig + cosTheta * pyOrig);	
}

const Double_t GeV_to_eV = 1.0E9;
const Double_t METERS_to_mm = 1E3;
//const Double_t mm_to_METERS = 1./METERS_to_mm;
const Double_t c_light = 299792458.; // [m/s]
const Double_t PI = acos(-1.);
const Double_t ABSURDLY_LARGE = pow(2., 71.);

// This code eschews TMath, as its safety checks are not needed

AllParticlePropagator::AllParticlePropagator() {}
AllParticlePropagator::~AllParticlePropagator() {}

void AllParticlePropagator::Init()
{
	// Find and import input arrays
	ExRootConfParam param = GetParam("InputArray");
	const Long_t size = param.GetSize();
	
	for(int i = 0; i < size; ++i)
	{
		fInputList.push_back(ImportArray(param[i].GetString()));
	}
	
	fRadius = GetDouble("Radius", 1.0) * METERS_to_mm; // Default = 1m
	fRadius2 = fRadius*fRadius;
	
	fHalfLength = GetDouble("HalfLength", 3.0) * METERS_to_mm; // Default = 3m
	
	fBz = GetDouble("Bz", 0.0); // Default = 0 Tesla (will cause fatal exception)
	
	fMinHelixAngle = GetDouble("MinHelixAngle", 5E-5); // Default ~= 1/360th of a degree, roughly the resolution of the CMS pixel tracker's middle layer
	
	fMinTrackLength = GetDouble("MinTrackLength", .7) * METERS_to_mm; // Default = 70 cm (high purity tracks, minimal double-counting)
	
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

void AllParticlePropagator::Finish() {/* No memory to clean up */}

// This function rotates the daughter's 4-momentum and sets its position
Candidate* AllParticlePropagator::Rotate(Candidate* const daughter, DaughterRotation const* const rotation, const TLorentzVector& mothersPosition)
{
	if(daughter)
	{
		daughter->Position = mothersPosition; // The daughter's creation vertex is the mother's decay vertex
		rotation->Rotate(daughter); // Rotate the daughter's 4-momentum
	}
	
	return daughter;
}

// Since AllParticlePropagator is intended to be used with longer decay chains than are 
// required to get to the cylinder, once a particle makes it to the cylinder's edge, 
// all of its decendents should be recursively flagged as null and propagated
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
bool AllParticlePropagator::Propagate(Candidate* const candidate,	DaughterRotation const* rotation) 
{
	if(not candidate) // test to make sure we can dereference the candidate
	{
		if(rotation)
		{
			// Propagate is only called by Propagate() or Process(). To get here with an 
			// invalid candidate pointer, but a valid rotation, means that the calling instance of 
			// Propagate() looked up a daughter index that didn't correspond to a valid Candidate
			throw runtime_error("(AllParticlePropagator::Propagate): Daughter index does not correspond to valid Candidate*! Did you re-index?");
			// For more information, see #6 in the class description in the header file
		}
		else return false;
	}
	
	// candidate->TrackLength < 0 is a flag that indicates a candidate was not yet propagated
	if(candidate->TrackLength >= 0.)
	{
		return true; // If candidate was already propagated, return immedietely
		// This is supposed to happen; see the long comment in Process().
	}
	
	// newRotation will be used as a flag to indicate that this function should "delete rotation" at the end
	bool newRotation = false;
	
	// Whenever there is an active rotation, all descendents of this particle will be
	// propagated recursively. This can create a lot of nested calls to Propagate. 
	// In order to keep memory space tidy, we'll keep all temporary/unneccessary 
	// in a scope that will close before the recursive call.
	
	// Begin temporary variable scope
	{
		// Assume the particle's proper lifetime is stored in [mm]. If the particle is
		// final state (Status == 1), make cTau absurdly large. This will help later,
		// when we find the shortest of two distances (decay or cylinder strike). It also
		// allows for the propagation of MASSLESS final state particles, whose cTau = 0.
		Double_t cTau = ((candidate->Status == 1) ? ABSURDLY_LARGE : candidate->CTau);  //[mm]
		
		if(cTau <= 0.)
		{
			// The particle decayed instantaneously (or has an invalid cTau)
			candidate->TrackLength = 0.; // Set the flag that says it's been processed
			
			const TLorentzVector& position = candidate->Position;
			candidate->Area = position; // It didn't go anywhere, but give it a valid initial position
			
			// Set CreationRadius, then check that the particle was created inside the cylinder
			if(((candidate->CreationRadius = hypot(position.X(), position.Y())) > fRadius)
				or (abs(position.Z()) > fHalfLength))
			{
				throw runtime_error("(AllParticlePropagator::Propagate): Particle created outside the cylinder! Did you include the entire decay chain? Did you re-index?");
				// For more information, see #6 in the class description in the header file
			}			
		}
		else
		{
			const TLorentzVector& momentum = candidate->Momentum; // 4-momentum
			TLorentzVector& position = candidate->Position; // creation vertex
	
			// Store the particle's creation vertex in its Area field (see #9 in class descrption in header file)
			(candidate->Area) = position;
		
			// Store local values of kinematic information, for simpler math expressions
			const Double_t
				// Assume 4-momentum & mass is stored in [GeV]
				energy = momentum.E(),
				
				// Assume 4-position is stored in [mm]
				x = position.X(),
				y = position.Y(),
				z = position.Z(),		
					
				// Assume charge is stored in elementary charged [|e|]
				q = candidate->Charge;
				
			// Set CreationRadius, then check that the particle was created inside the cylinder
			if(((candidate->CreationRadius = hypot(x, y)) > fRadius)
				or (abs(z) > fHalfLength))
			{
				throw runtime_error("(AllParticlePropagator::Propagate): Particle created outside the cylinder! Did you include the entire decay chain? Did you re-index?");
				// For more information, see #6 in the class description in the header file
			}
			
			const Double_t
				betaX = momentum.Px()/energy, // For double precision, machine epsilon ~= 2.2E-16. Thus, we don't 
				betaY = momentum.Py()/energy, // have to worry about these betas losing resolution until we're 
				betaZ = momentum.Pz()/energy; // dealing with gamma > 10^7  (e.g. a 5 TeV e+ or a 1.4 PeV Pi+)
								
			const Double_t
				Beta_T = hypot(betaX, betaY), // Transverse beta (used throughout instead of pT)
				LzOverPt = (x*betaY - y*betaX)/Beta_T, // The angular momentum (as a length) will come in handy several times
				ctDecay = cTau * energy / (candidate->Mass); // The total lab propagation time, in [mm], until decay
				// For massless particles, ctDeacy will be (inf). That's OK, because only the MINIMUM ct is used later on.
			
			// Figure out how long till the particle hits the endcap:
			const Double_t ctEndcap = (copysign(fHalfLength, betaZ) - z) / betaZ; // copysign(fHalfLength, betaZ)  ==  sign(betaZ)*fHalfLength
			// By using copysign, we can ensure that, if (betaZ == 0) or (betaZ == -0), we'll get (ctEndcap == inf, since we already ensured abs(z) < fHaflLength
			
			// The following variables will be set inside the propagation routines, 
			// depending on which route is taken.
			Double_t 
				x_final = x,
				y_final = y,
				ctProp = 0.;
			bool puncturedCylinder = false;
			
			// Next, find any particle which CAN propagate helicly
			//    1. |q| != 0
			//    2. |Beta_T| != 0
			if((q * Beta_T) not_eq 0.) 
			{
				// In order to model helical propagation, we'll need to know the gyration frequency:
				// 
				//     gyration frequency:    omega = -q fBz /(gamma m)     [rad/s]
				//
				//          Here, we use the omega sign convention of a RH xyz coordinate system, 
				//          where we are sitting at z^ and looking DOWN at the xy plane:
				//
				//              x^ = cos(0)   and    y^ = sin(Pi/2)   with     z^ || B>
				//          
				//          Thus, the LH rotation of a positively charged (q>0) particle is parameterized
				//          by a steadily  *decreasing*  angle (a negative frequency).
				//
				//     Since we'll be working with (ct) instead of (t), it will actually
				//     be easier to work with omega/c:
				//
				//          deltaPhi == omega*t == (omega/c)*ct
				//
				//          omega / c == -q fBz /(gamma m c) == -q fBz c / energy    [rad/m]
				//          omega / c == -q fBz c / energy / METERS_to_mm           [rad/mm]
				//				
				const Double_t omegaOverC = (-q * fBz / energy) * (c_light / (GeV_to_eV * METERS_to_mm)); // [rad/mm]
				
				// Next, we want to see if there is any point in undertaking the full helicle propagation.
				// For very short lived particles (or a weak magnetic field), the straight-line approximation is good enough
				
				if(ctDecay * omegaOverC > fMinHelixAngle)	// If the total bending angle is significant enough, do helical propagation
				{
					// 1. Now calculate how long to reach the barrel. We will use a helical coordinate system, where:
					// 
					//    r> == r>_center + r>helix ...................................................... from beamline to particle
					//
					//    r>_center == [x_center, y_center] = R_center * [cos(alpha), sin(alpha)] ........ from beamline to helix center
					//
					//    r>_helix  == R_helix * [cos(omega/c * ct + phi0), sin(omega/c * ct + phi0)] .... from helix center to particle
					//              == R_helix * [cos(phi(ct)), sin(phi(ct))]
									
					// 1a. First we need to find the center of the helix in the xy plane.
					//     We know that the particle's magnetic angular momentum is:
					// 
					//        l> == I * omega * z^ == (r>_helix x p>_T) == (r>_helix x beta>_T) * (energy/c)
					//     
					//     crossing in beta>_T from the left, and rearranging terms, we get: 
					// 
					//        r>_helix == (beta>_T x z^) / (omega / c)
					//
					//     Thus, to get from the particle's initial position to the center of the helix, we use -(r>_helix):
					const Double_t 
						x_center = x - betaY / omegaOverC, // [mm]
						y_center = y + betaX / omegaOverC; // [mm]
					
					// 1b. Now let's find 3 more important quantities:
					//        R_helix (the radius of the helix) .................. .......R_helix = (Beta_T * c) / abs(omega)  [mm]
					//        alpha (the angle from the origin to the center of the helix)
					//        R_center (the distance from the beam-line to the center of the helix)  [mm]
					const Double_t 
						R_helix = Beta_T / abs(omegaOverC), // [mm]
						alpha = atan2(y_center, x_center),
						R_center = hypot(x_center, y_center); // [mm]
						
					// 1c. The final piece we need is the particle's initial helix position (phi0).
					//     Calculating this requires some thought about the barrel exit solution.
					//     Even though the particle may not intersect the barrel (which we will check
					//     for later), we will lay out the barrel solution now.
					//
					//     Where does the particle cross the barrel (phiBarrel)? Well, if (r>_barrel) is the vector from
					//     the beamline to the particle's barrel intersect position, then the following equation must hold:
					//
					//        (r>_center + r>_helix)**2 == r>_barrel**2      
					//
					//     Expanding the dot product binomials (and using [cos(a)cos(b) + sin(a)sin(b) == cos(a-b)]) we get:
					// 
					//        epsilon == (phiBarrel - alpha) = acos((fRadius**2 - R_center**2 - R_helix**2) / (2 * R_helix * R_center))
					//
					//     Of course, since acos() only outputs the principle (or positive) angle, there are two solutions:
					//
					//        phiBarrel == alpha +/- epsilon
					//
					//     The question then becomes, which solution to use?
					//
					//     Well we obviously want the "closer" one, but this becomes tricky since "closer"
				   //     is defined in terms of the direction of propagation. However, think about a positively
					//     charged particle. It moves in a LH arc, constantly sweeping to smaller phi (omega < 0).
					//     Hence, the only way that   (phiBarrel = alpha - epsilon)  is the "closer" solution to phi0 
					//     is if phi0 starts BETWEEN the two solutions. But this would mean that the particle starts
					//     OUTSIDE the barrel (this is easiest to see via a sketch). Since the opposite argument
					//     works for the negatively charged particle, we find:
					// 
					//        phiBarrel == alpha + |q|*epsilon
					//  
					//     Thus, to find the time to intercept, we compute
					//
					//        ctBarrel = (phiBarrel - phi0)/omegeOverC
					//
					//     But in order for ctBarrel to be positive, we NEED:
					//
					//        sign(phiBarrel - phi0) == sign(omega) == -sign(q)
					//   
					//     Which (because phi0 can't start between the two barrel solutions) becomes:
					//     
					// 	    sign(alpha - phi0) == sign(alpha +/- epsilon - phi0).
					//
					//     Now that we know how to check that phi0 is in the correct place, we need to 
					//     calculate what it is. phi0 can be determined from a particles charge and initial
					//     momentum. The center of a positively charged particle's helix is always to its 
					//     right (i.e. -Pi/2 from it's direction of motion). Hence, from the perspective
					//     of the helix center, the positively charged particle's initial position is 90
					//     degress GREATER than it's direction of momentum (and vice versa for q < 0).
					//
					//     To correctly compare alpha and phi0, we also need phi0 to be in [-Pi, Pi].
					//     This is most easily accomplished by applying the active 90 deg rotation
					//     to p>_T *before* calling atan2.
					//
					//           [ cos(siqnQ*Pi/2)  -sin(signQ*Pi/2) ] 
					//  [px, py].[ sin(signQ*Pi/2)   cos(signQ*Pi/2) ] == signQ*(-py, px)
					
					const Double_t signQ = copysign(1., q); // It is probably OK to assume that q = +/-1,
					// since anything with q = +/-2 decays immedietely, but this is safer.
					// copysign() itself is safe since q == +/-0 won't make it past the first helix check
					
					Double_t phi0 = atan2(signQ*betaX, -signQ*betaY);
					
					// Vefify that phi0 is in the correct place 
					if(signQ*(alpha - phi0) > 0)
							phi0 += signQ*2*PI;
										
					// 2. Find the time to exit the cylinder
					Double_t ctCylinder;
					
					if(R_center + R_helix < fRadius)
						ctCylinder = ctEndcap; // The particle cannot cross the barrel
					else
					{
						const Double_t epsilon = acos( (fRadius2 - R_center*R_center - R_helix*R_helix) / (2 * R_helix * R_center) );
					
						const Double_t phi_exit = alpha + signQ * epsilon;
			
						const Double_t ctBarrel = (phi_exit - phi0)/omegaOverC;
						
						if(ctBarrel < 0.)
							throw runtime_error("(AllParticlePropagtor::Propagate): Keith, you fat fuck; your math is all wrong!");
			
						// Now we see which time is smaller, barrel or endcap
						ctCylinder = std::min(ctBarrel, ctEndcap);
					}
					
					// 4. Find the coordinates of close approach between the track circle and the beamline.
					
					// 4a. In the xy plane, this is found easily by drawing the helix circle about its center.
					//     The point of closest approach will lie on the line drawn through the origin and helix center,
					//     either in front of behind the beamline (depending on if R_helix > R_center)
					{
						Double_t closestApproach = R_center - R_helix;
						
						// The "impact parameter" is the total transverse close approach distance, 
						// with the convention that its sign is the same as the particle's angular momentum about the beamline
						candidate->Dxy = copysign(closestApproach, LzOverPt);
						
						// Scale the closest approach distance by R_center to use x_center and 
						// y_center to project out the close approach positions (avoiding trig functions)
						closestApproach /= R_center;
						candidate->Xd = x_center * closestApproach;
						candidate->Yd = y_center * closestApproach;
					}
						
					//  4c. To get the z coordinate when the particle is at close approach in the xy plane,
					//      add/subtract Pi from alpha (depending on the sign of the particle's omega) to 
					//      find the close approach angle (remember, we already forced phi0 to be in the
					//      correct place). This angle should usually be behind the particle, but can be
					//      in front of it if the particle is created moving backwards towards the beamline.
					{
						const Double_t phiClosestApproach = alpha + copysign(PI, signQ);
						const Double_t ctClosestApproach = (phiClosestApproach - phi0) / omegaOverC;
						
						candidate->Zd = z + betaZ * ctClosestApproach;
					}
		
					// 5. Now we can get the final propagated position
					
					// 5a. Did the particle decay before puncturing the cylinder?
					puncturedCylinder = ctCylinder < ctDecay;
					ctProp = puncturedCylinder ? ctCylinder : ctDecay;
							
					// 5b. Given ctProp, we can calculate the final position
					const Double_t phi_prop = omegaOverC * ctProp;
					x_final = (x_center + R_helix * cos(phi_prop + phi0)),
					y_final = (y_center + R_helix * sin(phi_prop + phi0));
					
					if(not puncturedCylinder)
					{
						// If the particle doesn't puncture the cylinder before it decays,
						// we'll need to rotate its daughters.
						newRotation = true;
						
						// Accumulate on top of existing rotation (if any)
						rotation = new DaughterRotation(rotation, phi_prop);
					}
				}// End helical is "useful"
			}// End helical is possible
			
			// 6. Propagate particle in a straight line 
			if(not (newRotation || puncturedCylinder)) // One of these will be set if the particle already propagated
			{
				// 6a. Once again, we'll use a transverse vector equation to solve for the propagation time.
				//     
				//         r>_barrel = r>_creation + beta>_T * ctBarrel
				//
				//     Because angular momentum is consvered:
				// 
				//         |r>_creation x beta^_T|**2 == |r>_barrel x beta^_T|**2
				//
				//     Which becomes:
				//
				//         LzOverPt**2 = fRadius**2 * (1 - (r^_barrel.beta^_T)**2)
				// 
				//     Which rearranges to:
				//     
				//         fRadius**2 - LzOverPt**2 = (r>_creation.beta^_T + ctBarrel * Beta_T)**2
				//
				//     So the final solution is: 
				//     
				//         ctBarrel = ( -r>_creation.beta^_T +/- sqrt(fRadius**2 - LzOverPt**2) )/Beta_T
				//
				//     Which we can nickname:
				// 
				//         ctBarrel = dotTerm +/- sqrtTerm
			
				// 6b. We don't need to worry about sqrtTerm being complex, since LzOverPt cannot be
				//     longer than the particle's creation vector (which we verified is inside the barrel).
				//     And, since the particle is in the barrel, there will only be one POSITIVE solution
				//     (the other must be the BACKWARD in time). Can we tell which is positive, without checking?
				// 
				//     Well, we know that sqrtTerm is manifestly positive, thus to prove that 
				//     sqrtTerm + dotTerm > 0, all we have to do is prove that |sqrtTerm| > |dotTerm|:
				//
				//        |sqrt(fRadius**2 - |LzOverPt|**2)| > |-(betaX*x + betaY*y)/Beta_T)|
				//
				//     which (upon squaring both sides) simplifies to:  
				//
				//        fRadius**2 > (x**2 + y**2)
				const Double_t dotTerm = -(betaX*x + betaY*y)/(Beta_T * Beta_T);
				const Double_t ctBarrel = dotTerm + sqrt(fRadius2 - LzOverPt*LzOverPt)/Beta_T;
				const Double_t ctCylinder = std::min(ctBarrel, ctEndcap);
				
				puncturedCylinder = ctCylinder < ctDecay;
				ctProp = puncturedCylinder ? ctCylinder : ctDecay;
				
				x_final = x + betaX * ctProp;
				y_final = y + betaY * ctProp;
				
				if(q not_eq 0.0)
				{
					// If fMinimumAngle is high enough, an extremely energetic, minimally-bending
					// particle could be tracked; we should set its coordinates of close approach.
					
					// l> == (r>_closesApproach x beta>_T). But since (r>_closesApproach.beta>_T == 0):
					candidate->Dxy = LzOverPt;
					
					// To project out the close approach components without resorting to a 
					// costly trig or sqrt function, let's use one last vector equation:
					// 
					//    r>_creation + ctCloseApproach * beta>_T == r>_closeApproach
					// 
					// Squaring both sides, then moving the beta>_T term to the RHS of the equation
					// and squaring again (relying once ag-ain on r>_closesApproach.beta>_T == 0),
					// gives a solution for ctCloseApproach.
					//
					//    ctCloseApproach = - r>_creation.beta>_T / Beta_T**2
					//
					// In fact, we've seen this term before. It's the "dotTerm" from the barrel solution.
					candidate->Xd = x + betaX * dotTerm;
					candidate->Yd = y + betaY * dotTerm;
					candidate->Zd = z + betaZ * dotTerm;					
				}
			}//End straight propagation
			
			// All propagation is now complete
			
			// Set the final position
			position.SetXYZT(x_final, y_final,
				z + betaZ * ctProp,
				position.T() + ctProp);
			
			// Set the track length and use it to do a simple pre-sort of track candidates,
			// so that the charged hardron output array doesn't get filled with a bunch of invisible tracks
			{
				const Double_t trackLength = ctProp * sqrt(Beta_T*Beta_T + betaZ*betaZ);
				
				candidate->TrackLength = trackLength;
				
				if(q not_eq 0.) // Charged particles, add to tracked particle output arrays
				{
					Int_t absID = abs(candidate->PID);
					
					if(absID == 13) // Always add muons regardless of track length
					{
						// Even if they were created in last mm of tracker, they still might show up in the muon chamber
						fMuonOutputArray->Add(candidate);
					}
					else if(trackLength >= fMinTrackLength) // Add other particles only with sufficient track length
					{
						if(absID == 11)
							fElectronOutputArray->Add(candidate);
						else // Don't add muons here because they are always added below
							fChargedHadronOutputArray->Add(candidate);
					}					
				}
			}
			
			if(puncturedCylinder)
			{
				// The candidate punctured the cylinder 
				candidate->Status = 1; // Make sure it is now considered final state
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

