/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef PropagatorAndPixelTracker_h
#define PropagatorAndPixelTracker_h

/** \class PropagatorAndPixelTracker
 *
 *  1.   This is a rewrite (by K. Pedersen) of AllParticlePropagator (by K. Pedersen),
 *  which was a rewrite of ParticlePropagator (by P. Demin). The motivation for
 *  AllParticlePropagator is explained more fully in its header file, but is
 *  also summarized here in #2.
 *
 *  2.   PropagatorAndPixelTracker (henceforth known as ProPxT) will propagate an *entire* 
 *  decay chain in the magnetic field. The propagation only begins once a decay chain 
 *  contains a particle with a non-zero actual lifetime (i.e. it skips any initial 
 *  QCD or QED decays). Once a "long-lived" particle is found, it is propagated, 
 *  and all of its descendents are specially handeled (i.e. even if they are neutral,
 *  they are rotated according to the cumulative rotation of their ancestors).
 *
 *  3.   The purpose of this module is to study the exact configuration of hits
 *  in the pixel detector during some events of interest. As such, ALL particles 
 *  need to be propagated helicly, even if the rotation is incredily tiny. This
 *  differs from AllParticlePropagator, where particles were only propagated helicly 
 *  if their expected rotation was above some minimum threshold (and below which, they 
 *  were propagated in a straight line).
 *
 *  4.   ProPxT does not currently simulate ANY energy loss (e.g. multiple scattering,
 *  bremsstrahlung, or cyclotron radiation) or other material effects (gamma -> e+e-).
 *  As such, there is no "range-out" for loopers. 
 *  A simple but efficient parameterization of these effects (especially range-out ) 
 *  would be a welcome addition.
 *
 *  5.   When simulating hits in the pixel detector, all particles are propagated
 *  EXACTLY (i.e. there is no eta cut controlling detection; a particle will be detected 
 *  in the pixel detector if its helix intersects the pixel detector). This means that 
 *  loopers can deposit multiple hits in pixel layers. However, since all particles are
 *  propagated IDEALLY (as mentioned in #3), loopers will not range-out. 
 *
 *  6.   Being detected in the pixel detector relies on several efficiencies:
 *    6a. When a particle physically strikes a pixel detector layer, a 4-vector describing
 *        its hit is only stored IFF it passes the layer's *intrinsic* efficiency :
 *          * Intrinsic efficiency = The probability that a particle will deposit enough charge to be detected.
 *          * There are seperate intrinsic efficiencies for muons and everything else (hadrons/electrons).
 *          * This does not include the probability that a hit won't be reconstructed due
 *            to cuts (e.g. bad cluster geometry), so the intrinsic efficiency should be
 *            independent of particle flux.
 *        There is a limitations to this efficiency model :
 *          * Since efficiency is calculated independently for each particle, a particle which
 *            deposits a below threshold charge cannot become a hit if the same pixel is struck
 *            by a second particle.
 *          * The intrinsic efficiency has no dependence on energy/eta
 *          * Barrel layers are modelled as infinitely thin; there is no accounting for grazing
 *            hits by loopers, which might streak the pixels. 
 *     6b. After all the hits in a layer have been collected, they will be binned in certain 
 *         regions of interest. Before this is accomplished, bins will be deactivated 
 *         based on the layers *online* efficiency (to simulate dead pixels). 
 *           @ Online efficiency = The probability that any given pixel is functioning properly. 
 *           @ There are seperate online efficiencies for pixel and endcap layers
 *         There are also limitations with this model:
 *           @ The location of dead pixels is randomly determined before each binning.
 *           @ Besides the seperate efficiency for barrel and endcap, no eta
 *             dependence is simulated.
 *
 *  7.  Because the entire decay chain is needed, including dereferencable daughter
 *  indices, PileupMerger cannot be used (in its 3.2.0 implementation). The inclusion
 *  of Pileup is natively supported by "add"-ing an InputArray containing pileup 
 *  in the configuration card (i.e. like a Merger module). Such arrays MUST
 *  include the entire decay chain (with functional/re-indexed daughter indices),
 *  though they don't need all the colored precursors. 
 *     This will likely require a new Pileup picking module (i.e. one that outputs 
 *  not only final state particles, but the entire, pertinent decay chain). This module
 *  is on an To-do list, but an example re-indexing function is supplied.
 *
 *  7.  ProPxT absolutely NEEDS to know each particles proper lifetime, so this field
 *  has been added to the Candidate. Two other fields have also been added, in the 
 *  hopes that more accurate parametric tracking efficiencies can be developed which
 *  rely on a particle's production radius and its track length (though this would
 *  also require the expansion of DelphesFormula):
 * 
 *  8.  New Candidate fields
 *      * Float_t  CTau               [in mm].
 *                              This is the particles proper lifetime (its randomly
 *                              chosen ACTUAL proper lifetime, not the AVERAGE 
 *                              proper lifetime for its species).
 *                              Pythia8 stores c*tau as Particle::Tau.
 *                              When a particle with descendents strikes the cylinder, 
 *                              terminating its decay chain, it's original lifetime
 *                              is replaced by its lifetime to the edge of the cylinder.
 *      * Float_t  CreationRadius     [in mm]. 
 *                              This is the distance, in the xy plane, from the beamline
 *                              to the particle's creation vertex. It can be used to 
 *                              determine the quality of a particle's seed (i.e. pixel seed?).
 *                              The original z-position can be accessed directly from
 *                              the 4-position of the creation vertex.
 *                              ATTENTION: This field is used as a flag for ProPxT, to 
 *                              determine when a Candidate has already been propagated. 
 *                              When it is negative, the particle has not been propagated
 *                              yet. Thus, it should be initialized to -1. 
 *      * Float_t  TrackLength       [in mm]. 
 *                              This is the total length of the particle's track as 
 *                              it moves throught the tracker. If is
 *                                   length = cTau * gamma * |p|/E = cTau * |p|/m = cTau*sqrt(gamma^2-1)
 *                              It can be viewed (roughly) as the number of track hits.
 *                              WARNING: No eta cut is applied by ProPxT; all propagating
 *                              particles will have a non-zero TrackLength, even if they 
 *                              are at eta = 10, and do not actually strike any sensors.
 * 
 *  9.  Instead of storing each Candidate's original, un-propagted Position by cloning the 
 *  the original Candidate and storing the clone in the fArray of the propagated Candidate,
 *  ProPxT stores the original position in the "Area" 4-vector (which is 
 *  unused by all Candidates except jets, which should not be propagated).
 *

 
 *
 *  \author K. Pedersen - Illinois Institute of Technology - https://github.com/keith-pedersen
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

// If the shower/hadronization/decay program does not include a magentic field,
// then placing a decay chain inside a magnetic field requires correcting the
// 4-momenta and creation vertices of all particles.

// When a charged particle bends in a magnetic field, it's stored 4-momentum 
// does not need to change (since its observed 4-momentum is determined from the 
// location of its track hits and the bending radius). The 4-momentum of its daughters,
// on the other hand, MUST change according to the angle of their mother's rotation.
// Furthermore, the daughters creation vertices are also wrong, and should be shifted. 

// The cumulative rotation and of a decay chain is kept in a DaughterRotation object.
// The position of each daughter is set to be identical to the propagated position of its mother.

// The DaughterShifter class rotates a particles 4-momentum and shifts its position. 
// It is preferrable to TLorentzRotation because we ONLY need XY rotations

class Candidate;
class TObjArray;
class TLorentzVector;
class VecXY;
class RotationXY;

class PropagatorAndPixelTracker: public DelphesModule
{
	public:
		PropagatorAndPixelTracker();
		~PropagatorAndPixelTracker();

		void Init();
		void Process();
		void Finish();

	private:
		Candidate* Rotate(Candidate* const daughter, RotationXY const* const rotation, const TLorentzVector& mothersPosition);
		void NullifyDaughters(Candidate* const candidate, const bool nullifyThisParticle = false);
		bool Propagate(Candidate* const candidate, RotationXY const* rotation);

		bool PropagateHelicly(Candidate* const candidate, bool decays,
		                      const VecXY& r0, const Double_t z0,
		                      const VecXY& r0Beta, const Double_t R0Beta, const Double_t z0Beta,
		                      const Double_t omegaOverC, 
		                      RotationXY const*& rotation,
		                      Double_t& ctProp, VecXY& rFinal);
		                      
		TObjArray* NewPixelArray(const char* name);
			                                       
		// Parameters
		Double_t fRadius, fRadius2;           //    assumed supplied in [meters], converted to [mm] by Init()
		                                      // Radius of full tracker cylinder (and squared radius, for re-use)
		
		Double_t fHalfLength;                 //    assumed supplied in [meters], converted to [mm] by Init()
		                                      // Half the length of the full tracker cylinder
		
		Double_t fBz;                         //    assumed supplied in [Tesla]
		                                      // Strengh of perflectly homogenous magnetic field
		                                      // alligned with the +z direction
		                                
		Double_t fMinTrackLength;             //    assumed supplied in [meters], converted to [mm] by Init()
		                                      // Used for a quick pre-sort of tracks which are entirely
		                                      // too small, so that the charged OutputArrays aren't
		                                      // filled with tiny, invisible tracks.
			
		TFolder* fPixelFolder;  // A TFolder for storing the TObjArrays of pixel objects
		ExRootTreeBranch* fPixelArrays; // An factory for creating TObjArrays of pixel objects
		
		std::vector<TObjArray*> fInputList; //! use "add InputArray @$%!" in configuration card
		TObjArray* currentInputArray; //!
		
		// Output arrays
		TObjArray* fOutputArray; //! Particles that intersect the cylinder
		TObjArray* fChargedHadronOutputArray; //! Charged hadron track Candidates
		TObjArray* fElectronOutputArray; //! Electron track Candidates
		TObjArray* fMuonOutputArray; //! Muon track candidates
		
		TObjArray* fBarrelLayers; //! A list of PixelBarrel objects
		TIterator* fItBarrel; 
		
		TObjArray* fEndcapLayers; //! A list of PixelEndcap objects
		TIterator* fItEndcap;
		
		ClassDef(PropagatorAndPixelTracker, 1)
};


// VecXY is a minimal class for working in the transverse plane
// 
class VecXY // Double checked 1/31/15
{
	public: 
		Double_t x, y;
		
		VecXY();    // Initialize to [0,0]
		VecXY(int); // Create uninitialized vector
		VecXY(const Double_t x_in, const Double_t y_in);
		
		VecXY&   operator += (const VecXY& otherVec);
		VecXY&   operator -= (const VecXY& otherVec);
		VecXY&   operator *= (const Double_t scalar);
		VecXY&   operator /= (const Double_t scalar);
		
		VecXY&   operator ~ (); // Parity flip (reverse coordinates)
		
		VecXY    operator + (const VecXY& otherVec) const;
		VecXY    operator + (VecXY&& otherVec)      const;
		VecXY    operator - (const VecXY& otherVec) const;
		VecXY    operator - (VecXY&& otherVec)      const;
		VecXY    operator * (const Double_t scalar) const;
		VecXY    operator / (const Double_t scalar) const;
		
		Double_t Norm()                             const;
		Double_t Norm2()                            const; // Norm**2
		Double_t Dot(const VecXY& otherVec)         const;
		Double_t Cross(const VecXY& otherVec)       const; // this> x other>

		VecXY    CrossZHat()                        const; // this> x z^
};

// RotationXY is a minimal class for performing ACTIVE rotations in the transverse plane
// A PASSIVE rotation is accomplished via the ReverseRotate() function
//
class RotationXY // Double checked 1/31/15
{
	protected:
		Double_t cosine, sine;
		
	public:
		RotationXY(const Double_t angle);
		RotationXY(const Double_t cos_in, const Double_t sin_in); // Does NOT check for unitarity
		RotationXY(const VecXY& vec);
		
		RotationXY&  Add(RotationXY const* const other);
		
		void         ForwardRotateDaughterMomentum(Candidate* const daughter) const;
		
		VecXY&       ForwardRotate(VecXY& vec)                                const;
		VecXY        ForwardRotate(VecXY&& vec)                               const;
	
		VecXY&       ReverseRotate(VecXY& vec)                                const;
		VecXY        ReverseRotate(VecXY&& vec)                               const;
		
		VecXY        ForwardRotateCopy(const VecXY& vec)                      const;
		VecXY        ReverseRotateCopy(const VecXY& vec)                      const;
		
		Double_t     CalculateAngle()                                         const;
};

// PixelHit stores the crucial information about a hit in the pixel detector: 
//        The particle, its instantaneous position, and instantaneous velocity
// The velocity can be combined with local information about the pixel detector
// to determine impact angle.
class PixelHit
{
	public: 
		Candidate* const particle;
		const VecXY r, 
			rBeta;
		const Double_t
			ct, z, 
			zBeta;
			
		PixelHit(Candidate* const particle_in,
		         const Double_t ct_in, const VecXY& r_in, const Double_t z_in, 
			      const VecXY& rBeta_in, const Double_t zBeta_in);	
};

// PixelBarrel stores the crucial information about a perfect, cylindrical pixel barrel layer:
//    radius, halfLength, and a vector of PixelHits
class PixelBarrel : public TObject
{
	private:
		std::vector<PixelHit> hits;
	
	public:
		const Double_t
			radius, radius2,
			halfLength,
			rPhiWidth, zWidth, thickness,
			lorentzAngle;			

		PixelBarrel(const Double_t radius_in, const Double_t halfLength_in,
			         const Double_t rPhiWidth_in, const Double_t zWidth_in, const Double_t thickness_in, 
			         const Double_t lorentzAngle_in);

		void AddHit(Candidate* const particle,
		            const Double_t ct, const VecXY& r, const Double_t z, 
						const VecXY& rBeta, const Double_t zBeta);
		void Clear();
		
		const std::vector<PixelHit>&  GetHits() const;
};

// PixelEndcap stores the crucial information about a pair of perfect disks, 
// orthogonal to the z-axis, which exist at (+/-)zPosition. PixelEndcap
// is not designed to be used to simulate the "blading" of a real detector;
// instead, all hits are placed at a single z-position. 
class PixelEndcap : public TObject
{
	private:
		std::vector<PixelHit> forward;
		std::vector<PixelHit> backward;
	
	public:		
		const Double_t
			innerRadius2,
			outerRadius2,
			zPosition,
			rPhiWidth, rWidth, thickness,
			lorentzAngle, bladeAngle;
		
		PixelEndcap(const Double_t innerRadius_in, const Double_t outerRadius_in, const Double_t zPosition_in,
			         const Double_t rPhiWidth_in, const Double_t rWidth_in, const Double_t thickness_in, 
			         const Double_t lorentzAngle_in, const Double_t bladeAngle_in);
		
		// Unlike PixelBarrel::AddHit, these don't take the z position, as it is known
		void AddForwardHit(Candidate* const particle,
		                   const Double_t ct, const VecXY& r, 
			                const VecXY& rBeta, const Double_t zBeta);
		void AddBackwardHit(Candidate* const particle,
		                    const Double_t ct, const VecXY& r,
			                 const VecXY& rBeta, const Double_t zBeta);			
		
		const std::vector<PixelHit>&  GetForwardHits() const;
		const std::vector<PixelHit>&  GetBackwardHits() const;
};

#endif
