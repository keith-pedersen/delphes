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

#ifndef AllParticlePropagator_h
#define AllParticlePropagator_h

/** \class AllParticlePropagator
 *
 *  1.   This is a re-write (by K. Pedersen) of ParticlePropagator (by P. Demin),
 *  with more complete propagtion in the magentic field and tracker material.
 *  Actual code from ParticlePropagator was not used, but its underlying math
 *  was used extensively in the crafting of this module.
 *
 *  2.   AllParticlePropagator (henceforth known as APPro) can propagate an *entire*
 *  decay chain in the magnetic field. The propagation only begins once a decay chain
 *  contains a particle with a non-zero actual lifetime (i.e. it skips any initial
 *  QCD or QED decays). Once a "long-lived" particle is found, it is propagated,
 *  and all of its descendents are specially handeled. Given the importance of
 *  ParticleFlow (i.e. eFlow) in event reconstruction, the author determined
 *  that such a module could be useful (especially since programs like Pythia generally
 *  do not simulate magnetic field effects).
 *
 *  3.  Imagine a 1 GeV K+ which decays early, after 20 cm, as (K+ -> mu+ nu_mu).
 *  The original ParticlePropagator will only propagate the muon in the magnetic field,
 *  as if the kaon was neutral. APPro makes both bend in the magnetic field (changing
 *  the creation vertex and 4-momentum of daughters based on the bending of the kaon).
 *  A cumalative bending (e.g. 3 deg left, 5 deg right, ...) is mantained as the
 *  decay chain propagates, along with a cumulative shifting of creation vertices.
 *      Such a decay is rare, but this module was created with the intention of
 *  studying common muons in rare b-jets vs. rare muons in overwhelming light jet background.
 *
 *  4.  APPro is intended to be used with the following Pythia options
 *        ParticleDecays:limitRadius = on
 *        ParticleDecays:rMax = 20000  ! 20m of decay (change as desired)
 *        211:mayDecay = true  !let charged pions decay
 *        321:mayDecay = true  !let charged kaons decay
 *        130:mayDecay = true  !let K_L0 decay
 *        13:mayDecay = false  !don't let muons decay (too rare at normal energies)
 *  These options turn on the decays of interest, then limit all decays
 *  to a sphere of 20m.
 *     Without decay limitations, (and with no "stable" particles anymore)
 *  EVERYTHING will decay to protons, neutrons, electrons, and neutrinos
 *  (cluttering up the event record). IMHO, ParticleDecays:limitRadius is preferrable to
 *  ParticleDecays:limitCylinder because the radius better describes how much
 *  matter the particle will traverse (i.e. low-eta loopers).
 *     Even with the decay limit, many decay chains will be significantly longer than
 *  they need to be. Once a particle is propagated to the edge of the cylinder,
 *  its decay chain is terminated (so any descendents it may have are thrown away,
 *  though its siblings are unaffeteced).
 *
 *  5.  Because the entire decay chain is needed, including dereferencable daughter
 *  indices, PileupMerger cannot be used (in its 3.2.0 implementation). The inclusion
 *  of Pileup is natively supported by "add"-ing an InputArray containing pileup
 *  in the configuration card (i.e. like a Merger module). Such arrays MUST
 *  include whole decay chains (with functional/re-indexed daughter indices),
 *  though they don't need all the colored precursors.
 *     This will likely require a new Pileup picking module (i.e. one that outputs
 *  not only final state particles, but the entire, pertinent decay chain). This module
 *  is on an To-do list, but a static class that strips whole decay chains from a
 *  Pythia event record is supplied.
 *
 *  6.  APPro absolutely NEEDS to know each particles proper lifetime, so this field (CTau)
 *  has been added to the Candidate. Two other fields have also been added, in the
 *  hopes that more accurate parametric tracking efficiencies can be developed which
 *  rely on a particle's CreationRadius and its TrackLength (though this would
 *  also require the expansion of DelphesFormula):
 *
 *  7.  New Candidate fields
 *      * Float_t  CTau               [in mm].
 *                              This is the particles proper lifetime (its randomly
 *                              chosen ACTUAL proper lifetime, not the 1 e-fold
 *                              average proper lifetime for its species).
 *                              Pythia8 stores c*tau as Particle::Tau.
 *      * Float_t  CreationRadius     [in mm].
 *                              This is the distance, in the xy plane, from the beamline
 *                              to the particle's creation vertex. It can be used to
 *                              determine the quality of a particle's track seed.
 *                              The original z-position can be accessed directly from
 *                              the 4-position of the creation vertex. *
 *      * Float_t  TrackLength       [in mm].
 *                              This is the total length of the particle's track as
 *                              it moves throught the tracker:
 *                                   TrackLength = ctPropagation * beta
 *                              It can be viewed (roughly) as the number of track hits,
 *                              and is used as a quick pre-sort to keep tiny, invisible
 *                              tracks from proliferating in the track output arrays.
 *                              WARNING: No eta cut is applied by APPro; all propagating
 *                              particles will have a non-zero TrackLength, even if they
 *                              are at eta = 10, and cannot actually strike any sensors.
 *                              ATTENTION: This field is used as a flag for APPro, to
 *                              determine when a Candidate has already been propagated.
 *                              When it is negative, the particle has not been proccessed
 *                              by APPro yet. Thus, it should be initialized to -1.
 *
 *  8.  Instead of storing each Candidate's original, un-propagted Position by cloning the
 *  original Candidate and storing the clone in the fArray of the propagated Candidate,
 *  APPro stores the original position in the "Area" 4-vector (which is
 *  unused by all Candidates except jets, which should not be propagated).
 *
 *  9. APPro does not currently simulate any energy loss from multiple scattering,
 *  bremsstrahlung, or cyclotron radiation. A simple parameterization
 *  of these effects would be a welcome addition.

 *
 *  \author K. Pedersen - Illinois Institute of Technology - https://github.com/keith-pedersen
 *
 */

#include "classes/DelphesModule.h"
#include <vector>

class Candidate;
class TObjArray;
class TIterator;
class TLorentzVector;
class VecXY;
class RotationXY;

class AllParticlePropagator: public DelphesModule
{
	public:
		AllParticlePropagator();
		~AllParticlePropagator();

		void Init();
		void Process();
		void Finish();

	private:
		Candidate* Rotate(Candidate* const daughter, RotationXY const* const rotation, const TLorentzVector& mothersPosition);
		// Rotate's a daughter's 4-momentum and sets its position to (mothersPosition)

		void NullifyDaughters(Candidate* const candidate, const bool nullifyThisParticle = false);
		// When a particle hits the cylinder, the rest of its decay chain should be nullified.
		// NullifyDaughters(...) is called recurisively on all daughters.

		bool Propagate(Candidate* const candidate, RotationXY const* rotation);
		// Propagates a particle to its final position (decay or cylinder strike)

		bool PropagateHelicly(Candidate* const candidate, const bool decays,
		                      const VecXY& r0, const Double_t z0,
		                      const VecXY& r0Beta, const Double_t z0Beta,
		                      const Double_t omegaOverC,
		                      RotationXY const*& rotation,
		                      Double_t& ctProp, VecXY& rFinal);
		// Charged particles require helicle propagation. To keep Propagate(...) clean,
		// helicle propagation is accomplished in a separate function.

		void PropagateAndStorePileup(const std::string& pileupFileName);
		void FillPileup();

		Double_t fRadius, fRadius2;     //    supplied in [meters], converted to [mm] by Init()
		                                // Radius of full tracker cylinder (and squared radius, for re-use)

		Double_t fHalfLength;           //    supplied in [meters], converted to [mm] by Init()
		                                // Half the length of the full tracker cylinder

		Double_t fBz;                   //    supplied in [Tesla]
		                                // Strengh of perflectly homogenous magnetic field (B> = fBz z^)

		Double_t fMinTrackLength;       //    supplied as ratio of fRadius, converted to [mm] by Init()
		                                // Used for a quick pre-sort of tracks which are entirely
		                                // too small to show up in charged OutputArrays

		Double_t fMeanPileup;              //    supplied as Integer
		                                // Average number of pileup events. Use to parameterize a
		                                // Poisson distribution that randomly chooses number of
		                                // pileup per event.

		// It would be nice to parameterize the energy loss, but this is not supported yet.

		// Input
		std::vector<TObjArray*> fInputList; //! use "add InputArray @$%!" in configuration card
		TObjArray* currentInputArray; //!

		// Pileup
		std::vector<std::vector<std::vector<Candidate*> > > pileupStore;
		std::vector<TObjArray*> pileupArraysToStore;
		ExRootTreeBranch* pileupFactory;

		// Output
		TObjArray* fOutputArray; //! Particles that intersect the cylinder
		TObjArray* fChargedHadronOutputArray; //! Charged hadron track Candidates
		TObjArray* fElectronOutputArray; //! Electron track Candidates
		TObjArray* fMuonOutputArray; //! Muon track candidates

		// Randomness test
		//
		// I noticed that there was SEVERE lack of uniformity in pt and eta distributions
		// for ALL jets (not for the top jets, which are still smooth). 
		// I attributed this to pileup, and to the small sample from which it is drawn.
		// But with a sample of 2000 and a <mu>=40, I couldn't account for the severity
		// of the spikes.
		// Then, I realized that whatever spikes exist in the pileup particles is ENHANCED
		// by the anti-kt clustering, so it looks much worse in jets. 
		// std::vector<UInt_t> selectionTracker;

		ClassDef(AllParticlePropagator, 2)
};

// If the shower/hadronization/decay program does not include a magentic field,
// then placing a decay chain inside a magnetic field requires correcting the
// 4-momenta and creation vertices of all particles.

// When a charged particle bends in a magnetic field, it's stored 4-momentum
// does not need to change (since its observed 4-momentum is determined from the
// location of its track hits and the bending radius). The 4-momentum of its DAUGHTERS,
// on the other hand, MUST change according to the angle of their mother's rotation.
// Furthermore, the daughters creation vertices are also wrong (i.e. because there
// mother was propagated in a straight line), and should also be shifted.

// In order to fascilliate the the cumulative rotation of a decay chain, two classes were
// created. These classes likely duplicated the efforts of other classes (e.g. Root), but
// they are more compact. For example, since we ONLY need rotations in the XY plane, using
// TLorentzRotation is overkill.

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
//		VecXY    operator + (VecXY&& otherVec)      const;
		VecXY    operator - (const VecXY& otherVec) const;
		VecXY    operator -() const;
//		VecXY    operator - (VecXY&& otherVec)      const;
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

// DecayChainExtractor is NOT thread safe
// DecayChainExtractor takes a Pythia event record and strips out the entire decay
// chains that are needed by APPro.

namespace Pythia8 {class Particle; class Event;};

class DecayChainExtractor
{
	private:
		static const Int_t PROCESSED_STATUS;

		static Pythia8::Event* event;
		static DelphesFactory* factory;
		static TObjArray* outputArray;
		static Int_t outputArrayCurrentSize;
		static void (*PythiaParticleToDelphesCandidate)(Pythia8::Particle const&, Candidate* const);

		static void AddDaughters(const Int_t particleIndex, Int_t pythiaIndexShift, Int_t motherThenLastDaughterIndex);

	public:
		static void AddFullDecayChain(Pythia8::Event& event_in, DelphesFactory* const factory_in,
			TObjArray* const outputArray_in, void (* const PythiaParticleToDelphesCandidate_in) (Pythia8::Particle const&, Candidate*));
};

#endif
