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
 *  4.  Helical propagation is skipped (in favor of a straight-line approximation)
 *  whenever the bending would be too small to be meaningful (1/20th of a degree is 
 *  the default threshold). This cuts down on the work done by this module, limiting
 *  the actual full propagations to relatively few particles in the actual chain.  
 * 
 *  5.  APPro is intended to be used with the following Pythia options
 *        ParticleDecays:limitRadius = on
 *        ParticleDecays:rMax = 20000 ! 20m of decay (change as desired)
 *        211:mayDecay = true !let charged pions decay
 *        321:mayDecay = true !let charged kaons decay
 *        130:mayDecay = true !let K_L0 decay
 *        13:mayDecay = false !don't let muons decay (too rare at normal energies)
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
 *  6.  Because the entire decay chain is needed, including dereferencable daughter
 *  indices, PileupMerger cannot be used (in its 3.2.0 implementation). The inclusion
 *  of Pileup is natively supported by "add"-ing an InputArray containing pileup 
 *  in the configuration card (i.e. like a Merger module). Such arrays MUST
 *  include the entire decay chain (with functional/re-indexed daughter indices),
 *  though they don't need all the colored precursors. 
 *     This will likely require a new Pileup picking module (i.e. one that outputs 
 *  not only final state particles, but the entire, pertinent decay chain). This module
 *  is on an To-do list, but an example re-indexing function is supplied.
 *
 *  7.  APPro absolutely NEEDS to know each particles proper lifetime, so this field
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
 *                              ATTENTION: This field is used as a flag for APPro, to 
 *                              determine when a Candidate has already been propagated. 
 *                              When it is negative, the particle has not been propagated
 *                              yet. Thus, it should be initialized to -1. 
 *      * Float_t  TrackLength       [in mm]. 
 *                              This is the total length of the particle's track as 
 *                              it moves throught the tracker. If is
 *                                   length = cTau * gamma * |p|/E = cTau * |p|/m = cTau*sqrt(gamma^2-1)
 *                              It can be viewed (roughly) as the number of track hits.
 *                              WARNING: No eta cut is applied by APPro; all propagating
 *                              particles will have a non-zero TrackLength, even if they 
 *                              are at eta = 10, and do not actually strike any sensors.
 * 
 *  9.  Instead of storing each Candidate's original, un-propagted Position by cloning the 
 *  the original Candidate and storing the clone in the fArray of the propagated Candidate,
 *  APPro stores the original position in the "Area" 4-vector (which is 
 *  unused by all Candidates except jets, which should not be propagated).
 *
 *  10. APPro does not currently simulate any energy loss from multiple scattering,
 *  bremsstrahlung, or cyclotron radiation. A simple parameterization
 *  of these effects would be a welcome addition.
 
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

class TLorentzVector;
class Candidate;

class DaughterRotation
{
	// This class
	// (i) ACTIVE-ly rotates (i.e. keeps spatial axes the same, moves the vector)
	// a Candidate's 4-momentum in the XY plane (energy and z-momentum is unaffected).
	//  * If the right hand curls from x to y, and the thumb points to z (a RH xyz frame),
	//    then an infinitessimal RH rotation is created by an infinitessimaly positive angle
	
	public:
		DaughterRotation(const Double_t theta_in);
		DaughterRotation(DaughterRotation const* const grandmother, const Double_t deltaTheta);
		
		void Rotate(Candidate* const candidate) const;
		
	private:
		const Double_t 
			theta,
			cosTheta,
			sinTheta;
};

class TObjArray;
class TIterator;

class AllParticlePropagator: public DelphesModule
{
	public:
		AllParticlePropagator();
		~AllParticlePropagator();

		void Init();
		void Process();
		void Finish();

	private:
		Candidate* Rotate(Candidate* const daughter, DaughterRotation const* const rotation, const TLorentzVector& mothersPosition);
		void NullifyDaughters(Candidate* const candidate, const bool nullifyThisParticle = false);
		bool Propagate(Candidate* const candidate, DaughterRotation const* rotation);

		Double_t fRadius, fRadius2;     //    assumed supplied in [meters], converted to [mm] by Init()
		                                // Radius of full tracker cylinder (and squared radius, for re-use)
		
		Double_t fHalfLength;           //    assumed supplied in [meters], converted to [mm] by Init()
		                                // Half the length of the full tracker cylinder
		
		Double_t fBz;                   //    assumed supplied in [Tesla]
		                                // Strengh of perflectly homogenous magnetic field
		                                // alligned with the +z direction
		
		Double_t fMinHelixAngle;        //    assumed supplied in [radians]
		                                // If a particle's total rotation, about its helical center,
		                                // is less than fMinimumHelixAngle (should be small),
		                                // the Candidate is propagated with a straight
		                                // line approximation
		                                
		Double_t fMinTrackLength;       //    assumed supplied in [meters], converted to [mm] by Init()
		                                // Used for a quick pre-sort of tracks which are entirely
		                                // too small, so that the charged OutputArrays aren't
		                                // filled with tiny, invisible tracks.

		// It would be nice to parameterize the energy loss, but this is not supported yet.
				
		// Input
		std::vector<TObjArray*> fInputList; //! use "add InputArray @$%!" in configuration card
		TObjArray* currentInputArray; //!
		
		// Output
		TObjArray* fOutputArray; //! Particles that intersect the cylinder
		TObjArray* fChargedHadronOutputArray; //! Charged hadron track Candidates
		TObjArray* fElectronOutputArray; //! Electron track Candidates
		TObjArray* fMuonOutputArray; //! Muon track candidates
				
		ClassDef(AllParticlePropagator, 1)
};

#endif
