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

#ifndef BoostedBTagger_h
#define BoostedBTagger_h

/** \class BoostedBTagger
 *
 *  Tags high-pt jets with a heavy flavor tag (based on arXiv:1307.1820)
 *
 *  1. Look for jets passing a pt cut.
 *       1a.    (ptJet >= fMinJetPt)
 *  2. Look inside for at least one muon passing a pt cut. This helps ensure
 *     that the the muon is well reconstructed.
 *       2a.    (ptMu >= fMinMuonPt)
 *  3. Tag jets where the muon resides with deltaR of the jet's centroid
 *       3a.    (deltaR(mu, jet) <= fMaxDeltaR)
 *
 *  \author K. Pedersen - Illinois Institute of Technology
 *
 */

#include "classes/DelphesModule.h"

class TObjArray;
class TIterator;
//class Candidate;
//class TFile;

class BoostedBTagger: public DelphesModule
{
	public:

		BoostedBTagger();
		~BoostedBTagger();

		void Init();
		void Process();
		void Finish();

	private:
		Int_t fBitNumber;
		Int_t fMaxJetRank;
		Double_t fMinJetPt;
		Double_t fMinMuonPt;
                Double_t fMaxDeltaR;

		TObjArray const* fJetInputArray; //!
		TIterator* fItJetInputArray; //!
                
                TObjArray const* fAllParticles;
                
                UInt_t static Flavor(Int_t PID);

		ClassDef(BoostedBTagger, 1)                
                
};

//bool SortCandidatePt_Low2High(Candidate const* one, Candidate const* two);

#endif
