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

#ifndef PixelPicture_h
#define PixelPicture_h

/** \class PixelPicture
 *
 *  Outputs a picture of the pixel clusters near the center of high-pt jets
 *
 *  \author K. Pedersen - Illinois Institute of Technology
 *
 */

#include "classes/DelphesModule.h"
#include <fstream>

class TObjArray;
class TLorentzVector;

class PixelPicture: public DelphesModule
{
	public:
		PixelPicture();
		~PixelPicture();

		void Init();
		void Process();
		void Finish();

	private:
		TObjArray* ImportPixelArray(const char* name);
		void PrintNewJet(const Double_t pt, const Double_t eta, const Double_t mass, const TLorentzVector& jetVec4);
		static int HadronFlavor(int pid);

		Double_t fMinPt;
		Double_t fMaxAbsEta;
		Double_t fMaxMass;
		Double_t fDeltaR;
				
		TObjArray* fJetInputArray;
		TIterator* fItJetInputArray;
		
		TObjArray* fBarrelLayers;
		TIterator* fItBarrelLayers;
		
		TObjArray* fEndcapLayers;
		TIterator* fItEndcapLayers;
		
		std::fstream coreSnapshotFile;

		ClassDef(PixelPicture, 1) 
};

#endif
