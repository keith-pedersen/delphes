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
 *  Being detected in the pixel detector relies on several efficiencies. These 
 *  are not currently imlpemented.
 *    a. When a particle physically strikes a pixel detector layer, a 4-vector describing
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
 *    b. After all the hits in a layer have been collected, they will be binned in certain 
 *         regions of interest. Before this is accomplished, bins will be deactivated 
 *         based on the layers *online* efficiency (to simulate dead pixels). 
 *           @ Online efficiency = The probability that any given pixel is functioning properly. 
 *           @ There are seperate online efficiencies for pixel and endcap layers
 *         There are also limitations with this model:
 *           @ The location of dead pixels is randomly determined before each binning.
 *           @ Besides the seperate efficiency for barrel and endcap, no eta
 *             dependence is simulated.
 *
 *  \author K. Pedersen - Illinois Institute of Technology
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
