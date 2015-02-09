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


/** \class PixelPicture
 *
 *  Outputs a picture of the pixel clusters near the center of high-pt jets
 *
 *  \author K. Pedersen - Illinois Institute of Technology
 *
 */

#include "modules/PixelPicture.h"
#include "modules/PropagatorAndPixelTracker.h"
#include "classes/DelphesClasses.h"

#include <cmath>
#include <vector>
#include <sstream>
#include <iostream>

#include "TObjArray.h"
#include "TLorentzVector.h"

//#include "classes/DelphesFactory.h"
//#include "classes/DelphesFormula.h"

//#include "ExRootAnalysis/ExRootResult.h"
//#include "ExRootAnalysis/ExRootFilter.h"
//#include "ExRootAnalysis/ExRootClassifier.h"

//#include "TMath.h"
//#include "TString.h"
//#include "TFormula.h"
//#include "TRandom3.h"

//#include "TDatabasePDG.h"

//#include <algorithm> 
//#include <stdexcept>
//#include <iostream>


using namespace std;

//------------------------------------------------------------------------------

PixelPicture::PixelPicture():
	fMinPt(0.), fMaxAbsEta(0.), fMaxMass(0.), fDeltaR(0.), 
	fJetInputArray(0), fItJetInputArray(), fBarrelLayers(0), fItBarrelLayers(),
	fEndcapLayers(0), fItEndcapLayers(), coreSnapshotFile()
{}

//------------------------------------------------------------------------------

PixelPicture::~PixelPicture()
{}

//------------------------------------------------------------------------------

void PixelPicture::Init()
{
	fMinPt = abs(GetDouble("MinPt", 1000.));
	fMaxAbsEta = abs(GetDouble("MaxAbsEta", 1.2));
	fMaxMass = abs(GetDouble("MaxMass", 100.));
	fDeltaR = abs(GetDouble("DeltaR", 0.05));
	
			
	fJetInputArray = ImportArray(GetString("JetInputArray", "JetFinder/jets"));
	fItJetInputArray = fJetInputArray->MakeIterator();
	
	fBarrelLayers = ImportPixelArray(GetString("PixelBarrelInputArray", "PropagatorAndPixelTracker/barrelLayers"));
	fItBarrelLayers = fBarrelLayers->MakeIterator();
	
	fEndcapLayers = ImportPixelArray(GetString("PixelEndcapInputArray", "PropagatorAndPixelTracker/endcapLayers"));
	fItEndcapLayers = fEndcapLayers->MakeIterator();
	
	coreSnapshotFile.open(GetString("CoreSnapshotFile", "jetCoreSnapshot.txt"), 
		                        std::fstream::out | std::fstream::trunc);
		                        
	if(not coreSnapshotFile.is_open())
		throw runtime_error("Cannot open CoreSnapShotFile");
}

//------------------------------------------------------------------------------

TObjArray* PixelPicture::ImportPixelArray(const char* name)
{
	TObjArray *object = static_cast<TObjArray *>(GetObject(Form("Pixel/%s", name), TObjArray::Class()));
	
	if(not object)
	{
		stringstream message;
		message << "can't access input list '" << name;
		message << "' in module '" << GetName() << "'";
		throw runtime_error(message.str());
	}

	return object;
}

//------------------------------------------------------------------------------

void PixelPicture::Finish()
{
	delete fItJetInputArray;
	delete fItBarrelLayers;
	delete fItEndcapLayers;
	coreSnapshotFile.close();
}

//------------------------------------------------------------------------------

void PixelPicture::Process()
{
	fJetInputArray->Sort(); // Sort by pT, high to low
	fItJetInputArray->Reset(); // Reset iterator to front of jet  list
	Candidate const* jetOfInterest;
	
	const string inBetweenLayers = "  |";
	const int BUFFER_SIZE = 64;
	
	while((jetOfInterest = static_cast<Candidate const*>(fItJetInputArray->Next())))
	{
		const TLorentzVector& jetVec4 = jetOfInterest->Momentum;
	
		const Double_t pt = jetVec4.Pt();
		if(pt >= fMinPt)
		{
			const Double_t eta = jetVec4.Eta();
			if(abs(eta) <= fMaxAbsEta)
			{
				const Double_t mass = jetVec4.M();
				if(mass <= fMaxMass)
				{
					PrintNewJet(pt, eta, mass, jetVec4);
				
					// Assume the energy, on average, travels in a straight line
					// Find the normalized (cylindrically) momentum vector
					const TVector3 jetP3N = (jetVec4.Vect())*=(1./pt);
					const VecXY rIntercept(jetP3N.x(), jetP3N.y());
					const RotationXY jetRotation(rIntercept.x, rIntercept.y);
											
					// Loop through barrel layers	
					fItBarrelLayers->Reset();
					PixelBarrel const* barrelLayer;
				
					// Assume that layers are processed from outside in
					// Plot z on the x-axis
					// Plot rPhi on the y-axis
				
					vector<string> lineBuffer;
					bool lineBufferUninitialized = true;
					int maxRows = 0;
			
					while((barrelLayer = static_cast<PixelBarrel*>(fItBarrelLayers->Next())))
					{
						const Double_t lorentzDrift = -(barrelLayer->thickness)*tan(barrelLayer->lorentzAngle);		
						const Double_t zIntercept = jetP3N.z() * (barrelLayer->radius);						
					
						// Find number of rPhi pixels to include
						const Double_t radiansPerPixel = (barrelLayer->rPhiWidth) / (barrelLayer->radius);
						const int rPhiPixelDomain = 2*(floor(fDeltaR/radiansPerPixel + 0.5) - 0.5) +.1; // Add .1 for rounding safety
						const Double_t deltaPhi = rPhiPixelDomain * (radiansPerPixel / 2);
						const Double_t cosDeltaPhiR = cos(deltaPhi) * (barrelLayer->radius);
					
						if(lineBufferUninitialized)
						{
							lineBufferUninitialized = false;
							maxRows = rPhiPixelDomain + 2;
							lineBuffer.assign(maxRows, "");
						}
					
						const int rowShift = (maxRows - rPhiPixelDomain)/2;
					
						// Find number of z pixels to include
						const Double_t 
							zUpperLimit = (floor(((barrelLayer->radius)*sinh(eta + fDeltaR) - zIntercept)/(barrelLayer->zWidth) + 0.5) - 0.5) * (barrelLayer->zWidth) + zIntercept,
							zLowerLimit = (ceil(((barrelLayer->radius)*sinh(eta - fDeltaR) - zIntercept)/(barrelLayer->zWidth) - 0.5) + 0.5) * (barrelLayer->zWidth) + zIntercept;
						const int zPixelDomain = (zUpperLimit - zLowerLimit)/(barrelLayer->zWidth) + .1;  // Add .1 for rounding safety
					
						// hitBin[rPhi][z]
						int hitBin[rPhiPixelDomain + 2*BUFFER_SIZE][zPixelDomain + 2*BUFFER_SIZE]; // Add a buffer bin on each side, just in case
						for(int rPhi = 0; rPhi < rPhiPixelDomain; ++rPhi)
						{
							for(int z = 0; z < zPixelDomain; ++z)
								hitBin[rPhi + BUFFER_SIZE][z + BUFFER_SIZE] = 0;
						}
					
						// Collect all nearby hits
						const std::vector<PixelHit>& hitsInLayer = barrelLayer->GetHits();
						for(std::vector<PixelHit>::const_iterator itHit = hitsInLayer.begin(); itHit != hitsInLayer.end(); ++itHit)
						{
							const PixelHit& hit = *itHit;
							if((hit.z < zUpperLimit) &&
								 (hit.z > zLowerLimit) && 
								 (rIntercept.Dot(hit.r) >= cosDeltaPhiR))
							{
								const Double_t deltaPhi = asin((jetRotation.ReverseRotateCopy(hit.r)).y / (barrelLayer->radius));
								
								const Double_t phiStart = (fDeltaR + deltaPhi);
								const Double_t phiEnd = phiStart + ((barrelLayer->thickness)*((hit.r).Cross(hit.rBeta))/((hit.r).Dot(hit.rBeta)) + lorentzDrift)/(barrelLayer->radius);
								const Double_t zStart = hit.z - zLowerLimit;
								const Double_t zEnd = zStart + (barrelLayer->thickness) * (hit.zBeta / (hit.rBeta).Norm());
								
								int rPhiIndexStart = phiStart/radiansPerPixel + 0.5;
								int rPhiIndexEnd = phiEnd/radiansPerPixel + 0.5;
								int zIndexStart = zStart/(barrelLayer->zWidth) + 0.5;
								int zIndexEnd = zEnd/(barrelLayer->zWidth) + 0.5;
								
								if(rPhiIndexStart > rPhiIndexEnd) std::swap(rPhiIndexStart, rPhiIndexEnd);
								if(zIndexStart > zIndexEnd) std::swap(zIndexStart, zIndexEnd);
																								
								const int flavor = HadronFlavor((hit.particle)->PID);
								int fillValue = 1;
								if(flavor == 13)
									fillValue = 1000000;
								else if(flavor == 5)
									fillValue = 10000;
								else if(flavor == 4)
									fillValue = 100;
									
								for(int rPhiIndex = rPhiIndexStart; rPhiIndex <= rPhiIndexEnd; ++rPhiIndex)
								{
									for(int zIndex = zIndexStart; zIndex <= zIndexEnd; ++zIndex)
										hitBin[rPhiIndex + BUFFER_SIZE][zIndex + BUFFER_SIZE] += fillValue;
								}
							}
						}
					
						// All hits have been binned (no Lorentz Angle yet)
						const string horizontalBuffer = (string(zPixelDomain, '~') += "|");
					
						lineBuffer[rowShift - 1] += (inBetweenLayers + horizontalBuffer);
					
						for(int rPhi = 0; rPhi < rPhiPixelDomain; ++rPhi)
						{
							string& line = lineBuffer[rPhi + rowShift];
						
							line += inBetweenLayers;
							for(int z = 0; z < zPixelDomain; ++z)
							{
								const int pixelValue = hitBin[rPhi + BUFFER_SIZE][z + BUFFER_SIZE];
								if(pixelValue >= 1000000)
									line += "M";
								else if(pixelValue >= 10000)
									line += "B";
								else if(pixelValue >= 100)
									line += "C";
								else if(pixelValue >= 10)
									line += "%";
								else if (pixelValue > 0)
									line += to_string(pixelValue);
								else
									line += " ";
							}
							line += "|";
						}
					
						lineBuffer[rowShift + rPhiPixelDomain] += (inBetweenLayers + horizontalBuffer);
					}
				
					for(int row = 0; row < maxRows; ++row)
					{
						coreSnapshotFile << lineBuffer[row] << "\n";
					}
				
					coreSnapshotFile << "\n\n\n";
				}
			}
		}
	}
}


//------------------------------------------------------------------------------

void PixelPicture::PrintNewJet(const Double_t pt, const Double_t eta, const Double_t mass, const TLorentzVector& jetVec4)
{
	// Create a hard break and give basic information about the jet
	char buffer[128];
	
	coreSnapshotFile << "=========================================================================================================\n";

	sprintf(buffer, "| %7s | %7s | pt(GeV) | mass(GeV) |\n", "eta", "phi");
	coreSnapshotFile << buffer;
	
	sprintf(buffer, "| %7.3f | %7.3f | %7.1f | %9.1f |\n\n\n", eta, jetVec4.Phi(), pt, mass);
	coreSnapshotFile << buffer;
}

//------------------------------------------------------------------------------

int PixelPicture::HadronFlavor(int pid)
{
	const int absPID = abs(pid);
	
	if(pid == 13)
		return 13;
	else if ((pid > 1000) or (pid < 10000))
		return absPID/100;
	else
		return (absPID%1000)/100;
}
