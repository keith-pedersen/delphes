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


/** \class Isolation
 *
 *  Sums transverse momenta of isolation objects (tracks, calorimeter towers, etc)
 *  within a DeltaR cone around a candidate and calculates fraction of this sum
 *  to the candidate's transverse momentum. outputs candidates that have
 *  the transverse momenta fraction within (PTRatioMin, PTRatioMax].
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "modules/Isolation.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"

#include "TMath.h"
#include "TString.h"
#include "TFormula.h"
#include "TRandom3.h"
#include "TObjArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

class IsolationClassifier : public ExRootClassifier
{
public:

  IsolationClassifier() {}

  Int_t GetCategory(TObject *object);

  Double_t fPTMin;
};

//------------------------------------------------------------------------------

Int_t IsolationClassifier::GetCategory(TObject *object)
{
  Candidate *track = static_cast<Candidate*>(object);
  const TLorentzVector &momentum = track->Momentum;

  if(momentum.Pt() < fPTMin) return -1;

  return 0;
}

//------------------------------------------------------------------------------

Isolation::Isolation() :
  fClassifier(0), fFilter(0),
  fItIsolationInputArray(0), fItCandidateInputArray(0),
  fItRhoInputArray(0)
{
  fClassifier = new IsolationClassifier;
}

//------------------------------------------------------------------------------

Isolation::~Isolation()
{
}

//------------------------------------------------------------------------------

void Isolation::Init()
{
  const char *rhoInputArrayName;

  fDeltaRMax = GetDouble("DeltaRMax", 0.5);
  fDeltaRMax2 = fDeltaRMax*fDeltaRMax;
  
  fPTRatioMax = GetDouble("PTRatioMax", 0.1);

  fPTSumMax = GetDouble("PTSumMax", 5.0);

  fUsePTSum = GetBool("UsePTSum", false);

  fClassifier->fPTMin = GetDouble("PTMin", 0.5);

  // import input array(s)

  fIsolationInputArray = ImportArray(GetString("IsolationInputArray", "Delphes/partons"));
  fItIsolationInputArray = fIsolationInputArray->MakeIterator();

  fFilter = new ExRootFilter(fIsolationInputArray);

  fCandidateInputArray = ImportArray(GetString("CandidateInputArray", "Calorimeter/electrons"));
  fItCandidateInputArray = fCandidateInputArray->MakeIterator();

  rhoInputArrayName = GetString("RhoInputArray", "");
  if(rhoInputArrayName[0] != '\0')
  {
    fRhoInputArray = ImportArray(rhoInputArrayName);
    fItRhoInputArray = fRhoInputArray->MakeIterator();
  }
  else
  {
    fRhoInputArray = 0;
  }

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "electrons"));
}

//------------------------------------------------------------------------------

void Isolation::Finish()
{
  if(fItRhoInputArray) delete fItRhoInputArray;
  if(fFilter) delete fFilter;
  if(fItCandidateInputArray) delete fItCandidateInputArray;
  if(fItIsolationInputArray) delete fItIsolationInputArray;
}

//------------------------------------------------------------------------------

void Isolation::Process()
{
  Candidate *candidate, *isolation, *object;
  TObjArray *isolationArray;
  Double_t sum, ratio;
  Int_t counter;
  Double_t absEta = 0.0;
  Double_t rho = 0.0;

  if(fRhoInputArray && fRhoInputArray->GetEntriesFast() > 0)
  {
    candidate = static_cast<Candidate*>(fRhoInputArray->At(0));
    rho = candidate->Momentum.Pt();
  }

  // select isolation objects
  fFilter->Reset();
  isolationArray = fFilter->GetSubArray(fClassifier, 0);

  if(isolationArray == 0) return;

  TIter itIsolationArray(isolationArray);
  
  //KDP - Strip all neccessary information from IsolationArray into a vector
  //      of IsolatingObjects
  std::vector<IsolatingObject> newIsolationArray;
  
  while((isolation = static_cast<Candidate*>(itIsolationArray.Next())))
  {
    newIsolationArray.push_back(IsolatingObject(isolation));
  }

  // loop over all input jets
  fItCandidateInputArray->Reset();
  while((candidate = static_cast<Candidate*>(fItCandidateInputArray->Next())))
  {
    const TLorentzVector &candidateMomentum = candidate->Momentum;
    absEta = TMath::Abs(candidateMomentum.Eta());
   
    // loop over all input tracks
    sum = 0.0;
    counter = 0;
    for(std::vector<IsolatingObject>::const_iterator itIsolation = newIsolationArray.begin();
      itIsolation not_eq newIsolationArray.end(); ++itIsolation)
    {
    	if((itIsolation->DeltaR2(candidateMomentum) <= fDeltaRMax2)
    	   and (candidate->GetUniqueID() not_eq itIsolation->GetUniqueID()))
    	{
        sum += itIsolation->Pt();
        ++counter;
      }
    }
    
    // find rho
    rho = 0.0;
    if(fRhoInputArray)
    {
      fItRhoInputArray->Reset();
      while((object = static_cast<Candidate*>(fItRhoInputArray->Next())))
      {
        if(absEta >= object->Edges[0] && absEta < object->Edges[1])
        {
          rho = object->Momentum.Pt();
        }
      }
    }

    // correct sum for pile-up contamination
    sum = sum - rho*fDeltaRMax*fDeltaRMax*TMath::Pi();

    ratio = sum/candidateMomentum.Pt();
    if((fUsePTSum && sum > fPTSumMax) || ratio > fPTRatioMax) continue;

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------

IsolatingObject::IsolatingObject(Candidate const* const isolatingObject):
	uniqueID(isolatingObject->GetUniqueID()), 
	pt((isolatingObject->Momentum).Pt()), 
	eta((isolatingObject->Momentum).Eta()),
	phi((isolatingObject->Momentum).Phi())
{}

//------------------------------------------------------------------------------

const Double_t PI = acos(-1.);

Double_t IsolatingObject::DeltaR2(const TLorentzVector& candidateMomentum) const
{
	// This way of calculating deltaPhi doesn't require branches, so it is slightly faster 
	// than the conventional way (finding the absolute difference, then subtracting Pi 
	// if it's greater than Pi)
	const Double_t 
		dPhi = PI - abs(PI - abs(phi - candidateMomentum.Phi())),
		dEta = eta - candidateMomentum.Eta();
		
	return dPhi*dPhi + dEta*dEta;
}

//------------------------------------------------------------------------------

Double_t IsolatingObject::Pt() const {return pt;}

//------------------------------------------------------------------------------

UInt_t IsolatingObject::GetUniqueID() const {return uniqueID;}

//------------------------------------------------------------------------------
