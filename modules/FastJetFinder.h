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

#ifndef FastJetFinder_h
#define FastJetFinder_h

/** \class FastJetFinder
 *
 *  Finds jets using FastJet library.
 *
 *  \author P. Demin - UCL, Louvain-la-Neuve
 *
 */

#include "classes/DelphesModule.h"

#include <vector>

class TObjArray;
class TIterator;

namespace fastjet {
  class JetDefinition;
  class AreaDefinition;
  class JetMedianBackgroundEstimator;
  namespace contrib {
    class NjettinessPlugin;
  }
}

class FastJetFinder: public DelphesModule
{
public:

  FastJetFinder();
  ~FastJetFinder();

  void Init();
  void Process();
  void Finish();

private:

  void *fPlugin; //!
  void *fRecomb; //!
  fastjet::contrib::NjettinessPlugin *fNjettinessPlugin; //!

  fastjet::JetDefinition *fDefinition; //!

  Int_t fJetAlgorithm;
  Double_t fParameterR;
  Double_t fJetPTMin;
  Double_t fConeRadius;
  Double_t fSeedThreshold;
  Double_t fConeAreaFraction;
  Int_t fMaxIterations;
  Int_t fMaxPairSize;
  Int_t fIratch;
  Int_t fAdjacencyCut;
  Double_t fOverlapThreshold;

  //-- N (sub)jettiness parameters --

  Bool_t fComputeNsubjettiness;
  Double_t fBeta;
  Int_t fAxisMode;
  Double_t fRcutOff;
  Int_t fN ;

  //-- Trimming parameters --
  
  Bool_t fComputeTrimming;
  Double_t fRTrim;
  Double_t fPtFracTrim;
  
  //-- Pruning parameters --

  Bool_t fComputePruning;
  Double_t fZcutPrun;
  Double_t fRcutPrun;
  Double_t fRPrun;

  //-- SoftDrop parameters --

  Bool_t fComputeSoftDrop;
  Double_t fBetaSoftDrop;
  Double_t fSymmetryCutSoftDrop;
  Double_t fR0SoftDrop;

  // --- FastJet Area method --------

  fastjet::AreaDefinition *fAreaDefinition;
  Int_t fAreaAlgorithm;
  Bool_t  fComputeRho;

  // -- ghost based areas --
  Double_t fGhostEtaMax;
  Int_t fRepeat;
  Double_t fGhostArea;
  Double_t fGridScatter;
  Double_t fPtScatter;
  Double_t fMeanGhostPt;

  // -- voronoi areas --
  Double_t fEffectiveRfact;

  // -- Strategy parameter for runs which crash due to bug in FastJet tiling
  Int_t fStrategy;
  /*
  /// New in FJ3.1
  N2MHTLazy9AntiKtSeparateGhosts   = -10,
  /// only looks into a neighbouring tile for a particle's nearest
  /// neighbour (NN) if that particle's in-tile NN is further than the
  /// distance to the edge of the neighbouring tile. Uses tiles of
  /// size R and a 3x3 tile grid around the particle.
  /// New in FJ3.1
  N2MHTLazy9   = -7,
  /// Similar to N2MHTLazy9, but uses tiles of size R/2 and a 5x5 tile
  /// grid around the particle.
  /// New in FJ3.1
  N2MHTLazy25   = -6,
  /// Like to N2MHTLazy9 but uses slightly different optimizations,
  /// e.g. for calculations of distance to nearest tile; as of
  /// 2014-07-18 it is slightly slower and not recommended for
  /// production use. To considered deprecated.
  /// New in FJ3.1
  N2MHTLazy9Alt   = -5,
  /// faster that N2Tiled above about 500 particles; differs from it
  /// by retainig the di(closest j) distances in a MinHeap (sort of
  /// priority queue) rather than a simple vector.
  N2MinHeapTiled   = -4,
  /// fastest from about 50..500
  N2Tiled     = -3,
  /// legacy
  N2PoorTiled = -2,
  /// fastest below 50
  N2Plain     = -1,
  /// worse even than the usual N^3 algorithms
  N3Dumb      =  0,
  /// automatic selection of the best (based on N), including
  /// the LazyTiled strategies that are new to FJ3.1
  Best        =  1,
  /// best of the NlnN variants -- best overall for N>10^4.
  /// (Does not work for R>=2pi)
  NlnN        =  2,
  /// legacy N ln N using 3pi coverage of cylinder.
  /// (Does not work for R>=2pi)
  NlnN3pi     =  3,
  /// legacy N ln N using 4pi coverage of cylinder
  NlnN4pi     =  4,
  /// Chan's closest pair method (in a variant with 4pi coverage),
  /// for use exclusively with the Cambridge algorithm.
  /// (Does not work for R>=2pi)
  NlnNCam4pi   = 14,
  /// Chan's closest pair method (in a variant with 2pi+2R coverage),
  /// for use exclusively with the Cambridge algorithm.
  /// (Does not work for R>=2pi)
  NlnNCam2pi2R = 13,
  /// Chan's closest pair method (in a variant with 2pi+minimal extra
  /// variant), for use exclusively with the Cambridge algorithm.
  /// (Does not work for R>=2pi)
  NlnNCam      = 12, // 2piMultD
  /// the automatic strategy choice that was being made in FJ 3.0
  /// (restricted to strategies that were present in FJ 3.0)
  BestFJ30     =  21,
  /// the plugin has been used...
  plugin_strategy = 999
  */

#if !defined(__CINT__) && !defined(__CLING__)
  struct TEstimatorStruct
  {
    fastjet::JetMedianBackgroundEstimator *estimator;
    Double_t etaMin, etaMax;
  };

  std::vector< TEstimatorStruct > fEstimators; //!
#endif

  TIterator *fItInputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!
  TObjArray *fRhoOutputArray; //!

  ClassDef(FastJetFinder, 1)
};

#endif
