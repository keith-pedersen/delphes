#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  AllParticlePropagator

  MuonTrackingEfficiency
  MuonMomentumSmearing

  Calorimeter
  CalorimeterMerger
  
  TowerJetFinder
  
  JetEnergyScale

  MissingETandScalarHT

  HighPtBTagger

  TreeWriter
}

#################################
# Propagate particles in cylinder
#################################

module AllParticlePropagator AllParticlePropagator {
  add InputArray Delphes/allParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius_m 1.15
  # half-length of the magnetic field coverage, in m
  set HalfLength_m 3.51
  
  set MinTrackLengthRatio .7

  # magnetic field
  set Bz_Tesla 2.0
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray AllParticlePropagator/muons
  set OutputArray muons

  # Parameterization of reconstruction efficiency for ATLAS STANDALONE muons, as a function of eta and pt
  # It is ASSUMED here that the detector will trigger on Ht (di-jets), not muon pt, so we do not include trigger efficiency. 
  # The paper "Recent performance results with the ATLAS Muon Spectrometer" (G.Z.D. Porta) indicates that, above 10 GeV
  # reconstruction efficiency is totally flat (95%). This is corroborated by other sources. This will allow us to use 
  # <https://twiki.cern.ch/twiki/bin/view/AtlasPublic/ApprovedPlotsMuon> fig. 5 (efficiency vs. eta for 100 GeV muon)
  # to characterize the eta dependence of the efficency. It shows that detector services degrade reconstruction efficiency
  # up until |eta| > 0.12. The efficiency ramp up is fairly linear. There are also slight dips in efficiency at |eta| = .48
  # and |eta| = 78, which are small enough to ignore. In the barrel endcap transition region there is a large dip which can be
  # approximated as a double line. 
  # The total efficiency is broken into two parts (eta sum) * (pt sum). Each part is the sum of mutually exlusive 
  # elements, so that only one term contributes from each sum for any given muon. Their product is the total efficiency. 

  # tracking efficiency formula for muons, connect the dots linearly
  #
  #    eta (fig. 5)
  # {0.00, 0.00} -> {0.12, 1.00} -> {1.10, 1.00} -> {1.22, 0.65} -> {1.32, 1.00} -> {2.65, 1.00} -> {2.70, 0.92}
  #
  #    pt (fig. 4)
  # {2.0, 0.0} -> {5.0, 0.82} -> {10.0, 0.95} -> {100, 0.99} -> {infinity, 0.99}
  
  set EfficiencyFormula {  (                    (abs(eta) <= 0.12) * (abs(eta) / 0.12)               + 
                             (abs(eta) > 0.12 && abs(eta) <= 1.10) *   1.0                           +
                             (abs(eta) > 1.10 && abs(eta) <= 1.22) * ( 4.20833 - 2.91667 * abs(eta)) +
                             (abs(eta) > 1.22 && abs(eta) <= 1.32) * (-3.62    +     3.5 * abs(eta)) +  
                             (abs(eta) > 1.32 && abs(eta) <= 2.65) *   1.0                           +
                             (abs(eta) > 2.65 && abs(eta) <= 2.70) * ( 5.24    -     1.6 * abs(eta)) +
                             (abs(eta) > 2.70)                     * 0.0                               )
                          *(             (pt <= 2.0) * 0.0                            +
                             (pt > 2.0 && pt <= 5.0) * (-0.546667 +    0.273333 * pt) + 
                             (pt > 5.0 && pt <= 10.) * ( 0.69     +       0.026 * pt) + 
                             (pt > 10. && pt <= 100.) * ( 0.945556 + 0.000444444 * pt) + 
                             (pt > 100.0)             * 0.99 ) }
}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons
  
  set MassHypothesis .106

  # Parameterization of muon pt resolution for STANDALONE muons, based on
  # <https://twiki.cern.ch/twiki/bin/view/AtlasPublic/ApprovedPlotsMuon>
  # To get the overall resolution, figs. 53(c) and 53(d) were used. Since the resolution (with fit params)
  # was split between large and small sectors, and large sectors make up ~60% of the azimuthal angle, 
  # the parameters were averaged 60/40 between large and small parameters. Fig. 1 indicates that pt
  # resoution is approximately flat at all eta, except for a bump to twice resolution at the transition
  # region. Instead of fully replicating the jaggedness of the jump, it was treated as a plateau.
  # Fig. 1 also indicates that we can use the same resolution distribution in both barrel and endcap.

  # resolution formula for muons, connect the dots linearly (assume upper bound set by tracking efficiency)
  #
  #    eta (fig. 1)
  # {0.00, 1.0} -> {1.10, 1.0} -> {1.20, 2.25} -> {1.65, 2.25} -> {1.90, 1.00} -> {infinity, 1.00}
  # 
  #    pt  (fig.  53(c), 53(d))
  # error = sqrt((p0/pt)**2 + (p1**2) + (p2*pt)**2)
  #     
  #      Large (60%)  Small (40%)  Average
  # p0 |  3.132e-1  |  4.280e-1  |  3.706e-1
  # p1 |  3.882e-2  |  3.887e-2  |  3.884e-2
  # p2 |  2.482e-4  |  1.136e-4  |  1.809e-4
  set ResolutionFormula {  (                    (abs(eta) <= 1.10) *   1.0                      + 
                             (abs(eta) > 1.10 && abs(eta) <= 1.20) * (-12.75 + 12.5 * abs(eta)) +
                             (abs(eta) > 1.20 && abs(eta) <= 1.65) *   2.25                     +
                             (abs(eta) > 1.65 && abs(eta) <= 1.90) * ( 10.5  -   5. * abs(eta)) +  
                             (abs(eta) > 1.90)                     *   1.0
                           )
                             *
                           sqrt(pow(.3706 / pt, 2.) + 1.509e-3 + pow(1.809e-4 * pt, 2.)) }
}

#############
# Calorimeter
#############

module Calorimeter Calorimeter {
  set ParticleInputArray AllParticlePropagator/stableParticles
  set TrackInputArray MuonMomentumSmearing/muons

  set TowerOutputArray towers
  set PhotonOutputArray photons
  
  set EFlowTrackOutputArray eflowTracks
  set EFlowPhotonOutputArray eflowPhotons
  set EFlowNeutralHadronOutputArray eflowNeutralHadrons

  set ECalEnergyMin 0.5
  set HCalEnergyMin 1.0

  set ECalEnergySignificanceMin 1.0
  set HCalEnergySignificanceMin 1.0

  set SmearTowerCenter false
  set FindCalPhotons false
  
  set pi [expr {acos(-1)}]

  # “add EtaPhiBins singleEta listOfPhi”
  # 1. Automatically sorted because each add commmand goes to a std::map<eta, set<phi>>
  #    1a. First eta in map is most negative edge of most negative eta bin
  #    1b. Last eta in map is most positive edge of most positive eta bin
  # 2. Each listOfPhi MUST start with -Pi and end with Pi, otherwise angle will be lost
  # 3. The phi edges of a given bin determined by the eta of its UPPER EDGE
  
  # Here I give the HCAL the same segmentation of the ECAL, to use extra pointing
  # information of the ECAL. This is mostly correct, since the ECAL and HCAL distinction 
  # is not really meant to differentiate the ECAL and HCAL as seperate calorimeters, but 
  # to differentiate the response of the *combined* system to EM vs. Hadronic particles.
  # Most hadrons start showering in the ECAL, and the quoted HCAL resolutions account for this.
  # This method does give too much angular resolution to neutral hadrons (neutrons and K0_L),
  # which infrequently start showering in the ECAL. More importantly, however, it 
  # applies the HCAL energy resolution to too small a grid.
  
  # .025 x .025 towers
  # [-3.175 to 3.2]
  set SquareWidth .025
  set EtaMin -3.175
  set EtaMax 3.2
  
  set PhiBins {}
  set HalfNumPhiBins [expr {int(round($pi/$SquareWidth))}]
  for {set i -$HalfNumPhiBins} {$i <= $HalfNumPhiBins} {incr i} {add PhiBins [expr {$i * $pi/$HalfNumPhiBins}] }
    
  set EtaBins {}
  set NumEtaBins [expr {int(round(($EtaMax - $EtaMin)/$SquareWidth))}]
  for {set i 0} {$i <= $NumEtaBins} {incr i} {add EtaBins [expr {$EtaMin + $i * $SquareWidth}] }
  
  foreach eta $EtaBins {add EtaPhiBins $eta $PhiBins}
  
  
  # .1 x .1 towers
  # [-4.9 to -3.2] && [3.3 to 4.9]
  set SquareWidth .1
  set EtaMin -4.9
  set EtaMax -3.2
  
  set EtaBins {}
  set NumEtaBins [expr {int(round(($EtaMax - $EtaMin)/$SquareWidth))}]
  for {set i 0} {$i <= $NumEtaBins} {incr i} {add EtaBins [expr {$EtaMin + $i * $SquareWidth}] }
  
  set EtaMin 3.3
  set EtaMax 4.9
  
  set NumEtaBins [expr {int(round(($EtaMax - $EtaMin)/$SquareWidth))}]
  for {set i 0} {$i <= $NumEtaBins} {incr i} {add EtaBins [expr {$EtaMin + $i * $SquareWidth}] }
  
  set PhiBins {}
  set HalfNumPhiBins [expr {int(round($pi/$SquareWidth))}]
  for {set i -$HalfNumPhiBins} {$i <= $HalfNumPhiBins} {incr i} {add PhiBins [expr {$i * $pi/$HalfNumPhiBins}] }
  
  foreach eta $EtaBins {add EtaPhiBins $eta $PhiBins}
  
  

  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  # Used primarily by neutral hadrons
  add EnergyFraction {0} {0.0 1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0 0.0}
  add EnergyFraction {22} {1.0 0.0}
  add EnergyFraction {111} {1.0 0.0}
  # energy fractions for muon
  add EnergyFraction {13} {0.0 0.0}
  # energy fractions for neutrinos and neutralinos
  add EnergyFraction {12} {0.0 0.0}
  add EnergyFraction {14} {0.0 0.0}
  add EnergyFraction {16} {0.0 0.0}
  add EnergyFraction {1000022} {0.0 0.0}
  add EnergyFraction {1000023} {0.0 0.0}
  add EnergyFraction {1000025} {0.0 0.0}
  add EnergyFraction {1000035} {0.0 0.0}
  add EnergyFraction {1000045} {0.0 0.0}
  # K0short and Lambda should be properly decayed/propagated by AllParticlePropagator
  # If they hit they're a hadron, don't need a special rule
  # add EnergyFraction {310} {0.3 0.7}
  # add EnergyFraction {3122} {0.3 0.7}


  # set ECalResolutionFormula {resolution formula as a function of eta and energy}
  # http://arxiv.org/pdf/physics/0608012v1 jinst8_08_s08003
  # http://villaolmo.mib.infn.it/ICATPP9th_2005/Calorimetry/Schram.p.pdf
  # http://www.physics.utoronto.ca/~krieger/procs/ComoProceedings.pdf
  set ECalResolutionFormula {                  (abs(eta) <= 3.2) * sqrt(energy^2*0.0017^2 + energy*0.101^2) +
                             (abs(eta) > 3.2 && abs(eta) <= 4.9) * sqrt(energy^2*0.0350^2 + energy*0.285^2)}

  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  # http://arxiv.org/pdf/hep-ex/0004009v1
  # http://villaolmo.mib.infn.it/ICATPP9th_2005/Calorimetry/Schram.p.pdf
  set HCalResolutionFormula {                  (abs(eta) <= 1.7) * sqrt(energy^2*0.0302^2 + energy*0.5205^2 + 1.59^2) +
                             (abs(eta) > 1.7 && abs(eta) <= 3.2) * sqrt(energy^2*0.0500^2 + energy*0.706^2) +
                             (abs(eta) > 3.2 && abs(eta) <= 4.9) * sqrt(energy^2*0.09420^2 + energy*1.00^2)}
}

module Merger CalorimeterMerger {
	add InputArray Calorimeter/towers
	add InputArray Calorimeter/photons
	add InputArray MuonMomentumSmearing/muons
	
	set OutputArray clusteringObjects
}

############
# Jet finder
############

module FastJetFinder TowerJetFinder {
  set InputArray CalorimeterMerger/clusteringObjects

  set OutputArray jets

  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 0

  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 20.0
}

##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray TowerJetFinder/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.0}
}

###################
# Missing ET merger
###################

module Merger MissingETandScalarHT {
# add InputArray InputArray
  add InputArray Calorimeter/towers
  add InputArray Calorimeter/photons
  add InputArray MuonMomentumSmearing/muons

  set MomentumOutputArray momentum
  set EnergyOutputArray energy
}

###########
# b-tagging
###########

module HighPtBTagger HighPtBTagger {
  set JetInputArray JetEnergyScale/jets

  set BitNumber             3
  set MinJetPt            500
  set MinMuonPt            10
  set MinTowerPtRatio     0.05
  set CoreAntiktR         0.04
  set MinCorePtRatio      0.5
  set CoreMassHypothesis  3.0
  set MinFinalMass        0.
  set MaxFinalMass        1000.
  set PercentSolidAngle   .9
  
}

##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
#  add Branch Delphes/allParticles Particle GenParticle

#  add Branch TrackMerger/tracks Track Track
  add Branch Calorimeter/towers Tower Tower

#  add Branch Calorimeter/eflowTracks EFlowTrack Track
#  add Branch Calorimeter/eflowPhotons EFlowPhoton Tower
#  add Branch Calorimeter/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch JetEnergyScale/jets Jet Jet
  add Branch MuonMomentumSmearing/muons Muon Muon
  add Branch MissingETandScalarHT/momentum MissingET MissingET
  add Branch MissingETandScalarHT/energy ScalarHT ScalarHT
}
