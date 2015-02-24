#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  AllParticlePropagator

  MuonTrackingEfficiency
  MuonMomentumSmearing

  Calorimeter
  
  TowerJetFinder
  
  JetEnergyScale

  MissingET

  HighPtBTagger

  ScalarHT

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
                           sqrt(pow(.3706 / pt, 2.) + .00185 + pow(.0001809 * pt, 2.)) }
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
  
  set ECalSquareGranularity 1
  set ECalAbsEtaMax 3.2
  set HCalAbsEtaMax 4.9
  
  set pi [expr {acos(-1)}]

  # “add EtaPhiBins singleEta listOfPhi”
  # Automatically sorted because each add commmand goes to a std::map<eta, set<phi>>
  # First eta in map is most negative edge of most negative eta bin
  # Last eta in map is most positive edge of most positive eta bin
  # Each listOfPhi MUST start with -Pi and end with Pi, otherwise angle will be lost
  # The phi edges of a given bin determined by the eta of its UPPER EDGE
  
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
  # K0short should be properly decayed/propagated by AllParticlePropagator
  # If it hits the Cal it's a neutral hadron, don't need a special rule
  # add EnergyFraction {310} {0.0 1.0}
  # According to arXiv:0610128, about 80% of pions start showering before the ECAL
  # Of those that shower inside the ecal, they split their energy about equally
  # Assume this is the same for all other long-lived charged particles
  # (and lamba, which will also be properly decayed/propagated by AllParticlePropagator)
  add EnergyFraction {211} {0.4 0.6}
  add EnergyFraction {321} {0.4 0.6}
  add EnergyFraction {2212} {0.4 0.6}  
  add EnergyFraction {3122} {0.4 0.6}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0 0.0}
  add EnergyFraction {14} {0.0 0.0}
  add EnergyFraction {16} {0.0 0.0}
  add EnergyFraction {1000022} {0.0 0.0}
  add EnergyFraction {1000023} {0.0 0.0}
  add EnergyFraction {1000025} {0.0 0.0}
  add EnergyFraction {1000035} {0.0 0.0}
  add EnergyFraction {1000045} {0.0 0.0}
  

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

##########################
# Track pile-up subtractor
##########################

module TrackPileUpSubtractor TrackPileUpSubtractor {
# add InputArray InputArray OutputArray
  add InputArray Calorimeter/eflowTracks eflowTracks
  add InputArray ElectronEnergySmearing/electrons electrons
  add InputArray MuonMomentumSmearing/muons muons

  # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution
  # Z vertex resolution in m
  set ZVertexResolution 0.0001
}

####################
# Neutral tower merger
####################

module Merger NeutralTowerMerger {
# add InputArray InputArray
  add InputArray Calorimeter/eflowPhotons
  add InputArray Calorimeter/eflowNeutralHadrons
  set OutputArray eflowTowers
}

##################################
# Energy flow merger (all tracks)
##################################

module Merger EFlowMergerAllTracks {
# add InputArray InputArray
  add InputArray TrackMerger/tracks
  add InputArray Calorimeter/eflowPhotons
  add InputArray Calorimeter/eflowNeutralHadrons
  set OutputArray eflow
}


####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray Calorimeter/eflowTracks
  add InputArray Calorimeter/eflowPhotons
  add InputArray Calorimeter/eflowNeutralHadrons
  set OutputArray eflow
}

#############
# Rho pile-up
#############

module FastJetGridMedianEstimator Rho {

  set InputArray Calorimeter/towers
  set RhoOutputArray rho

  # add GridRange rapmin rapmax drap dphi
  # rapmin - the minimum rapidity extent of the grid
  # rapmax - the maximum rapidity extent of the grid
  # drap - the grid spacing in rapidity
  # dphi - the grid spacing in azimuth

  add GridRange -5.0 -2.5 1.0 1.0
  add GridRange -2.5 2.5 0.5 0.5
  add GridRange 2.5 5.0 1.0 1.0

}


#####################
# Neutrino Filter
#####################

module PdgCodeFilter NeutrinoFilter {

  set InputArray Delphes/stableParticles
  set OutputArray filteredParticles

  set PTMin 0.0

  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}

}

#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.6

  set JetPTMin 20.0
}

############
# Jet finder
############

module FastJetFinder FastJetFinder {
  set InputArray Calorimeter/towers

  set OutputArray jets

  # area algorithm: 0 Do not compute area, 1 Active area explicit ghosts, 2 One ghost passive area, 3 Passive area, 4 Voronoi, 5 Active area
  set AreaAlgorithm 5

  # jet algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.6

  set JetPTMin 20.0
}

###########################
# Jet Pile-Up Subtraction
###########################

module JetPileUpSubtractor JetPileUpSubtractor {
  set JetInputArray FastJetFinder/jets
  set RhoInputArray Rho/rho

  set OutputArray jets

  set JetPTMin 20.0
}

##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray JetPileUpSubtractor/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.0}
}

###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray Calorimeter/photons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for photons
  set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 10.0)  * (0.85) +
                         (abs(eta) > 2.5)                                   * (0.00)}
}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowMerger/eflow
  set RhoInputArray Rho/rho

  set OutputArray photons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.1
}

#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray TrackPileUpSubtractor/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
  set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 10.0)  * (0.85) +
                         (abs(eta) > 2.5)                                   * (0.00)}
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowMerger/eflow
  set RhoInputArray Rho/rho

  set OutputArray electrons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.1
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
  set InputArray TrackPileUpSubtractor/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency as a function of eta and pt}

  # efficiency formula for muons
  set EfficiencyFormula {                                      (pt <= 10.0) * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.7) * (pt > 10.0)  * (0.85) +
                         (abs(eta) > 2.7)                                   * (0.00)}
}

################
# Muon isolation
################

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
  set IsolationInputArray EFlowMerger/eflow
  set RhoInputArray Rho/rho

  set OutputArray muons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.1
}

###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMergerAllTracks/eflow
  set MomentumOutputArray momentum
}


##################
# Scalar HT merger
##################

module Merger ScalarHT {
# add InputArray InputArray
  add InputArray UniqueObjectFinder/jets
  add InputArray UniqueObjectFinder/electrons
  add InputArray UniqueObjectFinder/photons
  add InputArray UniqueObjectFinder/muons
  set EnergyOutputArray energy
}

###########
# b-tagging
###########

module BTagging BTagging {
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set BitNumber 0

  set DeltaR 0.5

  set PartonPTMin 1.0

  set PartonEtaMax 2.5

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
  # gluon's PDG code has the lowest priority

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}

  # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormula {4} {                                      (pt <= 15.0) * (0.000) +
                                                (abs(eta) <= 1.2) * (pt > 15.0) * (0.2*tanh(pt*0.03 - 0.4)) +
                              (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 15.0) * (0.1*tanh(pt*0.03 - 0.4)) +
                              (abs(eta) > 2.5)                                  * (0.000)}

  # efficiency formula for b-jets
  add EfficiencyFormula {5} {                                      (pt <= 15.0) * (0.000) +
                                                (abs(eta) <= 1.2) * (pt > 15.0) * (0.5*tanh(pt*0.03 - 0.4)) +
                              (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 15.0) * (0.4*tanh(pt*0.03 - 0.4)) +
                              (abs(eta) > 2.5)                                  * (0.000)}
}

module TauTagging TauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 2.5

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

module UniqueObjectFinder UniqueObjectFinder {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScale/jets jets
}

##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle

#  add Branch TrackMerger/tracks Track Track
  add Branch Calorimeter/towers Tower Tower

#  add Branch Calorimeter/eflowTracks EFlowTrack Track
#  add Branch Calorimeter/eflowPhotons EFlowPhoton Tower
#  add Branch Calorimeter/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch GenJetFinder/jets GenJet Jet
  add Branch UniqueObjectFinder/jets Jet Jet
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/photons Photon Photon
  add Branch UniqueObjectFinder/muons Muon Muon
  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
  add Branch Rho/rho Rho Rho
  add Branch PileUpMerger/vertices Vertex Vertex

}
