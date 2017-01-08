// Basic functions used in the
// dijet-hadron and jet-hadron
// correlation analysis
// Nick Elsey

// STL 
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <limits.h>
#include <unistd.h>

// fastjet 3
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TObjArray.h"
#include "TString.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "TBranch.h"
#include "TMath.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TStopwatch.h"

// TStarJetPico
#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"

#include "ktTrackEff.hh"

#ifndef CORRFUNCTIONS_HH
#define CORRFUNCTIONS_HH


namespace jetHadron {
  
  // Declare ahead of time
  // For correlation functions
  class histograms;
  
	// IO/OS MANIP Functions
  
	// Helper to build the TChain, used to decide which input format
	bool HasEnding (std::string const &full_string, std::string const &ending);
  
  // Other string manipulations
  // Used to check if a string begins with a substring
	bool BeginsWith (std::string const &full_string, std::string const &beginning);
  
  // Used to pull the current directory from its absolute path
  std::string GetDirFromPath( std::string path );
  
	// Used to find the path to current working directory
	// Used to make this all relatively machine independent
	std::string getPWD();
	
	// Returns mphi-vphi, in [ -pi, pi ]
	// Not used anymore- everything is done with
	// Fastjet::PseudoJets
	double GetdPhi ( double mphi, double vphi );
	
	// Checks there are the proper number of jets
	// Then calculates Aj
	double CalcAj ( std::vector<fastjet::PseudoJet>& jets );
  
  // Used to get reference centrality from gRefMult
  int GetReferenceCentrality( int gRefMult );
  
  // Used to get the inverse definition of reference centality ( 8->0, 7->1, etc)
  int GetReferenceCentralityAlt( int RefCent );
  
  // Find the vertex Z bin that corresponds to each Vz
  // Returns -1 if Vz outside of accepted range
  int GetVzBin( double Vz );
  
  // Converts TStarJetPicoVectors into PseudoJets
  void ConvertTStarJetVector( TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles, bool ClearVector = true, double towerScale = 1.0 );
  // applies an effective 90% relative efficiency compared to auau
  void ConvertTStarJetVectorPP( TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles, ktTrackEff& eff, unsigned seed, bool ClearVector = true, double towerScale = 1.0 );
  
  // Used in pp to convert either all embedding tracks or
  // only hard embedding tracks ( > 2.0 GeV )
  void ConvertTStarJetVectorPPEmbedded( TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles, bool allTracks = false, double towerScale = 1.0 );
  
  // Finds the triggers and saves them, if requireTrigger == True
  void GetTriggers( bool requireTrigger, TClonesArray* triggerObjs, std::vector<fastjet::PseudoJet> & triggers );
  
  // For the pp data where the trigger objects dont seem to be working
  void GetTriggersPP( bool requireTrigger, std::vector<fastjet::PseudoJet> ppParticles, std::vector<fastjet::PseudoJet>& triggers );
	
	// Summary of initial settings for dijet-hadron correlation
	void BeginSummaryDijet ( double jetRadius, double leadJetPtMin, double subLeadJetPtMin, double jetMaxPt, double hardJetConstPt, double softJetConstPt, int nVzBins, double VzRange, std::string dijetFile, std::string corrFile );
	
	// Summary of initial settings for jet-hadron correlation
	void BeginSummaryJet ( double jetRadius, double jetPtMin, double jetPtMax, double jetConstPt, int nVzBins, double VzRange, std::string jetFile, std::string corrFile );
	
	// Efficiency information post-analysis for dijet-hadron correlation
	void EndSummaryDijet ( int ntotal, int nviable, int nused, double time );
	
	// Efficiency information post-analysis for jet-hadron correlation
	void EndSummaryJet ( int ntotal, int nused, double time );
  
  // Initializes the TStarJetPicoReader, so we dont have
  // To have all that code hanging around in the analysis
  // Collision Type is 'AuAu' or 'pp'
  void InitReader( TStarJetPicoReader & reader, TChain* chain, std::string collisionType, std::string triggerString, double softwareTrigger, int nEvents );
  
  // Use this to decide if there are 2 dijets for dijet analysis in the proper pt ranges
  // Or for jet analysis if there is a single jet
  bool CheckHardCandidateJets( std::string analysisType, std::vector<fastjet::PseudoJet> & HiResult, double leadJetPtMin, double subJetPtMin );
  
  // Use this to select either one or two jets depending on the analysis type
  std::vector<fastjet::PseudoJet> BuildHardJets( std::string analysisType, std::vector<fastjet::PseudoJet> & HiResult );
  
  // Use this to match either hard single jets or hard dijets to full event jets
  // analysisType: dijet or jet
  // hardJets: jet(s) corresponding to analysisType
  // LoResult: the jets from clustering with the whole event
  // requireTrigger: whether or not the jets need to be matched to a trigger
  // triggers: the list of triggers
  // jetRadius: radius used for clustering
  std::vector<fastjet::PseudoJet> BuildMatchedJets( std::string analysisType, std::vector<fastjet::PseudoJet> & hardJets, std::vector<fastjet::PseudoJet> & LoResult, bool requireTrigger, std::vector<fastjet::PseudoJet> & triggers, double jetRadius = 0.4 );
  
  // Finally, correlation function -
  // It correlates leading and subleading jets
  // With the associated particle given the associated weight
  // Checks to make sure it is charged, has a proper efficiency,
  // And that the associated track is within our eta range
  
  // First, to check that the track makes all cuts
  bool useTrack( fastjet::PseudoJet& assocTrack, double efficiency );
  
  // Correlate Leading
  bool correlateLeading( std::string analysisType, int vzBin, int centBin, histograms* histogram, fastjet::PseudoJet& leadJet, fastjet::PseudoJet& assocTrack, double efficiency, double aj );
  
  // Correlate Subleading
  bool correlateSubleading( std::string analysisType, int vzBin, int centBin, histograms* histogram, fastjet::PseudoJet& subJet, fastjet::PseudoJet& assocTrack, double efficiency, double aj );
  
  // Correlate for jet-hadron
  bool correlateTrigger( std::string analysisType, int vzBin, int centBin, histograms* histogram, fastjet::PseudoJet& triggerJet, fastjet::PseudoJet& assocTrack, double efficiency );
	
	// FastJet functionality
	
  // Build default jet definitions, using kt and anti-kt for background and analysis, respectively
  fastjet::JetDefinition	AnalysisJetDefinition( double R );
  fastjet::JetDefinition	BackgroundJetDefinition( double R );
  
  // Constituent selectors for high and low constituent pt clustering
  fastjet::Selector	SelectLowPtConstituents( double trackMaxEta, double trackMinPtLow );
  fastjet::Selector SelectHighPtConstituents( double trackMaxEta, double trackMinPtHigh );
  
  // Candidate jet selectors
  fastjet::Selector SelectJetCandidates( double trackMaxEta, double jetRadius, double jetMinPt, double JetMaxPt );
  
  // Removes 2 hard jets from background estimation
  fastjet::Selector SelectBkgEstimator( double maxTrackRap, double jetRadius );
  
  // Ghosted area spec, used for background
  // Estimation and subtraction
  fastjet::GhostedAreaSpec GhostedArea( double trackMaxEta, double jetRadius );
  
  // Definition for the area estimation
  fastjet::AreaDefinition  AreaDefinition( fastjet::GhostedAreaSpec ghostAreaSpec );
  
  // --------------------------
  // ------ Event Mixing ------
  // --------------------------
  
  // Pulls analysis variables from the directory name
  int GetVarsFromString( std::string& analysisType, std::string analysisString, double& leadPt, double& subPt, double& maxPt, double& jetRadius, double& hardPt, bool& useEff, bool& reqTrigger, unsigned& binsEta, unsigned& binsPhi );
  
  // Used to decide what max pt can be for a jet in
  // HT events to still be used in mixing
  // Used for both jet-hadron and dijet-hadron
  double GetMixEventJetPtMax( bool useMB, std::string analysisType, double leadJetPtMin );
  
  // Decides whether an event should be used
  // In mixing or not - logic depends on analysis type
  bool UseEventInMixing( std::string analysisType, bool isMB, std::vector<fastjet::PseudoJet>& highPtConsJets, int refMult, int vzBin );
  
}

#endif
