//hello
// Produces AuAu dijet-hadron
// Or monojet-hadron correlations
// Useage defined in submit/auau_correlation.csh
// Nick Elsey - 08.04.2016

// The majority of the jetfinding
// And correlation code is located in
// corrFunctions.hh
#include "corrFunctions.hh"

// All reader and histogram settings
// Are located in corrParameters.hh
#include "corrParameters.hh"

// Histogram class used for storing
// correlation histograms in bins of
// Aj, Vz and Centrality plus event info
#include "histograms.hh"

// ROOT is used for histograms and
// As a base for the TStarJetPico library
// ROOT Headers
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
#include "TRandom3.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TSystem.h"

// Make use of std::vector,
// std::string, IO and algorithm
// STL Headers
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <limits.h>
#include <unistd.h>

// Data is read in by TStarJetPico
// Library, we convert to FastJet::PseudoJet
// TStarJetPico headers
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

// The analysis is run on FastJet::PseudoJets
// We make use of the jetfinding tools
// And the convenient FastJet::Selectors
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

// Used for year 7 tracking efficiency corrections,
// if they are being used
#include "ktTrackEff.hh"

// -------------------------
// Command line arguments: ( Defaults
// Defined for debugging in main )
// [0]: Choose between dijet-hadron correlations or jet-hadron correlations
//      choices: dijet || jet
// [1]:	Choose to use particle efficiency corrections or not: true/false
// [2]: trigger coincidence: require a trigger with the leading jet
// [3]: software trigger: require a software trigger above this value in the event
// [4]: subleading jet pt min ( not used if doing jet-hadron correlations )
// [5]: leading jet pt min
// [6]: jet pt max
// [7]: jet radius, used in the jet definition
// [8]: hard pt cut
// [9]: number of bins for eta in correlation histograms
// [10]: number of bins for phi in correlation histograms
// [11]: output directory
// [12]: name for the correlation histogram file
// [13]: name for the dijet TTree file
// [14]: input file: can be a single .root or a .txt or .list of root files

// DEF MAIN()
int main ( int argc, const char** argv ) {
  
  // First check to make sure we're located properly
  std::string currentDirectory = jetHadron::getPWD( );
  
  // If we arent in the analysis directory, exit
  if ( !(jetHadron::HasEnding ( currentDirectory, "jet_hadron_corr" ) || jetHadron::HasEnding ( currentDirectory, "jet_hadron_correlation" )) ) {
    std::cerr << "Error: Need to be in jet_hadron_corr directory" << std::endl;
    return -1;
  }
  
  //Start a timer
  TStopwatch TimeKeeper;
  TimeKeeper.Start( );
  
  // Histograms will calculate gaussian errors
  // -----------------------------------------
  TH1::SetDefaultSumw2( );
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );
  
  
  // Read in command line arguments
  // ------------------------------
  // Defaults
  std::string 	executable    = "./bin/auau_correlation"; // placeholder
  std::string 	analysisType  = "dijet";									// choose dijet or jet ( decides value of requireDijets )
  bool					requireDijets	= true;											// this switches between dijet-hadron and jet-hadron
  bool					useEfficiency = true;										// choose to use particle-by-particle efficiency
  bool					requireTrigger= true;											// require leading jet to be within jetRadius of a trigger tower
  double        softwareTrig  = 6.0;                     // require there to be an offline trigger with E > softwareTrig
  double 				subJetPtMin   = 10.0;											// subleading jet minimum pt requirement
  double 				leadJetPtMin  = 20.0;											// leading jet minimum pt requirement
  double				jetPtMax			= 100.0;										// maximum jet pt
  double				jetRadius 		= 0.4;                      // jet radius for jet finding
  double        hardPtCut     = 2.0;                      // hard cut on constituent momentum for initial jet finding
  unsigned      binsEta       = 24;                       // default number of bins for eta for correlation histograms
  unsigned      binsPhi       = 24;                       // default number of bins for phi for correlation histograms
  std::string		outputDir 		= "tmp/";										// directory where everything will be saved
  std::string 	corrOutFile		= "corr.root";							// histograms will be saved here
  std::string		treeOutFile		= "jet.root";								// jets will be saved in a TTree here
  std::string	 	inputFile			= "/nfs/rhi/STAR/Data/CleanAuAuY7/Clean809.root";		// input file: can be .root, .txt, .list
  std::string 	chainName     = "JetTree";								// Tree name in input file
  
  // Now check to see if we were given modifying arguments
  switch ( argc ) {
    case 1: // Default case
      __OUT( "Using Default Settings" )
      break;
    case 16: { // Custom case
      __OUT( "Using Custom Settings" )
      std::vector<std::string> arguments( argv+1, argv+argc );
      // Set non-default values
      // ----------------------
      
      // Choose the analysis type
      if ( arguments[0] == "dijet" ) {
        analysisType 	= "dijet";
        requireDijets = true;
      }
      else if ( arguments[0] == "jet" ) {
        analysisType 	= "jet";
        requireDijets = false;
      }
      else {
        __ERR( "Unknown analysis type: Either dijet or jet " )
        return -1;
      }
      
      // Choose if we use particle-by-particle efficiency
      if ( arguments[1] == "true" ) 			{ useEfficiency = true; }
      else if ( arguments[1] == "false" ) 	{ useEfficiency = false; }
      else { __ERR( "argument 1 must be true or false" ) return -1; }
      
      // Choose if we require a trigger tower to be within
      // jetRadius of the leading jet
      if ( arguments[2] == "true" ) 				{ requireTrigger = true; }
      else if ( arguments[2] == "false" ) 	{ requireTrigger = false; }
      else { __ERR( "argument 2 must be true or false" ) return -1; }
      
      // require there to be a high tower offline
      // above this value
      softwareTrig = atof ( arguments[3].c_str() );
      
      // jet kinematics
      subJetPtMin 	= atof ( arguments[4].c_str() );
      leadJetPtMin 	= atof ( arguments[5].c_str() );
      jetPtMax 			= atof ( arguments[6].c_str() );
      jetRadius 		= atof ( arguments[7].c_str() );
      hardPtCut     = atof ( arguments[8].c_str() );
      
      // bins in eta and phi
      binsEta = atoi ( arguments[9].c_str() );
      binsPhi = atoi ( arguments[10].c_str() );
      
      // output and file names
      outputDir 		= arguments[11];
      corrOutFile		= arguments[12];
      treeOutFile		= arguments[13];
      inputFile 		= arguments[14];
      
      break;
    }
    default: { // Error: invalid custom settings
      __ERR( "Invalid number of command line arguments" )
      return -1;
      break;
    }
  }
  
  
  // Announce our settings
  if ( requireDijets ) { jetHadron::BeginSummaryDijet ( jetRadius, leadJetPtMin, subJetPtMin, jetPtMax, hardPtCut, jetHadron::trackMinPt, jetHadron::binsVz, jetHadron::vzRange, treeOutFile, corrOutFile ); }
  else { jetHadron::BeginSummaryJet ( jetRadius, leadJetPtMin, jetPtMax, hardPtCut, jetHadron::binsVz, jetHadron::vzRange, treeOutFile, corrOutFile ); }
  
  // We know what analysis we are doing now, so build our output histograms
  jetHadron::histograms* histograms = new jetHadron::histograms( analysisType, binsEta, binsPhi );
  histograms->Init();
  
  std::cout<<"histogram analysis type: "<<histograms->GetAnalysisType()<<std::endl;
  
  // Build our input now
  TChain* chain = new TChain( chainName.c_str() );
  // Check to see if the input is a .root file or a .txt
  bool inputIsRoot = jetHadron::HasEnding( inputFile.c_str(), ".root" );
  bool inputIsTxt  = jetHadron::HasEnding( inputFile.c_str(), ".txt"  );
  bool inputIsList = jetHadron::HasEnding( inputFile.c_str(), ".list" );
  
  // If its a recognized file type, build the chain
  // If its not recognized, exit
  if ( inputIsRoot )		 	{ chain->Add( inputFile.c_str() ); }
  else if ( inputIsTxt )  { chain = TStarJetPicoUtils::BuildChainFromFileList( inputFile.c_str() ); }
  else if ( inputIsList )  { chain = TStarJetPicoUtils::BuildChainFromFileList( inputFile.c_str() ); }
  else 										{ __ERR("data file is not recognized type: .root, .list or .txt only.") return -1; }
  
  // Intialize the reader and set the chain
  // All analysis parameters are located in
  // corrParameters.hh
  // --------------------------------------
  TStarJetPicoReader reader;
  jetHadron::InitReader( reader, chain, "auau", jetHadron::triggerHT, softwareTrig, jetHadron::allEvents );
  
  // Data classes
  TStarJetVectorContainer<TStarJetVector>* container;
  TStarJetVector* sv; // TLorentzVector* would be sufficient
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  TClonesArray* triggerObjs;
  
  
  // Build fastjet selectors, containers and definitions
  // ---------------------------------------------------
  
  // Particle container
  std::vector<fastjet::PseudoJet> particles;
  // Trigger container - used to match
  // leading jet with trigger particle
  std::vector<fastjet::PseudoJet> triggers;
  
  // clustering definitions
  // First: used for the analysis - anti-kt with radius jetRadius
  fastjet::JetDefinition 	analysisDefinition = jetHadron::AnalysisJetDefinition( jetRadius );
  // Second: background estimation - kt with radius jetRadius
  fastjet::JetDefinition	backgroundDefinition = jetHadron::BackgroundJetDefinition( jetRadius );
  
  // Build Selectors for the jet finding
  // -----------------------------------
  // Constituent selectors
  fastjet::Selector selectorLowPtCons  = jetHadron::SelectLowPtConstituents( jetHadron::maxTrackRap, jetHadron::trackMinPt );
  fastjet::Selector selectorHighPtCons = jetHadron::SelectHighPtConstituents( jetHadron::maxTrackRap, hardPtCut );
  
  // Jet candidate selector
  fastjet::Selector	selectorJetCandidate;
  if ( requireDijets )
    selectorJetCandidate = jetHadron::SelectJetCandidates( jetHadron::maxTrackRap, jetRadius, subJetPtMin, jetPtMax );
  else
    selectorJetCandidate = jetHadron::SelectJetCandidates( jetHadron::maxTrackRap, jetRadius, leadJetPtMin, jetPtMax );
  
  // Create the Area definition used for background estimation
  fastjet::GhostedAreaSpec	areaSpec = jetHadron::GhostedArea( jetHadron::maxTrackRap, jetRadius );
  fastjet::AreaDefinition 	areaDef  = jetHadron::AreaDefinition( areaSpec );
  
  // selector used to reject hard jets in background estimation
  fastjet::Selector	selectorBkgEstimator	= jetHadron::SelectBkgEstimator( jetHadron::maxTrackRap, jetRadius );
  
  
  // When we do event mixing we need the jets, so save them
  // in trees
  // Tree to hold the di-jets
  TTree* correlatedDiJets;
  // To hold the TLorentzVectors for leading, subleading
  TLorentzVector leadingJet, subleadingJet;
  // Records centrality and vertex information for event mixing
  Int_t centralityBin, vertexZBin;
  Double_t dijetAj = 1.0;
  // Branchs to be written to file
  TBranch* CDJBranchHi, * CDJBranchLo;
  TBranch* CDJBranchCentralityBin;
  TBranch* CDJBranchVertexZBin;
  TBranch* CDJBranchAj;
  
  if ( requireDijets ) {
    correlatedDiJets = new TTree("dijets","Correlated Dijets" );
    CDJBranchVertexZBin = correlatedDiJets->Branch("vertexZBin", &vertexZBin );
    CDJBranchCentralityBin = correlatedDiJets->Branch("centralityBin", &centralityBin );
    CDJBranchHi = correlatedDiJets->Branch("leadJet", &leadingJet );
    CDJBranchLo = correlatedDiJets->Branch("subLeadJet", &subleadingJet );
    CDJBranchAj = correlatedDiJets->Branch("aj", &dijetAj );
  }
  else {
    correlatedDiJets = new TTree("jets","Correlated Jets" );
    CDJBranchHi = correlatedDiJets->Branch("triggerJet", &leadingJet );
    CDJBranchVertexZBin = correlatedDiJets->Branch("vertexZBin", &vertexZBin );
    CDJBranchCentralityBin = correlatedDiJets->Branch("centralityBin", &centralityBin );
  }
  
  // Finally, make ktEfficiency obj for pt-eta
  // Efficiency corrections
  ktTrackEff efficiencyCorrection( jetHadron::y7EfficiencyFile );
  
  // Now everything is set up
  // We can start the event loop
  // First, our counters
  int nEvents = 0;
  int nHardDijets = 0;
  int nMatchedHard = 0;
  
  try{
    while ( reader.NextEvent() ) {
      
      // Count the event
      nEvents++;
      
      // Print out reader status every 10 seconds
      reader.PrintStatus(10);
      
      // Get the event header and event
      event = reader.GetEvent();
      header = event->GetHeader();
      
      // Get the output container from the reader
      container = reader.GetOutputContainer();
      // and trigger objects
      triggerObjs = event->GetTrigObjs();
      
      // Find the reference centrality
      // for y14 it takes the corrected gRefMult and
      // corresponding reference centrality
      int gRefMult = 0;
      int refCent  = 0;
      if ( header->GetCorrectedGReferenceMultiplicity() ) {
        gRefMult = header->GetCorrectedGReferenceMultiplicity();
        refCent = header->GetGReferenceCentrality();
      }
      else {
        gRefMult = header->GetGReferenceMultiplicity();
        refCent  = jetHadron::GetReferenceCentrality( gRefMult );
      }
      // Define the opposite centrality index: 0->8, 1->7, 2->6...
      // Used for the histogram arrays, etc
      int refCentAlt = jetHadron::GetReferenceCentralityAlt( refCent );
      
      // Find vertex Z bin
      double vertexZ = header->GetPrimaryVertexZ();
      int VzBin = jetHadron::GetVzBin( vertexZ );
      
      // Check to see if we use those centralities
      if ( refCent < 0 )                      							 { continue; }
      if ( refCent < jetHadron::y7EfficiencyRefCentLower )   { continue; }
      if ( refCent > jetHadron::y7EfficiencyRefCentUpper )   { continue; }
      
      // Check to see if Vz is in the accepted range; if not, discard
      if ( VzBin == -1 )																				{ continue; }
      
      // Convert TStarJetVector to PseudoJet
      jetHadron::ConvertTStarJetVector( container, particles, true, 1 );
      
      // Get HT triggers
      jetHadron::GetTriggers( requireTrigger, triggerObjs, triggers );
      
      // If we require a trigger and we didnt find one, then discard the event
      if ( requireTrigger && triggers.size() == 0 ) 						{ continue; }
      
      // Start FastJet analysis
      // ----------------------
      
      // get our two sets of particles:
      // low: |eta| < maxTrackRap && pt > 0.2 GeV Used second when we've found viable hard dijets
      // high: |eta| < maxTrackRap && pt > 2.0 GeV Used first to find hard jets
      std::vector<fastjet::PseudoJet> lowPtCons = selectorLowPtCons( particles );
      std::vector<fastjet::PseudoJet> highPtCons = selectorHighPtCons( particles );
      
      // Find high constituent pT jets
      // NO background subtraction
      // -----------------------------
      // First cluster
      fastjet::ClusterSequence clusterSequenceHigh ( highPtCons, analysisDefinition );
      // Now first apply global jet selector to inclusive jets, then sort by pt
      std::vector<fastjet::PseudoJet> HiResult = fastjet::sorted_by_pt( selectorJetCandidate ( clusterSequenceHigh.inclusive_jets() ) );
      
      // Check to see if there are enough jets,
      // and if they meet the momentum cuts - if dijet, checks if they are back to back
      if ( !jetHadron::CheckHardCandidateJets( analysisType, HiResult, leadJetPtMin, subJetPtMin ) ) 	{ continue; }
      
      // count "dijets" ( monojet if doing jet analysis )
      nHardDijets++;
      
      // make our hard dijet vector
      std::vector<fastjet::PseudoJet> hardJets = jetHadron::BuildHardJets( analysisType, HiResult );
      
      // now recluster with all particles if necessary ( only used for dijet analysis )
      // Find corresponding jets with soft constituents
      // ----------------------------------------------
      std::vector<fastjet::PseudoJet> LoResult;
      if ( requireDijets ) {
        fastjet::ClusterSequenceArea ClusterSequenceLow ( lowPtCons, backgroundDefinition, areaDef ); // WITH background subtraction
        
        // Background initialization
        // -------------------------
        
        // Energy density estimate from median ( pt_i / area_i )
        fastjet::JetMedianBackgroundEstimator bkgdEstimator ( selectorBkgEstimator, backgroundDefinition, areaDef );
        bkgdEstimator.set_particles( lowPtCons );
        // Subtract A*rho from the original pT
        fastjet::Subtractor bkgdSubtractor ( &bkgdEstimator );
        LoResult = fastjet::sorted_by_pt( bkgdSubtractor( ClusterSequenceLow.inclusive_jets() ) );
      }
      
      // Get the jets used for correlations
      // Returns hardJets if doing jet analysis
      // it will match to triggers if necessary - if so, trigger jet is at index 0
      std::vector<fastjet::PseudoJet> analysisJets = jetHadron::BuildMatchedJets( analysisType, hardJets, LoResult, requireTrigger, triggers, jetRadius );
      
      // if zero jets were returned, exit out
      if ( analysisJets.size() == 0 )		{ continue; }
      nMatchedHard++;
      
      // now we have analysis jets, write the trees
      // for future event mixing
      vertexZBin = VzBin;
      centralityBin = refCent;
      if ( requireDijets ) {
        // leading jet
        leadingJet.SetPtEtaPhiE( analysisJets.at(0).pt(), analysisJets.at(0).eta(), analysisJets.at(0).phi_std(), analysisJets.at(0).E() );
        subleadingJet.SetPtEtaPhiE( analysisJets.at(1).pt(), analysisJets.at(1).eta(), analysisJets.at(1).phi_std(), analysisJets.at(1).E() );
        dijetAj = jetHadron::CalcAj( hardJets );
      }
      else {
        leadingJet.SetPtEtaPhiE( analysisJets.at(0).pt(), analysisJets.at(0).eta(), analysisJets.at(0).phi_std(), analysisJets.at(0).E() );
        // set a dummy value for event counting
        dijetAj = 0.05;
      }
      
      // now write
      correlatedDiJets->Fill();
      
      // Now we can fill our event histograms
      histograms->CountEvent( VzBin, refCent, dijetAj );
      histograms->FillGRefMult( gRefMult );
      histograms->FillVz( vertexZ );
      if ( requireDijets ) {
        histograms->FillAjHigh( jetHadron::CalcAj( hardJets ) );
        histograms->FillAjLow( jetHadron::CalcAj( analysisJets ) );
        histograms->FillAjDif( jetHadron::CalcAj( hardJets ), jetHadron::CalcAj( analysisJets ) );
        histograms->FillAjStruct( jetHadron::CalcAj( analysisJets ), analysisJets[0].pt(), refCent );
        histograms->FillLeadJetPt( analysisJets.at(0).pt() );
        histograms->FillLeadJetPtHard( hardJets.at(0).pt() );
        histograms->FillLeadEtaPhi( analysisJets.at(0).eta(), analysisJets.at(0).phi_std() );
        histograms->FillSubJetPt( analysisJets.at(1).pt() );
        histograms->FillSubJetPtHard( hardJets.at(1).pt() );
        histograms->FillSubEtaPhi( analysisJets.at(1).eta(), analysisJets.at(1).phi_std() );
      }
      else {
        histograms->FillJetPt( analysisJets.at(0).pt() );
        histograms->FillJetPtHard( hardJets.at(0).pt() );
        histograms->FillJetEtaPhi( analysisJets.at(0).eta(), analysisJets.at(0).phi_std() );
      }
      
      // Now we can perform the correlations
      for ( int i = 0; i < particles.size(); ++i ) {
        fastjet::PseudoJet assocParticle = particles.at(i);
        
        // if we're using particle - by - particle efficiencies, get it,
        // else, set to one
        double assocEfficiency = 1.0;
        if ( useEfficiency ) { assocEfficiency = efficiencyCorrection.EffAAY07( assocParticle.eta(), assocParticle.pt(), refCentAlt );
        }
        
        // now correlate it with leading and subleading jets
        if ( requireDijets ) {
          jetHadron::correlateLeading( analysisType, VzBin, refCent, histograms, analysisJets.at(0), assocParticle, assocEfficiency, dijetAj );
          jetHadron::correlateSubleading( analysisType, VzBin, refCent, histograms, analysisJets.at(1), assocParticle, assocEfficiency, dijetAj );
        }
        else {
          jetHadron::correlateTrigger( analysisType, VzBin, refCent, histograms, analysisJets.at(0), assocParticle, assocEfficiency );
        }
        
      }
      
      
    }
  }catch ( std::exception& e) {
    std::cerr << "Caught " << e.what() << std::endl;
    return -1;
  }
  
  if ( requireDijets )
    jetHadron::EndSummaryDijet ( nEvents, nHardDijets, nMatchedHard, TimeKeeper.RealTime() );
  else
    jetHadron::EndSummaryJet ( nEvents, nHardDijets, TimeKeeper.RealTime() );
  
  // write out the dijet/jet trees
  TFile*  treeOut   = new TFile( (outputDir + treeOutFile).c_str(), "RECREATE" );
  treeOut->cd();
  correlatedDiJets->Write();
  treeOut->Close();
  
  // write out the histograms
  TFile* histOut = new TFile( (outputDir + corrOutFile).c_str(), "RECREATE");
  histOut->cd();
  histograms->Write();
  histOut->Close();
  
  return 0;
}

