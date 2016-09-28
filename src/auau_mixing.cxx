// Produces AuAu dijet-hadron
// Or monojet-hadron mixed event correlations
// Useage defined in submit/auau_mixing.csh
// Nick Elsey - 08.29.2016

// The majority of the jetfinding
// And correlation code is located in
// corrFunctions.hh
#include "corrFunctions.hh"

// All reader and histogram settings
// Are located in corrParameters.hh
#include "corrParameters.hh"

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
// [0]: Input directory ( this should be the file name that includes all
//      the correlation parameters
// [1]:	Input tree file name
// [2]: Output file name
// [3]: Is mixing data min bias or HT? ( MB/mb or HT/ht )
// [4]: Mixing data list  ( should pass a list of all root files to be used for mixed events )
// [5]: Total number of events to consider in mixing data set
// [6]: Number of events to mix with each trigger
// [7]: the mixing data list

// DEF MAIN()
int main ( int argc, const char** argv) {
  
  // First check to make sure we're located properly
  std::string currentDirectory = corrAnalysis::getPWD( );
  
  // If we arent in the analysis directory, exit
  if ( !corrAnalysis::HasEnding ( currentDirectory, "jet_hadron_corr" ) ) {
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
 
  // analysis type - this is set later by GetVarsFromDir()
  std::string    analysisType  = "mix";
  bool           requireDijets = false;
  
  // Defaults
  // --------
  // placeholder (first cmd line argument)
  std::string    executable    = "./bin/auau_mixing";
  // home directory for that analysis
  std::string    inputDir      = "out/dijet/dijet_trigger_true_eff_false_lead_20.0_sub_10.0_max_100.0_rad_0.4";
  // file holding the jet tree
  std::string    inputTreeFile = "tree/tree_Clean809.root";
  // output file name
  std::string    outputFile    = "mixing.root";
  // is mixing data min bias or high tower
  std::string    mixDataType   = "MB" ;
  // is mixing data MB? (set by mixDataType)
  bool           isMixMB       = true;
  // number of events in mixing data to consider
  // choose -1 for all
  int            nMixTotal     = -1;
  // number of events to mix with each trigger
  unsigned       nEventsToMix  = 20;
  // the file/list of events to use in event mixing
  std::string    mixEventsFile = "test/mbAuAu.root";
  // Tree name in input file
  std::string 	 chainName     = "JetTree";
  
  // now check if we'll use the defaults or not
  switch ( argc ) {
    case 1: // Default case
      __OUT( "Using Default Settings" )
      break;
    case 9: { // Custom case
      __OUT( "Using Custom Settings" )
      std::vector<std::string> arguments( argv+1, argv+argc );
      // Set non-default values
      // ----------------------
      
      // load file names and directories
      inputDir = arguments[0];
      inputTreeFile = arguments[1];
      outputFile = arguments[2];
      
      // will it be mixed with MB or HT data
      if ( arguments[3] == "HT" ) {
        mixDataType = "HT";
        isMixMB = false;
      }
      else if ( arguments[3] == "MB" ) {
        mixDataType = "MB";
        isMixMB = true;
      }
      else {
        __ERR( "Unknown data type: Either MB or HT " )
        return -1;
      }

      // number of events to use during mixing
      nMixTotal = atoi ( arguments[4].c_str() );
      
      // number of events to mix with each trigger
      nEventsToMix = atoi ( arguments[5].c_str() );
      
      // the file list/ root file for mixing
      mixEventsFile = arguments[6];
      
      break;
    }
    default: { // Error: invalid custom settings
      __ERR( "Invalid number of command line arguments" )
      return -1;
      break;
    }
  }
  
  // Now we'll get the analysis variables from the directory name
  // ( if its in the proper format )
  std::string analysisString = corrAnalysis::GetDirFromPath( inputDir );

  // set kinematic variables the analysis
  // ----------------------------------
  double leadJetPtMin, subJetPtMin, jetPtMax, jetRadius;
  leadJetPtMin = subJetPtMin = jetPtMax = jetRadius = -999;
  bool useEfficiency, matchTrigger;
  
  // set the variables
  int extractResult = corrAnalysis::GetVarsFromString( analysisType, analysisString, leadJetPtMin, subJetPtMin, jetPtMax, jetRadius, useEfficiency, matchTrigger );
  if (  extractResult == 1 )
    __OUT("Successfully parsed analysis variables")
  else if ( extractResult == 0 )
    __OUT("Using default mixing variables")
  else if ( extractResult == -1 ) {
    __ERR("Submitting incorrect trees for the specified analysis: exit")
    return -1;
  }
  else {
    __ERR("Could not process string: exit")
    return -1;
  }
  
  // check to make sure the analysis type is valid
  // if it is, go ahead and make histogram class
  // else, exit
  if ( analysisType == "dijetmix" || analysisType == "jetmix" )
    __OUT("Performing AuAu mixing")
  else if ( analysisType == "ppdijetmix" || analysisType == "ppjetmix" )
    __OUT("Performing PP mixing")
  else {
      __ERR("unknown analysis type while parsing correlation variables: exiting")
      return -1;
    }
  
  // we need to pick a minimum jet pt in case
  // we use HT events
  double mixingJetPtMax = corrAnalysis::GetMixEventJetPtMin( isMixMB, analysisType, leadJetPtMin );
  
  // Initialize the chain
  // Build our input now
  TChain* chain = new TChain( chainName.c_str() );
  
  // check to see if the input is in .root or .list
  // Check to see if the input is a .root file or a .txt
  bool inputIsRoot = corrAnalysis::HasEnding( mixEventsFile.c_str(), ".root" );
  bool inputIsTxt  = corrAnalysis::HasEnding( mixEventsFile.c_str(), ".txt"  );
  bool inputIsList = corrAnalysis::HasEnding( mixEventsFile.c_str(), ".list" );
  
  // If its a recognized file type, build the chain
  // If its not recognized, exit
  if ( inputIsRoot )		 	{ chain->Add( mixEventsFile.c_str() ); }
  else if ( inputIsTxt )  { chain = TStarJetPicoUtils::BuildChainFromFileList( mixEventsFile.c_str() ); }
  else if ( inputIsList ) { chain = TStarJetPicoUtils::BuildChainFromFileList( mixEventsFile.c_str() ); }
  else 										{ __ERR("data file is not recognized type: .root, .list, .txt only.") return -1; }
  
  // Now we can initialize the reader for mixing events
  // Intialize the reader and set the chain
  // All analysis parameters are located in
  // corrParameters.hh
  // --------------------------------------
  TStarJetPicoReader reader;
  
  if ( corrAnalysis::BeginsWith( analysisType, "pp"  ) )
    corrAnalysis::InitReader( reader, chain, "pp", corrAnalysis::triggerAll, corrAnalysis::allEvents );
  else
    corrAnalysis::InitReader( reader, chain, "auau", corrAnalysis::triggerAll, corrAnalysis::allEvents );
  
  
  // Data classes
  TStarJetVectorContainer<TStarJetVector>* container;
  TStarJetVector* sv; // TLorentzVector* would be sufficient
  TStarJetPicoEventHeader* header;
  TStarJetPicoEvent* event;
  TClonesArray* triggerObjs;
  
  // Build a structure to hold all the event IDs in each mixing bin
  // Inner vector = centrality bins -
  // Outer vector = vz bins -
  std::vector<std::vector<std::vector<unsigned> > > mixing_events;
  for (int i = 0; i < corrAnalysis::binsVz; ++i ) {
    std::vector<std::vector<unsigned> > tmp_vec;
    for (int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      std::vector<unsigned> tmp_set;
      tmp_vec.push_back(tmp_set);
    }
    mixing_events.push_back(tmp_vec);
  }

  // And we'll count the total number of events
  unsigned useable_events = 0;
  unsigned total_events = 0;
  
  // Now build the fastjet framework for finding jets
  
  // Build fastjet selectors, containers and definitions
  // ---------------------------------------------------
  
  // Particle container
  std::vector<fastjet::PseudoJet> particles;
  
  // clustering definitions
  // First: used for the analysis - anti-kt with radius jetRadius
  fastjet::JetDefinition 	analysisDefinition = corrAnalysis::AnalysisJetDefinition( jetRadius );
  
  // Build Selectors for the jet finding
  // -----------------------------------
  // Constituent selectors
  fastjet::Selector selectorLowPtCons  = corrAnalysis::SelectLowPtConstituents( corrAnalysis::maxTrackRap, corrAnalysis::trackMinPt );
  fastjet::Selector selectorHighPtCons = corrAnalysis::SelectHighPtConstituents( corrAnalysis::maxTrackRap, corrAnalysis::hardTrackMinPt );

  
  return 0;
}
