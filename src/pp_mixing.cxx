// Produces PP dijet-hadron
// Or monojet-hadron mixed event correlations
// Useage defined in submit/auau_correlation.csh
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
// Defined for debugging in main() )
// [0]: Choose between dijet-hadron correlations or jet-hadron correlations
//      choices: dijet || jet
// [1]:	Choose to use particle efficiency corrections or not: true/false
// [2]: trigger coincidence: require a trigger with the leading jet
// [3]: subleading jet pt min ( not used if doing jet-hadron correlations )
// [4]: leading jet pt min
// [5]: jet pt max
// [6]: jet radius, used in the jet definition
// [7]: output directory
// [8]: name for the correlation histogram file
// [9]: name for the dijet TTree file
// [10]: input data: can be a single .root or a .txt or .list of root files
// [11]: MB AuAu event file for embedding, can be .root, .txt, .list

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
  
  
  
  
  
  
  return 0;
}
