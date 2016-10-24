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
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <random>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

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

// the grid does not have std::to_string() for some ungodly reason
// replacing it here. Simply ostringstream
namespace patch {
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}


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
  if ( !(corrAnalysis::HasEnding ( currentDirectory, "jet_hadron_corr" ) || corrAnalysis::HasEnding ( currentDirectory, "jet_hadron_correlation" )) ) {
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
  std::string    inputDir      = "out/dijet/dijet_trigger_true_eff_true_lead_20.0_sub_10.0_max_100.0_rad_0.4";
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
  std::string    mixEventsFile = "auau_list/grid_AuAuy7MB.list";
  // Tree name in input file
  std::string 	 chainName     = "JetTree";
  
  // now check if we'll use the defaults or not
  switch ( argc ) {
    case 1: // Default case
      __OUT( "Using Default Settings" )
      break;
    case 8: { // Custom case
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
  
  if ( analysisType == "dijetmix" || analysisType == "ppdijetmix" )
    requireDijets = true;
  std::cout<<"Use Trigger: "<< matchTrigger<<std::endl;
  std::cout<<"Use Efficiency: "<< useEfficiency<<std::endl;
  
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
  
  // initialize histogram container
  corrAnalysis::histograms* histograms = new corrAnalysis::histograms( analysisType );
  histograms->Init();
  
  // we need to pick a minimum jet pt in case
  // we use HT events
  double mixingJetPtMax = corrAnalysis::GetMixEventJetPtMax( isMixMB, analysisType, leadJetPtMin );
  
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
    corrAnalysis::InitReader( reader, chain, "pp", corrAnalysis::triggerAll, nMixTotal );
  else
    corrAnalysis::InitReader( reader, chain, "auau", corrAnalysis::triggerAll, nMixTotal );
  
  
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
  
  // initial jet selector to find if there is a high momentum jet
  // only used if the data is HT triggered
  // looks for jets with pt > 0.8*jetPtMin from the analysis
  fastjet::Selector selectorFindHTJet = corrAnalysis::SelectJetCandidates( corrAnalysis::maxTrackRap, jetRadius, mixingJetPtMax, jetPtMax );
  
  // make ktEfficiency obj for pt-eta
  // Efficiency corrections
  // Only used if efficiency corrections were
  // used in the correlation analysis
  ktTrackEff efficiencyCorrection( corrAnalysis::y7EfficiencyFile );
  
  // Now we need to load the dijet trees
  // -----------------------------------
  // The root file
  TFile treeFile( (inputDir+"/"+inputTreeFile).c_str(), "READ");
  
  // pull the tree
  TTree* jetTree;
  if (analysisType == "dijetmix" )
    jetTree = (TTree*) treeFile.Get( "dijets" );
  else if ( analysisType == "jetmix" )
    jetTree = (TTree*) treeFile.Get( "jets" );
  else if ( analysisType == "ppdijetmix" )
    jetTree = (TTree*) treeFile.Get( "pp_dijets" );
  else if ( analysisType == "ppjetmix" )
    jetTree = (TTree*) treeFile.Get( "pp_jets" );
  else {
    __ERR("analysis type unrecognized")
    return -1;
  }
  
  // Check to make sure the tree exists
  if ( !jetTree ) {
    __ERR("error retrieving tree from input file")
    return -1;
  }
  
  // Get the number of entries
  unsigned treeEntries = jetTree->GetEntries();
  
  if ( treeEntries == 0 ) {
    __ERR("no entries in jet tree - aborting")
    return -1;
  }
  
  std::string reportEntries = "total of " + patch::to_string( treeEntries ) + " trigger events";
  __OUT( reportEntries.c_str() )
  // Define our branches
  TLorentzVector *leadBranch = new TLorentzVector();
  TLorentzVector *subBranch = new TLorentzVector();
  int centBranch, vzBranch;
  
  // set the branch addresses for the tree
  if ( requireDijets ) {
    jetTree->SetBranchAddress( "leadJet", &leadBranch );
    jetTree->SetBranchAddress( "subLeadJet", &subBranch );
    jetTree->SetBranchAddress( "vertexZBin", &vzBranch );
    if ( analysisType == "dijetmix" )
      jetTree->SetBranchAddress( "centralityBin", &centBranch );
    
  }
  else {
    jetTree->SetBranchAddress( "triggerJet", &leadBranch );
    jetTree->SetBranchAddress( "vertexZBin", &vzBranch );
    if ( analysisType == "jetmix" )
      jetTree->SetBranchAddress( "centralityBin", &centBranch );

  }
  __OUT("loaded branches")
  
  
  // collect info on how many events there are to mix with
  TH2D* hCentVz = new TH2D( "cent_vz", "Mixing Event Count;centrality;vz", corrAnalysis::binsCentrality, -0.5, corrAnalysis::binsCentrality-0.5, corrAnalysis::binsVz, -0.5, corrAnalysis::binsVz-0.5 );

  
  // Now loop over events and store their IDs in the proper
  // Vz-Centrality bin
  // And we'll count the total number of events
  unsigned useable_events = 0;
  unsigned total_events = 0;
  try{
    while ( reader.NextEvent() ) {
      // Count the event
      total_events++;
      
      // Print out reader status every 10 seconds
      reader.PrintStatus(10);
      
      // Get the event header and event
      event = reader.GetEvent();
      header = event->GetHeader();
      
      // A few variables needed
      // ----------------------
      
      // Vz position and corresponding bin
      double vertexZ = header->GetPrimaryVertexZ();
      int vzBin = corrAnalysis::GetVzBin( vertexZ );
      
      // Get the centrality information
      // Find the reference centrality
      // (for pp, set to zero by default )
      int gRefMult = header->GetGReferenceMultiplicity();
      int refCentrality = corrAnalysis::GetReferenceCentrality( gRefMult );
      if ( corrAnalysis::BeginsWith( analysisType, "pp") )
        refCentrality = 0;
      
      // Now we need to check if it has a hard jet in it
      // Get the output container from the reader
      container = reader.GetOutputContainer();
      
      // Clear the PseudoJet container
      particles.clear();
      
      // Transform TStarJetVectors into (FastJet) PseudoJets
      // ---------------------------------------------------
      for ( int i=0; i < container->GetEntries() ; ++i ){
        sv = container->Get(i);
        
        fastjet::PseudoJet tmpPJ = fastjet::PseudoJet( *sv );
        tmpPJ.set_user_index( sv->GetCharge() );
        particles.push_back( tmpPJ );
      }
      
      // Apply selector to the full particle set
      // pHi = hard jet constituents
      std::vector<fastjet::PseudoJet> pHi = selectorHighPtCons( particles );
      
      // Find high constituent pT jets
      // NO background subtraction
      // -----------------------------
      // First cluster
      fastjet::ClusterSequence csaHi ( pHi, analysisDefinition );
      // Now first apply global jet selector to inclusive jets, then sort by pt
      std::vector<fastjet::PseudoJet> HiResult = fastjet::sorted_by_pt( selectorFindHTJet ( csaHi.inclusive_jets() ) );
      
      // check to see if the event needs to be discarded
      if ( !corrAnalysis::UseEventInMixing( analysisType, isMixMB, HiResult, gRefMult, vzBin ) )
        continue;
     
      // Now, we know its an event we will use, so fill it in
      useable_events++;
      unsigned tmp_event_id = reader.GetNOfCurrentEvent();
      mixing_events[vzBin][refCentrality].push_back( tmp_event_id );
      hCentVz->Fill( refCentrality, vzBin );
      
    }
  }catch ( std::exception& e) {
    std::cerr << "Caught " << e.what() << std::endl;
    return -1;
  }
  __OUT("Finished inital event binning in Vz-centrality")
  std::string finishEventCheck = "Out of " + patch::to_string(total_events) + ", " + patch::to_string(useable_events) + " will be used";
  __OUT( finishEventCheck.c_str() )
  // Quick check for the size of these arrays
  // Remove cent/vz bins that have too few events
  __OUT("Checking each Vz/centrality bin for the minimum number of entries")
  for ( int i = 0; i < corrAnalysis::binsVz; ++i )
    for ( int j = 6; j < corrAnalysis::binsCentrality; ++j ) {
      if ( mixing_events[i][j].size() < nEventsToMix*1.5 ) {
        mixing_events[i][j].clear();
        std::string outMessage = "Removing bin ";
        outMessage += patch::to_string(i);
        outMessage += " ";
        outMessage += patch::to_string(j);
        __OUT( outMessage.c_str() )
      }
    }
  __OUT("Done removing bins")
  
  // create a RNG for shuffling events
  std::random_device rd;
  std::mt19937 g(rd());
  g.seed( clock() );
  
  // Now we can run over all tree entries and perform the mixing
  __OUT("Starting to perform event mixing")
  for ( int i = 0; i < treeEntries; ++i ) {

    // Pull the next jet/dijet
    jetTree->GetEntry(i);
    
    if ( i % 20 == 0) {
      std::string eventOut = "Mixing tree entry: " + patch::to_string(i);
      __OUT( eventOut.c_str() )
    }
    
    // get the proper cent/vz bin
    std::vector< unsigned > randomizedEventID;
    if ( corrAnalysis::BeginsWith( analysisType, "pp") )
      randomizedEventID = mixing_events[vzBranch][0];
    else
      randomizedEventID = mixing_events[vzBranch][centBranch];
    
    // If the event list was set to zero earlier,
    // Then we will not be using that bin
    if ( randomizedEventID.size() == 0 )  { continue;}
    
    // then randomize the list
    std::shuffle( randomizedEventID.begin(), randomizedEventID.end(), g );
    
    // now use the first nEventsToMix
    for ( int i = 0; i < nEventsToMix; ++i ) {
      // get the event
      reader.ReadEvent( randomizedEventID[i] );
      
      // count event
      if ( corrAnalysis::BeginsWith( analysisType, "pp") )
        histograms->CountEvent( vzBranch );
      else
        histograms->CountEvent( centBranch, vzBranch );
      
      // get the reference centrality definition used by
      // the track efficiency class
      int refCentAlt = corrAnalysis::GetReferenceCentralityAlt( centBranch );
      
      // get event headers
      event = reader.GetEvent();
      header = event->GetHeader();
      container = reader.GetOutputContainer();
      
      // first convert to pseudojets
      particles.clear();
      corrAnalysis::ConvertTStarJetVector( container, particles );
      
      // now do the correlation
      if ( requireDijets ) {
        // make the trigger pseudojets
        fastjet::PseudoJet leadTrigger = fastjet::PseudoJet( *leadBranch );
        fastjet::PseudoJet subTrigger = fastjet::PseudoJet( *subBranch );
        
        histograms->FillLeadEtaPhi( leadTrigger.eta(), leadTrigger.phi_std() );
        histograms->FillSubEtaPhi( subTrigger.eta(), subTrigger.phi_std() );
        
        // loop over associated particles
        for ( int j = 0; j < particles.size(); ++j ) {
          fastjet::PseudoJet assocParticle = particles[j];
          
          // if we're using particle - by - particle efficiencies, get it,
          // else, set to one
          double assocEfficiency = 1.0;
          if ( useEfficiency ) { assocEfficiency = efficiencyCorrection.EffAAY07( assocParticle.eta(), assocParticle.pt(), refCentAlt );
          }
          
          corrAnalysis::correlateLeading( analysisType, vzBranch, centBranch, histograms, leadTrigger, assocParticle, assocEfficiency );
          corrAnalysis::correlateSubleading( analysisType, vzBranch, centBranch, histograms, subTrigger, assocParticle, assocEfficiency );
        }
        
      }
      else {
        // make the trigger pseudojets
        fastjet::PseudoJet leadTrigger = fastjet::PseudoJet( *leadBranch );
        
        histograms->FillJetEtaPhi( leadTrigger.eta(), leadTrigger.phi_std() );
        
        // loop over associated particles
        for ( int j = 0; j < particles.size(); ++j ) {
          fastjet::PseudoJet assocParticle = particles[j];
          
          // if we're using particle - by - particle efficiencies, get it,
          // else, set to one
          double assocEfficiency = 1.0;
          if ( useEfficiency ) assocEfficiency = efficiencyCorrection.EffAAY07( assocParticle.eta(), assocParticle.pt(), refCentAlt );
          
          corrAnalysis::correlateTrigger( analysisType, vzBranch, centBranch, histograms, leadTrigger, assocParticle, assocEfficiency );

        }
      }
    }
  }
  
  // create an output file
  TFile out((inputDir+"/"+outputFile).c_str(), "RECREATE");
  
  histograms->Write();

  hCentVz->Write();
  
  out.Close();
  
  return 0;
}
