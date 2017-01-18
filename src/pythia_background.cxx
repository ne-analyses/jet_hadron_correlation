// Produces an estimate of the
// background contribution of the underlying event
// in AuAu jetfinding using pythia
// Nick Elsey - 01.17.2017

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


// And we need pythia for embedding
#include "Pythia8/Pythia.h"

// used to convert pythia events to vectors of pseudojets
int convertToPseudoJet( Pythia8::Pythia& p, double max_rap, std::vector<fastjet::PseudoJet>& all ) {
  
  // get partons first
  // the initial protons are ids 1 & 2,
  // the hard scattering products are ids 5 & 6
  
  // now loop over all particles, and fill the vectors
  for ( int i = 0; i < p.event.size(); ++i ) {
    if ( p.event[i].isFinal() && p.event[i].isVisible() ) {
      fastjet::PseudoJet tmp( p.event[i].px(), p.event[i].py(), p.event[i].pz(), p.event[i].e() );
      tmp.set_user_index( p.event[i].charge() );
      
      // check to make sure its in our rapidity range
      if ( fabs( tmp.rap() ) > max_rap  )
        continue;
      
      all.push_back( tmp );
      
    }
  }
  
  return 1;
}



int main( int argc, const char** argv ) {
  
  // First check to make sure we're located properly
  std::string currentDirectory = jetHadron::getPWD( );
  
  // If we arent in the analysis directory, exit
  if ( !(jetHadron::HasEnding ( currentDirectory, "jet_hadron_corr" ) || jetHadron::HasEnding ( currentDirectory, "jet_hadron_correlation" )) ) {
    std::cerr << "Error: Need to be in jet_hadron_corr directory" << std::endl;
    return -1;
  }

  // Histograms will calculate gaussian errors
  // -----------------------------------------
  TH1::SetDefaultSumw2( );
  TH2::SetDefaultSumw2( );
  TH3::SetDefaultSumw2( );

  // Read in command line arguments
  // ------------------------------
  // Defaults
  std::string analysisType    = "dijet";
  std::string 	executable    = "./bin/pythia_background"; // placeholder
  double        softwareTrig  = 0.0;                     // require there to be an offline trigger with E > softwareTrig
  double 				subJetPtMin   = 10.0;											// subleading jet minimum pt requirement
  double 				leadJetPtMin  = 20.0;											// leading jet minimum pt requirement
  double				jetPtMax			= 100.0;										// maximum jet pt
  double				jetRadius 		= 0.4;                      // jet radius for jet finding
  double        hardPtCut     = 2.0;                      // hard cut on constituent momentum for initial jet finding
  unsigned      binsEta       = 22;                       // default number of bins for eta for correlation histograms
  unsigned      binsPhi       = 22;                       // default number of bins for phi for correlation histograms
  unsigned      nEvents       = 10e6;
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
    case 6: { // Custom case
      __OUT( "Using Custom Settings" )
      std::vector<std::string> arguments( argv+1, argv+argc );
      // Set non-default values
      // ----------------------
      
      // output and file names
      nEvents       = atoi ( arguments[0].c_str() );
      outputDir 		= arguments[1];
      corrOutFile		= arguments[2];
      treeOutFile		= arguments[3];
      inputFile 		= arguments[4];
      
      break;
    }
    default: { // Error: invalid custom settings
      __ERR( "Invalid number of command line arguments" )
      return -1;
      break;
    }
  }

  // build our two histograms
  TH3D* hCorrLead = new TH3D( "correlationsLead", "correlationsLead", binsEta, -2, 2, binsPhi, jetHadron::phiLowEdge, jetHadron::phiHighEdge, jetHadron::binsPt, jetHadron::ptLowEdge, jetHadron::ptHighEdge  );
  TH3D* hCorrSub = new TH3D( "correlationsSub", "correlationsSub", binsEta, -2, 2, binsPhi, jetHadron::phiLowEdge, jetHadron::phiHighEdge, jetHadron::binsPt, jetHadron::ptLowEdge, jetHadron::ptHighEdge  );
  TH1D* hCounter = new TH1D("events", "events", 1, -1, 1 );
  TH1D* hLead = new TH1D("leadPt", "leadPt", 50, 0, 50 );
  TH1D* hSub  = new TH1D("subPt", "subPt", 50, 0, 50 );
  TH1D* result = new TH1D("ptcount", "ptcount", jetHadron::binsPt, jetHadron::ptLowEdge, jetHadron::ptHighEdge);
  TH1D* resultSub = new TH1D("ptcountSub", "ptcountSub", jetHadron::binsPt, jetHadron::ptLowEdge, jetHadron::ptHighEdge);
  
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
  jetHadron::InitReader( reader, chain, "auau", jetHadron::triggerAll, softwareTrig, jetHadron::allEvents );
  
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
  // pure auau event
  std::vector<fastjet::PseudoJet> auauParticles;
  
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

    selectorJetCandidate = jetHadron::SelectJetCandidates( jetHadron::maxTrackRap, jetRadius, subJetPtMin, jetPtMax );

  
  // Create the Area definition used for background estimation
  fastjet::GhostedAreaSpec	areaSpec = jetHadron::GhostedArea( jetHadron::maxTrackRap, jetRadius );
  fastjet::AreaDefinition 	areaDef  = jetHadron::AreaDefinition( areaSpec );
  
  // selector used to reject hard jets in background estimation
  fastjet::Selector	selectorBkgEstimator	= jetHadron::SelectBkgEstimator( jetHadron::maxTrackRap, jetRadius );
  
  // Finally, make ktEfficiency obj for pt-eta
  // Efficiency corrections
  ktTrackEff efficiencyCorrection( jetHadron::y7EfficiencyFile );
  
  // finally, make a pythia generator
  Pythia8::Pythia pythia( "/wsu/home/dx/dx54/dx5412/software/pythia8219/share/Pythia8/xmldoc" );
  
  // settings for LHC pp at 13 TeV
  pythia.readString("Beams:eCM = 200");
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  pythia.readString("PhaseSpace:pTHatMin = 14.0");
  
  // initialize the pythia generator
  pythia.init();
  
  
  // Now everything is set up
  // We can start the event loop
  // First, our counters
  int eventCount = 0;
  int nHardDijets = 0;
  int nMatchedHard = 0;

  // Now we can do a loop over pythia events
  while ( eventCount < nEvents ) {
    eventCount++;
    
    pythia.next();
    
    // read in the next AuAu event
    // Print out reader status every 10 seconds
    reader.PrintStatus(10);
    
    // loop the HT events
    // If the reader runs out, start again
    if ( !reader.NextEvent() ) {
      reader.ReadEvent(0);
      std::cout<<"RESET AuAu events"<<std::endl;
    }
    
    // Get the event header and event
    event = reader.GetEvent();
    header = event->GetHeader();
    // Get the output container from the reader
    container = reader.GetOutputContainer();
    
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

    // Check to see if we use those centralities
    if ( refCent < 0 )                      								 	{ continue; }
    if ( refCent < jetHadron::y7EfficiencyRefCentLower )   { continue; }
    if ( refCent > jetHadron::y7EfficiencyRefCentUpper )   { continue; }
    
    // Find vertex Z bin
    double vertexZ = header->GetPrimaryVertexZ();
    int VzBin = jetHadron::GetVzBin( vertexZ );
    
    // Check to see if Vz is in the accepted range; if not, discard
    if ( VzBin == -1 )																				{ continue; }
    
    
    // we will make sure there is no high energy jet in the event
    bool goodBkgEvent = false;
    while ( goodBkgEvent != true ) {
      std::vector<fastjet::PseudoJet> tmpParticles;
      jetHadron::ConvertTStarJetVector( container, tmpParticles, true, 0 );
      
      std::vector<fastjet::PseudoJet> highPtCons = selectorHighPtCons( tmpParticles );
      
      // First cluster
      fastjet::ClusterSequence clusterSequenceHigh ( highPtCons, analysisDefinition );
      // Now first apply global jet selector to inclusive jets, then sort by pt
      std::vector<fastjet::PseudoJet> HiResult = fastjet::sorted_by_pt( selectorJetCandidate ( clusterSequenceHigh.inclusive_jets() ) );
      
      if ( HiResult.size() && HiResult[0].pt() > leadJetPtMin*0.8 ) {
        goodBkgEvent = true;
      }
      else {
        reader.NextEvent();
      }
    }

    
    // Convert TStarJetVector to PseudoJet
    jetHadron::ConvertTStarJetVector( container, particles, true, 0 );
    jetHadron::ConvertTStarJetVector( container, auauParticles, true, 0 );
    convertToPseudoJet( pythia, 1, particles );
    
    // now get the particles we will use
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
    
    // make our hard dijet vector
    std::vector<fastjet::PseudoJet> hardJets = jetHadron::BuildHardJets( analysisType, HiResult );
    
    // ok, we're using it, so count it
    hCounter->Fill(0);
    hLead->Fill( hardJets[0].pt() );
    hSub->Fill( hardJets[1].pt() );
    
    // now loop over all tracks in auau event
    for ( int j = 0; j < auauParticles.size(); ++j ) {
      fastjet::PseudoJet assocTrack = particles.at(j);
      
      if ( assocTrack.user_index() == 0 )
        continue;
      if ( fabs( assocTrack.eta() ) > jetHadron::maxTrackRap )
        continue;
      
      
      // if we're using particle - by - particle efficiencies, get it,
      // else, set to one
      double assocEfficiency = efficiencyCorrection.EffAAY07( assocTrack.eta(), assocTrack.pt(), refCentAlt );
      
      // get our correlations
      
      double deltaEta = hardJets[0].eta() - assocTrack.eta();
      double deltaPhi = hardJets[0].delta_phi_to( assocTrack );
      double deltaEtaSub = hardJets[1].eta() - assocTrack.eta();
      double deltaPhiSub = hardJets[1].delta_phi_to( assocTrack );
      double assocPt =	assocTrack.pt();
      double weight = 1.0/assocEfficiency;
      
      // aaaaaand plot if its about 2 GeV
      if ( assocPt > 2.0 ) {
        hCorrLead->Fill( deltaEta, deltaPhi, assocPt, weight );
        hCorrSub->Fill( deltaEtaSub, deltaPhiSub, assocPt, weight );
      }
      if ( hardJets[0].delta_R( assocTrack ) )
        result->Fill( assocPt, weight );
      if ( hardJets[1].delta_R( assocTrack ) )
        resultSub->Fill( assocPt, weight );
    }
    
  }
  
  // write out the dijet/jet trees
  TFile*  Out   = new TFile( (outputDir + corrOutFile).c_str(), "RECREATE" );
  Out->cd();
  hCorrLead->Write();
  hCorrSub->Write();
  hLead->Write();
  hSub->Write();
  hCounter->Write();
  result->Write();
  resultSub->Write();
  Out->Close();
  
  return 0;
}






