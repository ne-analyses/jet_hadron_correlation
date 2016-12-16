// macro to produce output of dijet hadron
// correlations with event mixing, efficiency
// correction and pt bin dependence.
// Nick Elsey 10/04/2016

// All reader and histogram settings
// Are located in corrParameters.hh
#include "corrParameters.hh"

// The majority of the jetfinding
// And correlation code is located in
// corrFunctions.hh
#include "corrFunctions.hh"

// The output functionality is
// located here
#include "outputFunctions.hh"

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
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
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
#include "TStyle.h"

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

// list all input files as arguments -
// 0 = corr1
// 1 = mix1
// 2 = analysis1 identifying string
// 3 = corr2
// 4 = mix2
// .......

int main( int argc, const char** argv) {
  
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);

  
  // First check to make sure we're located properly
  std::string currentDirectory = jetHadron::getPWD( );
  
  // If we arent in the analysis directory, exit
  if ( !(jetHadron::HasEnding ( currentDirectory, "jet_hadron_corr" ) || jetHadron::HasEnding ( currentDirectory, "jet_hadron_correlation" )) ) {
    std::cerr << "Error: Need to be in jet_hadron_corr directory" << std::endl;
    return -1;
  }
  
  std::vector<std::string> defaultCorrNames;
  defaultCorrNames.resize(1);
  defaultCorrNames[0] = "Dijet";
  
  // files and naming
  std::vector<TFile*> corrFiles;
  std::vector<TFile*> mixFiles;
  std::vector<std::string> analysisNames;
  
  // bin to split Aj on
  int ajSplitBin = 0;
  // jet radius
  double R = 0.4;
  // and the output location
  std::string outputDirBase;
  
  switch ( argc ) {
    case 1: { // Default case
      __OUT( "Using Default Settings" )
      
      // default files
      TFile* tmp = new TFile( "out/tmp/auau_corr_20_10.root", "READ" );
      TFile* tmpMix = new TFile( "out/tmp/auau_mix_20_10.root", "READ" );
      corrFiles.push_back( tmp );
      mixFiles.push_back( tmpMix );

      ajSplitBin = 5;
      analysisNames = defaultCorrNames;
      outputDirBase = "/output";
      
      break;
    }
    default: {
      if ( (argc-4)%3 != 0 ) {
        __ERR("Need correlation file, mixing file, and analysis name for each entry")
        return -1;
      }
      std::vector<std::string> arguments( argv+1, argv+argc );
      
      // number of correlations
      int nCorrFiles = ( argc - 3 )/3;
      
      analysisNames.resize( nCorrFiles );
      
      ajSplitBin = atoi ( arguments[ 0 ].c_str() );
      outputDirBase = arguments[ 1 ];
      R = atof( arguments[ 2 ].c_str() );
      
      for ( int i = 0; i < nCorrFiles; ++i ) {
        
        TFile* tmp = new TFile( arguments[ 3*i+3 ].c_str(), "READ" );
        TFile* tmpMix =  new TFile( arguments[ (3*i)+4 ].c_str(), "READ" );
        
        corrFiles.push_back( tmp );
        mixFiles.push_back( tmpMix );
        analysisNames[i] = arguments[ (3*i)+5 ];
      }
    }
  }

  int nFiles = analysisNames.size();
  
  // put full path for output directory
  outputDirBase = currentDirectory + "/" + outputDirBase;
  
  // Build our bin selector with default settings
  jetHadron::binSelector selector;
  selector.ChangeRadius( R );
  
  // Build our initial histogram holders
  std::vector<TH3F*> nEvents;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leadingCorrelationIn;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > subleadingCorrelationIn;
  // and event mixing
  std::vector<TH3F*> nEventsMixing;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leadingMixIn;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > subleadingMixIn;
  // reading in the histograms
  jetHadron::ReadInFiles( corrFiles, leadingCorrelationIn, subleadingCorrelationIn, nEvents, selector );
  jetHadron::ReadInFilesMix( mixFiles, leadingMixIn, subleadingMixIn, nEventsMixing, selector );
  
  // Find the pt bin center for future use
  std::vector<TH1F*> ptSpectra;
  std::vector<std::vector<double> > ptBinCenters = jetHadron::FindPtBinCenter( leadingCorrelationIn, ptSpectra, selector );
  
  // Now build the correlation histograms, both the total and the aj split
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > leadingCorrelation;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > subleadingCorrelation;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > correlationAjBalanced;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > correlationAjUnbalanced;
  
  jetHadron::BuildSingleCorrelation( leadingCorrelationIn, leadingCorrelation, selector );
  jetHadron::BuildSingleCorrelation( subleadingCorrelationIn, subleadingCorrelation, selector );
  jetHadron::BuildAjSplitCorrelation( leadingCorrelationIn, correlationAjUnbalanced, correlationAjBalanced, selector, ajSplitBin );
  
  // get averaged correlations
  std::vector<std::vector<TH2F*> > averagedSignal = jetHadron::AverageCorrelations( leadingCorrelation, selector );
  std::vector<std::vector<TH2F*> > averagedSignalBalanced = jetHadron::AverageCorrelations( correlationAjBalanced, selector, "balanced" );
  std::vector<std::vector<TH2F*> > averagedSignalUnbalanced = jetHadron::AverageCorrelations( correlationAjUnbalanced, selector, "unbalanced" );
  
  // Now build and scale the event mixing histograms.
  // we will use the averaged event mixing for now
  // First average over all centrality/vz/aj, project into pt
  std::vector<std::vector<TH2F*> > leadingMix =  jetHadron::RecombineMixedEvents( leadingMixIn, selector, "avg_mix_" );
  std::vector<std::vector<TH2F*> > subleadingMix = jetHadron::RecombineMixedEvents( subleadingMixIn, selector, "avg_mix_sub" );
  
  // Build mixed events that are still not averaged as well
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > leadingMixNotAveraged = jetHadron::BuildMixedEvents( leadingMixIn, selector, "not_avg_mix" );
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > subleadingMixNotAveraged = jetHadron::BuildMixedEvents( subleadingMixIn, selector, "not_avg_mix_sub");
  
  // And Scale so that MaxBinContent = 1
  jetHadron::ScaleMixedEvents( leadingMix );
  jetHadron::ScaleMixedEvents( subleadingMix );
  jetHadron::ScaleMixedEvents( leadingMixNotAveraged );
  jetHadron::ScaleMixedEvents( subleadingMixNotAveraged );
  
  // And get the event mixing corrected both averaged or not
  std::vector<std::vector<TH2F*> > averagedMixedEventCorrected = jetHadron::EventMixingCorrection( leadingCorrelation, leadingMix, selector, "leading_avg"  );
  std::vector<std::vector<TH2F*> > notAveragedMixedEventCorrected = jetHadron::EventMixingCorrection( leadingCorrelation, leadingMixNotAveraged, selector, "leading_not_avg" );
  std::vector<std::vector<TH2F*> > averagedMixedEventCorrectedSub = jetHadron::EventMixingCorrection( subleadingCorrelation, subleadingMix, selector, "subleading_avg"  );
  std::vector<std::vector<TH2F*> > notAveragedMixedEventCorrectedSub = jetHadron::EventMixingCorrection( subleadingCorrelation, subleadingMixNotAveraged, selector, "subleading_not_avg" );
  
  // ***************************
  // print out the 2d histograms
  // ***************************
  for ( int i = 0; i < nFiles; ++ i ) {
    jetHadron::Print2DHistograms( leadingMix[i], outputDirBase+"/mixing_"+analysisNames[i], analysisNames[i], selector );
    jetHadron::Print2DHistograms( averagedSignal[i], outputDirBase+"/uncorr_lead_"+analysisNames[i], analysisNames[i], selector );
    jetHadron::Print2DHistogramsEtaRestricted( averagedMixedEventCorrected[i], outputDirBase+"/avg_mix_corrected_lead_"+analysisNames[i], analysisNames[i], selector );
    jetHadron::Print2DHistogramsEtaRestricted( notAveragedMixedEventCorrected[i], outputDirBase+"/mix_corrected_lead_"+analysisNames[i], analysisNames[i], selector );
  }
  
  // **********************
  // get the 1D projections
  // **********************
  
  std::vector<std::vector<TH1F*> > uncorrected_dphi = jetHadron::ProjectDphi( averagedSignal, selector, "not_mixing_corrected", false );
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted = jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrected, selector, "mixing_corrected_near_far_sub_dphi", true );
  
  // *******************************
  // get 1D projections for Aj split
  // *******************************
  std::vector<std::vector<TH1F*> > aj_balanced_dphi = jetHadron::ProjectDphi( averagedSignalBalanced, selector, "aj_balanced_", false );
  std::vector<std::vector<TH1F*> > aj_unbalanced_dphi = jetHadron::ProjectDphi( averagedSignalUnbalanced, selector, "aj_unbalanced_", false );
  
  jetHadron::Normalize1DAjSplit( aj_balanced_dphi, nEvents, 1, ajSplitBin );
  jetHadron::Normalize1DAjSplit( aj_unbalanced_dphi, nEvents, ajSplitBin+1, 20 );
  
  return 0;
}


