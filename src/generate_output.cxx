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
  
  switch ( argc ) {
    case 1: { // Default case
      __OUT( "Using Default Settings" )
      
      // default files
      TFile* tmp =  = new TFile( "out/tmp/dijet_corr.root", "READ" );
      TFile* tmpMix = new TFile( "out/tmp/dijet_mix.root", "READ" );
      corrFiles.push_back( tmp );
      mixFiles.push_back( tmpMix );

      analysisNames = defaultCorrNames;
      
      break;
    }
    default: {
      if ( (argc-1)%3 != 0 ) {
        __ERR("Need correlation file, mixing file, and analysis name for each entry")
        return -1;
      }
      std::vector<std::string> arguments( argv+1, argv+argc );
      
      // number of correlations
      int nCorrFiles = ( argc - 1 )/3;
      
      analysisNames.resize( nCorrFiles );
      
      for ( int i = 0; i < nCorrFiles; ++i ) {
        
        TFile* tmp = new TFile( arguments[ 3*i ].c_str(), "READ" );
        TFile* tmpMix =  new TFile( arguments[ (3*i)+1 ].c_str(), "READ" );
        
        corrFiles.push_back( tmp );
        mixFiles.push_back( tmp );
        analysisNames[i] = arguments[ (3*i)+2 ];
      }
    }
  }
  
  int nFiles = analysisNames.size();
  
  // Build our bin selector with default settings
  jetHadron::binSelector selector;
  
  // Build our initial histogram holders
  std::vector<TH3F*> nEvents;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leadingCorrelation;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > subleadingCorrelation;
  // and event mixing
  std::vector<TH3F*> nEventsMixing;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leadingMix;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > subleadingMix;
  // reading in the histograms
  jetHadron::ReadInFiles( corrFiles, leadingCorrelation, subleadingCorrelation, nEvents, selector );
  jetHadron::ReadInFiles( mixFiles, leadingMix, subleadingMix, nEventsMixing, selector );
  
  // Find the pt bin center for future use
  std::vector<TH1F*> ptSpectra;
  std::vector<std::vector<double> > ptBinCenters = jetHadron::FindPtBinCenter( leadingCorrelation, ptSpectra, selector );
  
  // Now build the correlation histograms, both the total and the aj split
  
  
  
  return 0;
}


