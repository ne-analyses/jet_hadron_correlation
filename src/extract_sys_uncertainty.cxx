// macro to extract systematic uncertainties
// for pp correlations wrt auau
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
#include <random>
#include <chrono>


int main () {
  // we're going to extract the systematic errors for tower energy
  // scale and tracking efficiency
  
  // first, get systematic files
  
  // files and naming
  std::vector<TFile*> towFiles;
  std::vector<TFile*> trkFiles;
  std::vector<TFile*> mixFiles;
  std::vector<std::string> analysisNames;
  
  // event mixing file
  mixFiles.push_back( new TFile("out/tmp/pp_mix_6.root", "READ") );
  
  TString location = "out/tmp/pp/sys/";
  TString trigger = "trg5.6/"
  
  TString path = location + trigger;
  
  TString towLow = path + "trk_0_tow_-1.root";
  TString towHigh = path + "trk_0_tow_1.root";
  TString trkLow = path + "trk_-1_tow_0.root";
  TString trkHigh = path + "trk_1_tow_0.root";
  
  towFiles.push_back( new TFile( towLow, "READ" ) );
  towFiles.push_back( new TFile( towHigh, "READ" ) );
  trkFiles.push_back( new TFile( trkLow, "READ" ) );
  trkFiles.push_back( new TFile( trkHigh, "READ" ) );
  // Build our bin selector with default settings
  jetHadron::binSelector selector;
  
  // output file
  TFile* out = new TFile("out/tmp/pp/sys/sys.root","UPDATE");
  
  // read in files
  // *************
  
  // Build our initial histogram holders
  std::vector<TH3F*> nEventsTow;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > towCorrIn;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > towCorrInSub;
  
  std::vector<TH3F*> nEventsTrk;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > trkCorrIn;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > trkCorrInSub;
  
  // Build our initial histogram holders
  std::vector<TH3F*> nEventsMix;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > mixIn;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > mixInSub;
  
  // reading in the histograms
  jetHadron::ReadInFiles( towFiles, towCorrIn, towCorrInSub, nEventsTow, selector );
  jetHadron::ReadInFiles( trkFiles, trkCorrIn, trkCorrInSub, nEventsTrk, selector );
  jetHadron::ReadInFilesMix( mixFiles, mixIn, mixInSub, nEventsMix, selector );
  
  
  return 0;
}
