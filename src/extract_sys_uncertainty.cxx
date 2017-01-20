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
  
  // turning off info printing
  gErrorIgnoreLevel = kInfo;
  
  // first, get systematic files
  
  // files and naming
  std::vector<TFile*> towFiles;
  std::vector<TFile*> trkFiles;
  std::vector<TFile*> mixFiles;
  std::vector<std::string> analysisNames;
  
  TString location = "out/added/pp/";
  TString trigger = "trg6";
  
  TString path = location + trigger;
  
  // event mixing file
  TString mixFileName = path + "/mix.root";
  mixFiles.push_back( new TFile( mixFileName, "READ") );
  
  TString towLow = path + "/sys/tow-1trk0.root";
  TString towHigh = path + "/sys/tow1trk0.root";
  TString trkLow = path + "/sys/tow0trk-1.root";
  TString trkHigh = path + "/sys/tow0trk1.root";
  
  towFiles.push_back( new TFile( towLow, "READ" ) );
  towFiles.push_back( new TFile( towHigh, "READ" ) );
  trkFiles.push_back( new TFile( trkLow, "READ" ) );
  trkFiles.push_back( new TFile( trkHigh, "READ" ) );
  
  // Build our bin selector with default settings
  jetHadron::binSelector selector;
  
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
  jetHadron::ReadInFiles( towFiles, towCorrIn, towCorrInSub, nEventsTow, selector, "tow" );
  jetHadron::ReadInFiles( trkFiles, trkCorrIn, trkCorrInSub, nEventsTrk, selector, "trk" );
  jetHadron::ReadInFilesMix( mixFiles, mixIn, mixInSub, nEventsMix, selector );
  
  // Find the pt bin center for future use
  std::vector<TH1F*> ptSpectra;
  std::vector<std::vector<double> > ptBinCenters = jetHadron::FindPtBinCenter( towCorrIn, ptSpectra, selector );
  
  
  // building a pt bin error
  std::vector<std::vector<double> > zeros;
  zeros.resize( ptBinCenters.size() );
  for ( int i = 0; i < ptBinCenters.size(); ++i ) {
    zeros[i].resize( ptBinCenters[i].size() );
  }
  
  // setup the event mixing histograms
  std::vector<std::vector<TH2F*> > leadingMix =  jetHadron::RecombineMixedEvents( mixIn, selector, "avg_mix_" );
  std::vector<std::vector<TH2F*> > subleadingMix = jetHadron::RecombineMixedEvents( mixInSub, selector, "avg_mix_sub" );
  
  //scale the event mixing histograms
  jetHadron::ScaleMixedEvents( leadingMix );
  jetHadron::ScaleMixedEvents( subleadingMix );
  
  // since we might be using the same event mixing stuff...
  // hack it and just do a pointer copy if the second doesnt exist,
  // we can use the same mixed events
  if ( leadingMix.size() == 1 )
    leadingMix.push_back( leadingMix[0] );
  if ( subleadingMix.size() == 1 )
    subleadingMix.push_back( subleadingMix[0] );
  
  // now reorganize the correlations... get rid of the aj split
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > towCorr;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > towCorrSub;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > trkCorr;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > trkCorrSub;

  
  jetHadron::BuildSingleCorrelation( towCorrIn, towCorr, selector, "leadTow" );
  jetHadron::BuildSingleCorrelation( towCorrInSub, towCorrSub, selector, "subTow" );
  jetHadron::BuildSingleCorrelation( trkCorrIn, trkCorr, selector, "leadTrk" );
  jetHadron::BuildSingleCorrelation( trkCorrInSub, trkCorrSub, selector, "subTrk" );
  
  // now get the corrected histograms
  std::vector<std::vector<TH2F*> > correctedTow = jetHadron::EventMixingCorrection( towCorr, leadingMix, selector, "corTow"  );
  std::vector<std::vector<TH2F*> > correctedTowSub = jetHadron::EventMixingCorrection( towCorrSub, subleadingMix, selector, "corTowSub" );
  std::vector<std::vector<TH2F*> > correctedTrk = jetHadron::EventMixingCorrection( trkCorr, leadingMix, selector, "corTrk"  );
  std::vector<std::vector<TH2F*> > correctedTrkSub = jetHadron::EventMixingCorrection( trkCorrSub, subleadingMix, selector, "corTrkSub" );
  
  // get the projections
  // first Subtracted DPhi
  // *************************************

  
  std::vector<std::vector<TH1F*> > corrected_dphi_tow = jetHadron::ProjectDphiNearMinusFar( correctedTow, selector, "corTowDPhi" );
  std::vector<std::vector<TH1F*> > corrected_dphi_tow_sub = jetHadron::ProjectDphiNearMinusFar( correctedTowSub, selector,  "corTowDPhiSub" );
  std::vector<std::vector<TH1F*> > corrected_dphi_trk = jetHadron::ProjectDphiNearMinusFar( correctedTrk, selector,  "corTrkDPhi" );
  std::vector<std::vector<TH1F*> > corrected_dphi_trk_sub = jetHadron::ProjectDphiNearMinusFar( correctedTrkSub, selector, "corTrkDPhiSub"  );
  
  // normalize with 1/dijets 1/bin width
  jetHadron::Normalize1D( corrected_dphi_tow, nEventsTow );
  jetHadron::Normalize1D( corrected_dphi_tow_sub, nEventsTow );
  jetHadron::Normalize1D( corrected_dphi_trk, nEventsTrk );
  jetHadron::Normalize1D( corrected_dphi_trk_sub, nEventsTrk );
  
  // do final fitting
  std::vector<std::vector<TF1*> > corrected_dphi_tow_fit = jetHadron::FitDphi( corrected_dphi_tow, selector );
  std::vector<std::vector<TF1*> > corrected_dphi_tow_sub_fit = jetHadron::FitDphi( corrected_dphi_tow_sub, selector );
  std::vector<std::vector<TF1*> > corrected_dphi_trk_fit = jetHadron::FitDphi( corrected_dphi_trk, selector );
  std::vector<std::vector<TF1*> > corrected_dphi_trk_sub_fit = jetHadron::FitDphi( corrected_dphi_trk_sub, selector );

  // extract fit values
  std::vector<std::vector<double> > corrected_dphi_tow_fit_yield, corrected_dphi_tow_fit_width, corrected_dphi_tow_fit_width_err, corrected_dphi_tow_fit_yield_err;
  std::vector<std::vector<double> > corrected_dphi_tow_sub_fit_yield, corrected_dphi_tow_sub_fit_width, corrected_dphi_tow_sub_fit_width_err, corrected_dphi_tow_sub_fit_yield_err;
  std::vector<std::vector<double> > corrected_dphi_trk_fit_yield, corrected_dphi_trk_fit_width, corrected_dphi_trk_fit_width_err, corrected_dphi_trk_fit_yield_err;
  std::vector<std::vector<double> > corrected_dphi_trk_sub_fit_yield, corrected_dphi_trk_sub_fit_width, corrected_dphi_trk_sub_fit_width_err, corrected_dphi_trk_sub_fit_yield_err;
  
  jetHadron::ExtractFitVals( corrected_dphi_tow_fit, corrected_dphi_tow_fit_yield, corrected_dphi_tow_fit_width, corrected_dphi_tow_fit_yield_err, corrected_dphi_tow_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_dphi_tow_sub_fit, corrected_dphi_tow_sub_fit_yield, corrected_dphi_tow_sub_fit_width, corrected_dphi_tow_sub_fit_yield_err, corrected_dphi_tow_sub_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_dphi_trk_fit, corrected_dphi_trk_fit_yield, corrected_dphi_trk_fit_width, corrected_dphi_trk_fit_yield_err, corrected_dphi_trk_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_dphi_trk_sub_fit, corrected_dphi_trk_sub_fit_yield, corrected_dphi_trk_sub_fit_width, corrected_dphi_trk_sub_fit_yield_err, corrected_dphi_trk_sub_fit_width_err, selector  );
  
  // now build the error histograms
  std::vector<TH1F*> dphi_tow_err = jetHadron::BuildSystematicHistogram( corrected_dphi_tow[0], corrected_dphi_tow[1], selector, "dphi_tow_sys_err" );
  std::vector<TH1F*> dphi_trk_err = jetHadron::BuildSystematicHistogram( corrected_dphi_trk[0], corrected_dphi_trk[1], selector, "dphi_trk_sys_err" );
  std::vector<TH1F*> dphi_tow_err_sub = jetHadron::BuildSystematicHistogram( corrected_dphi_tow_sub[0], corrected_dphi_tow_sub[1], selector, "dphi_tow_sub_sys_err" );
  std::vector<TH1F*> dphi_trk_err_sub = jetHadron::BuildSystematicHistogram( corrected_dphi_trk_sub[0], corrected_dphi_trk_sub[1], selector, "dphi_trk_sub_sys_err" );

  // fix the damn bins
  jetHadron::FixTheDamnBins( dphi_tow_err );
  jetHadron::FixTheDamnBins( dphi_tow_err_sub );
  jetHadron::FixTheDamnBins( dphi_trk_err );
  jetHadron::FixTheDamnBins( dphi_trk_err_sub );
  
  // GOT HERE OK
  // DEBUG FROM HERE
  
  // now DEta
  // *************************************
  std::vector<std::vector<TH1F*> > corrected_deta_tow = jetHadron::ProjectDeta( correctedTow, selector, "corTowDEta");
  std::vector<std::vector<TH1F*> > corrected_deta_tow_sub = jetHadron::ProjectDeta( correctedTowSub, selector, "corTowDEtaSub" );
  std::vector<std::vector<TH1F*> > corrected_deta_trk = jetHadron::ProjectDeta( correctedTrk, selector, "corTrkDEta" );
  std::vector<std::vector<TH1F*> > corrected_deta_trk_sub = jetHadron::ProjectDeta( correctedTrkSub, selector, "corTrkDEtaSub" );
  
  // do background subtraction
  jetHadron::SubtractBackgroundDeta( corrected_deta_tow, selector );
  jetHadron::SubtractBackgroundDeta( corrected_deta_tow_sub, selector );
  jetHadron::SubtractBackgroundDeta( corrected_deta_trk, selector );
  jetHadron::SubtractBackgroundDeta( corrected_deta_trk_sub, selector );
  
  // normalize with 1/dijets 1/bin width
  jetHadron::Normalize1D( corrected_deta_tow, nEventsTow );
  jetHadron::Normalize1D( corrected_deta_tow_sub, nEventsTow );
  jetHadron::Normalize1D( corrected_deta_trk, nEventsTrk );
  jetHadron::Normalize1D( corrected_deta_trk_sub, nEventsTrk );
  
  // do final fitting
  std::vector<std::vector<TF1*> > corrected_deta_tow_fit = jetHadron::FitDeta( corrected_deta_tow, selector );
  std::vector<std::vector<TF1*> > corrected_deta_tow_sub_fit = jetHadron::FitDeta( corrected_deta_tow_sub, selector );
  std::vector<std::vector<TF1*> > corrected_deta_trk_fit = jetHadron::FitDeta( corrected_deta_trk, selector );
  std::vector<std::vector<TF1*> > corrected_deta_trk_sub_fit = jetHadron::FitDeta( corrected_deta_trk_sub, selector );
  
  // extract fit values
  std::vector<std::vector<double> > corrected_deta_tow_fit_yield, corrected_deta_tow_fit_width, corrected_deta_tow_fit_width_err, corrected_deta_tow_fit_yield_err;
  std::vector<std::vector<double> > corrected_deta_tow_sub_fit_yield, corrected_deta_tow_sub_fit_width, corrected_deta_tow_sub_fit_width_err, corrected_deta_tow_sub_fit_yield_err;
  std::vector<std::vector<double> > corrected_deta_trk_fit_yield, corrected_deta_trk_fit_width, corrected_deta_trk_fit_width_err, corrected_deta_trk_fit_yield_err;
  std::vector<std::vector<double> > corrected_deta_trk_sub_fit_yield, corrected_deta_trk_sub_fit_width, corrected_deta_trk_sub_fit_width_err, corrected_deta_trk_sub_fit_yield_err;
  
  jetHadron::ExtractFitVals( corrected_deta_tow_fit, corrected_deta_tow_fit_yield, corrected_deta_tow_fit_width, corrected_deta_tow_fit_yield_err, corrected_deta_tow_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_deta_tow_sub_fit, corrected_deta_tow_sub_fit_yield, corrected_deta_tow_sub_fit_width, corrected_deta_tow_sub_fit_yield_err, corrected_deta_tow_sub_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_deta_trk_fit, corrected_deta_trk_fit_yield, corrected_deta_trk_fit_width, corrected_deta_trk_fit_yield_err, corrected_deta_trk_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_deta_trk_sub_fit, corrected_deta_trk_sub_fit_yield, corrected_deta_trk_sub_fit_width, corrected_deta_trk_sub_fit_yield_err, corrected_deta_trk_sub_fit_width_err, selector  );
  
  // now build the error histograms
  std::vector<TH1F*> deta_tow_err = jetHadron::BuildSystematicHistogram( corrected_deta_tow[0], corrected_deta_tow[1], selector, "deta_tow_sys_err" );
  std::vector<TH1F*> deta_trk_err = jetHadron::BuildSystematicHistogram( corrected_deta_trk[0], corrected_deta_trk[1], selector, "deta_trk_sys_err" );
  std::vector<TH1F*> deta_tow_err_sub = jetHadron::BuildSystematicHistogram( corrected_deta_tow_sub[0], corrected_deta_tow_sub[1], selector, "deta_tow_sub_sys_err" );
  std::vector<TH1F*> deta_trk_err_sub = jetHadron::BuildSystematicHistogram( corrected_deta_trk_sub[0], corrected_deta_trk_sub[1], selector, "deta_trk_sub_sys_err" );
  
  // fix the damn bins
  jetHadron::FixTheDamnBins( dphi_tow_err );
  jetHadron::FixTheDamnBins( dphi_tow_err_sub );
  jetHadron::FixTheDamnBins( dphi_trk_err );
  jetHadron::FixTheDamnBins( dphi_trk_err_sub );
  
  // aaaand build the sum in quadrature
  std::vector<TH1F*> deta_err = jetHadron::AddInQuadrature( deta_tow_err, deta_trk_err, selector, "deta" );
  std::vector<TH1F*> deta_err_sub = jetHadron::AddInQuadrature( deta_tow_err_sub, deta_trk_err_sub, selector, "sub_deta");
  std::vector<TH1F*> dphi_err = jetHadron::AddInQuadrature( dphi_tow_err, dphi_trk_err, selector, "dphi");
  std::vector<TH1F*> dphi_err_sub = jetHadron::AddInQuadrature( dphi_tow_err_sub, dphi_trk_err_sub, selector, "sub_dphi");
  
  // fix the damn bins
  jetHadron::FixTheDamnBins( dphi_err );
  jetHadron::FixTheDamnBins( dphi_err_sub );
  jetHadron::FixTheDamnBins( dphi_err );
  jetHadron::FixTheDamnBins( dphi_err_sub );
  
  
  // ******************************************
  // we have our errors from the TH1s now
  // write them out
  // ******************************************
  
  
  // output file
  TString outPath = path + "/sys.root";
  TFile* out = new TFile( outPath,"RECREATE");
  
  for ( int i = 0; i < deta_err.size(); ++i ) {
    deta_err[i]->Write();
    deta_err_sub[i]->Write();
    dphi_err[i]->Write();
    dphi_err_sub[i]->Write();
    
  }
  
  for ( int i = 0; i < deta_err.size(); ++i ) {
    
    deta_tow_err[i]->Write();
    deta_tow_err_sub[i]->Write();
    deta_trk_err[i]->Write();
    deta_trk_err_sub[i]->Write();
  }
  
  // now
  out->Close();
  
  return 0;
}
