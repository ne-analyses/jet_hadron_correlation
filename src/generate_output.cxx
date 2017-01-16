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



// list all input files as arguments -
// 0 = aj split bin
// 1 = output directory
// 2 = include the lowest pt bin in graphs
// 3 = corr1
// 4 = mix1
// 5 = analysis1 identifying string
// 6 = corr2
// 7 = mix2
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
  defaultCorrNames.resize(2);
  defaultCorrNames[0] = "Dijet";
  defaultCorrNames[1] = "ppDijet";
  
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
  // include low pt bin or no?
  bool includeLowPt = false;
  
  switch ( argc ) {
    case 1: { // Default case
      __OUT( "Using Default Settings" )
      
      // default files
      TFile* tmp = new TFile( "out/added/auau/trg5.6/corr.root", "READ" );
      TFile* tmpMix = new TFile( "out/added/auau/trg5.6/mix.root", "READ" );
      corrFiles.push_back( tmp );
      mixFiles.push_back( tmpMix );
      
      tmp = new TFile( "out/added/pp/trg5.6/corr.root", "READ" );
      tmpMix = new TFile( "out/added/pp/trg5.6/mix.root", "READ");
      
      corrFiles.push_back( tmp );
      mixFiles.push_back( tmpMix );

      ajSplitBin = 5;
      analysisNames = defaultCorrNames;
      outputDirBase = "/results/jet_20_10_trig_5.6";
      
      break;
    }
    default: {
      if ( (argc-5)%3 != 0 ) {
        __ERR("Need correlation file, mixing file, and analysis name for each entry")
        return -1;
      }
      std::vector<std::string> arguments( argv+1, argv+argc );
      
      // number of correlations
      int nCorrFiles = ( argc - 4 )/3;
      
      analysisNames.resize( nCorrFiles );
      
      ajSplitBin = atoi ( arguments[ 0 ].c_str() );
      outputDirBase = arguments[ 1 ];
      R = atof( arguments[ 2 ].c_str() );
      includeLowPt = atoi(arguments[3].c_str() );
      
      for ( int i = 0; i < nCorrFiles; ++i ) {
        
        TFile* tmp = new TFile( arguments[ 3*i+4 ].c_str(), "READ" );
        TFile* tmpMix =  new TFile( arguments[ (3*i)+5 ].c_str(), "READ" );
        
        corrFiles.push_back( tmp );
        mixFiles.push_back( tmpMix );
        analysisNames[i] = arguments[ (3*i)+6 ];
      }
    }
  }

  int nFiles = analysisNames.size();
  
  // put full path for output directory
  outputDirBase = currentDirectory + "/" + outputDirBase;
  
  // Build our bin selector with default settings
  jetHadron::binSelector selector;
  selector.ChangeRadius( R );
  
  // and choose whether to plot the full range or remove the lowest bin
  int graphPtBinLow = 1;
  int graphPtBinHigh = 5;
  if ( includeLowPt ) {
    graphPtBinLow = 0;
  }
  
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
  
  // building a pt bin error
  std::vector<std::vector<double> > zeros;
  zeros.resize( ptBinCenters.size() );
  for ( int i = 0; i < ptBinCenters.size(); ++i ) {
    zeros[i].resize( ptBinCenters[i].size() );
  }
  
  // Now build the correlation histograms, both the total and the aj split
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > leadingCorrelation;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > subleadingCorrelation;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > correlationAjBalanced;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > correlationAjUnbalanced;
  
  jetHadron::BuildSingleCorrelation( leadingCorrelationIn, leadingCorrelation, selector );
  jetHadron::BuildSingleCorrelation( subleadingCorrelationIn, subleadingCorrelation, selector, "sublead_uncorrsplit" );
  jetHadron::BuildAjSplitCorrelation( leadingCorrelationIn, correlationAjUnbalanced, correlationAjBalanced, selector, ajSplitBin );
  
  // get averaged correlations
  std::vector<std::vector<TH2F*> > averagedSignal = jetHadron::AverageCorrelations( leadingCorrelation, selector );
  std::vector<std::vector<TH2F*> > averagedSignalSub = jetHadron::AverageCorrelations( subleadingCorrelation, selector, "uncorr_sub" );
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
    jetHadron::Print2DHistogramsMixing( leadingMix[i], outputDirBase+"/mixing_"+analysisNames[i], analysisNames[i], selector );
    jetHadron::Print2DHistograms( averagedSignal[i], outputDirBase+"/uncorr_lead_"+analysisNames[i], analysisNames[i], selector );
    jetHadron::Print2DHistogramsEtaRestricted( averagedMixedEventCorrected[i], outputDirBase+"/avg_mix_corrected_lead_"+analysisNames[i], analysisNames[i], selector );
    jetHadron::Print2DHistogramsEtaRestricted( notAveragedMixedEventCorrected[i], outputDirBase+"/mix_corrected_lead_"+analysisNames[i], analysisNames[i], selector );
  }
  
  
  // ***************************
  // print out the 1d dEta for
  // mixing histograms
  // ***************************

  // first do the 1D dEta for mixing
  std::vector<std::vector<TH1F*> > mixingProjection = jetHadron::ProjectDeta( leadingMix, selector, "mixing_deta" );
  
  // scale by #bins in phi
  for ( int i = 0; i < mixingProjection.size(); ++i ) {
    for ( int j = 0; j < mixingProjection[i].size(); ++j ) {
      
      mixingProjection[i][j]->Scale( 1.0 / mixingProjection[i][j]->GetYaxis()->GetNbins() );
      
    }
  }
  
  // and print
  for ( int i = 0; i < mixingProjection.size(); ++i ) {
    jetHadron::Print1DHistogramsDeta( mixingProjection[i], outputDirBase+"/mixing_deta_"+analysisNames[i], analysisNames[i], selector );
  }
  
  // *************************************
  // get the 1D projection for uncorrected
  // (no mixing) for dPhi
  // *************************************
  
  std::vector<std::vector<TH1F*> > uncorrected_dphi_lead = jetHadron::ProjectDphi( averagedSignal, selector, "not_mixing_corrected", false );
  std::vector<std::vector<TH1F*> > uncorrected_dphi_sub = jetHadron::ProjectDphi( averagedSignalSub, selector, "not_mixing_corrected_sub", false );
  
  // do background subtraction
  jetHadron::SubtractBackgroundDphi( uncorrected_dphi_lead, selector );
  jetHadron::SubtractBackgroundDphi( uncorrected_dphi_sub, selector );
  
  // normalize with 1/dijets 1/bin width
  jetHadron::Normalize1D( uncorrected_dphi_lead, nEvents );
  jetHadron::Normalize1D( uncorrected_dphi_sub, nEvents );
  
  // do final fitting
  std::vector<std::vector<TF1*> > uncorrected_dphi_lead_fit = jetHadron::FitDphi( uncorrected_dphi_lead, selector );
  std::vector<std::vector<TF1*> > uncorrected_dphi_sub_fit = jetHadron::FitDphi( uncorrected_dphi_sub, selector );
  
  std::vector<std::vector<double> > uncorrected_dphi_lead_fit_yield, uncorrected_dphi_lead_fit_width, uncorrected_dphi_lead_fit_width_err, uncorrected_dphi_lead_fit_yield_err;
  std::vector<std::vector<double> > uncorrected_dphi_sub_fit_yield, uncorrected_dphi_sub_fit_width, uncorrected_dphi_sub_fit_width_err, uncorrected_dphi_sub_fit_yield_err;
  
  jetHadron::ExtractFitVals( uncorrected_dphi_lead_fit, uncorrected_dphi_lead_fit_yield, uncorrected_dphi_lead_fit_width, uncorrected_dphi_lead_fit_yield_err, uncorrected_dphi_lead_fit_width_err, selector  );
  jetHadron::ExtractFitVals( uncorrected_dphi_sub_fit, uncorrected_dphi_sub_fit_yield, uncorrected_dphi_sub_fit_width, uncorrected_dphi_sub_fit_yield_err, uncorrected_dphi_sub_fit_width_err, selector  );
  
  // now overlay and save
  jetHadron::Print1DHistogramsOverlayedDphiWFit( uncorrected_dphi_lead, uncorrected_dphi_lead_fit, outputDirBase+"/uncorrected_dphi_lead"+analysisNames[0], analysisNames, selector );
  jetHadron::Print1DHistogramsOverlayedDphiWFit( uncorrected_dphi_sub, uncorrected_dphi_sub_fit, outputDirBase+"/uncorrected_dphi_sub"+analysisNames[0], analysisNames, selector );
  jetHadron::PrintGraphWithErrors( ptBinCenters, uncorrected_dphi_lead_fit_yield, zeros, uncorrected_dphi_lead_fit_yield_err, outputDirBase+"/uncorrected_dphi_lead_graph", analysisNames, "Trigger Jet Yields", selector,  graphPtBinLow, graphPtBinHigh );
  jetHadron::PrintGraphWithErrors( ptBinCenters, uncorrected_dphi_sub_fit_yield, zeros, uncorrected_dphi_sub_fit_yield_err, outputDirBase+"/uncorrected_dphi_sub_graph", analysisNames, "Recoil Jet Yields", selector, graphPtBinLow, graphPtBinHigh );
  // *************************************
  // Mixing corrected stuff -
  // first Subtracted DPhi
  // *************************************
  // define what "regions" we want the subtraction to be done in
  double subtractionRegions[4] = { -1.0, -0.6, 0.6, 1.0 };

  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted = jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrected, selector, subtractionRegions, "mixing_corrected_near_far_sub_dphi", true );
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted_sub = jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrectedSub, selector, subtractionRegions, "mixing_corrected_near_far_sub_dphi_sub", true  );
  
  // and to get the individual near and far histograms
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted_near, corrected_dphi_subtracted_far;
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted_sub_near, corrected_dphi_subtracted_sub_far;
  jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrected, corrected_dphi_subtracted_near, corrected_dphi_subtracted_far, selector, subtractionRegions, "mixing_corrected_near_far_sub_dphi", true );
  jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrectedSub, corrected_dphi_subtracted_sub_near, corrected_dphi_subtracted_sub_far, selector, subtractionRegions, "mixing_corrected_near_far_sub_dphi_sub", true  );
  
  
  // do background subtraction
  jetHadron::SubtractBackgroundDphi( corrected_dphi_subtracted, selector );
  jetHadron::SubtractBackgroundDphi( corrected_dphi_subtracted_sub, selector );
  
  // normalize with 1/dijets 1/bin width
  jetHadron::Normalize1D( corrected_dphi_subtracted, nEvents );
  jetHadron::Normalize1D( corrected_dphi_subtracted_sub, nEvents );
  
  // normalize
  jetHadron::Normalize1D( corrected_dphi_subtracted_near, nEvents );
  jetHadron::Normalize1D( corrected_dphi_subtracted_far, nEvents );
  jetHadron::Normalize1D( corrected_dphi_subtracted_sub_near, nEvents );
  jetHadron::Normalize1D( corrected_dphi_subtracted_sub_far, nEvents );
  
  // do final fitting
  std::vector<std::vector<TF1*> > corrected_dphi_subtracted_fit = jetHadron::FitDphiRestricted( corrected_dphi_subtracted, selector );
  std::vector<std::vector<TF1*> > corrected_dphi_subtracted_sub_fit = jetHadron::FitDphiRestricted( corrected_dphi_subtracted_sub, selector );
  
  std::vector<std::vector<double> > corrected_dphi_subtracted_fit_yield, corrected_dphi_subtracted_fit_width, corrected_dphi_subtracted_fit_width_err, corrected_dphi_subtracted_fit_yield_err;
  std::vector<std::vector<double> > corrected_dphi_subtracted_sub_fit_yield, corrected_dphi_subtracted_sub_fit_width, corrected_dphi_subtracted_sub_fit_width_err, corrected_dphi_subtracted_sub_fit_yield_err;
  
  jetHadron::ExtractFitVals( corrected_dphi_subtracted_fit, corrected_dphi_subtracted_fit_yield, corrected_dphi_subtracted_fit_width, corrected_dphi_subtracted_fit_yield_err, corrected_dphi_subtracted_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_dphi_subtracted_sub_fit, corrected_dphi_subtracted_sub_fit_yield, corrected_dphi_subtracted_sub_fit_width, corrected_dphi_subtracted_sub_fit_yield_err, corrected_dphi_subtracted_sub_fit_width_err, selector  );
  
  
  // now overlay and save
  jetHadron::Print1DHistogramsOverlayedDphiWFitRestricted( corrected_dphi_subtracted, corrected_dphi_subtracted_fit, outputDirBase+"/corrected_dphi_subtracted_lead"+analysisNames[0], analysisNames, selector );
  jetHadron::Print1DHistogramsOverlayedDphiWFitRestricted( corrected_dphi_subtracted_sub, corrected_dphi_subtracted_sub_fit, outputDirBase+"/corrected_dphi_subtracted_sub"+analysisNames[0], analysisNames, selector );
  jetHadron::PrintGraphWithErrors( ptBinCenters, corrected_dphi_subtracted_fit_yield, zeros, corrected_dphi_subtracted_fit_yield_err, outputDirBase+"/corrected_dphi_subtracted_graph", analysisNames, "Trigger Jet Yields", selector, graphPtBinLow, graphPtBinHigh );
  jetHadron::PrintGraphWithErrors( ptBinCenters, corrected_dphi_subtracted_sub_fit_yield, zeros, corrected_dphi_subtracted_sub_fit_yield_err, outputDirBase+"/corrected_dphi_subtracted_sub_graph", analysisNames, "Recoil Jet Yields", selector, graphPtBinLow, graphPtBinHigh );
  
  // overlay and save near/far
  for ( int i = 0; i < nFiles; ++i ) {
    std::vector<std::string> tmpVec;
    tmpVec.push_back("near");
    tmpVec.push_back("far");
    jetHadron::Print1DHistogramsOverlayedDphiOther( corrected_dphi_subtracted_far[i], corrected_dphi_subtracted_near[i], outputDirBase+"/near_overlay_"+analysisNames[i], tmpVec[0], tmpVec[1], selector );
    jetHadron::Print1DHistogramsOverlayedDphiOther( corrected_dphi_subtracted_sub_far[i], corrected_dphi_subtracted_sub_near[i], outputDirBase+"/far_overlay_"+analysisNames[i], tmpVec[0], tmpVec[1], selector );
  }

  
  // Now we will do not subtracted projections
  // and dEta
  std::vector<std::vector<TH1F*> > corrected_dphi_lead = jetHadron::ProjectDphi( averagedMixedEventCorrected, selector, "mixing_corrected_dphi", true );
  std::vector<std::vector<TH1F*> > corrected_dphi_sub = jetHadron::ProjectDphi( averagedMixedEventCorrectedSub, selector, "mixing_corrected_dphi_sub", true );
  std::vector<std::vector<TH1F*> > corrected_deta_lead = jetHadron::ProjectDeta( averagedMixedEventCorrected, selector, "mixing_corrected_deta", true );
  std::vector<std::vector<TH1F*> > corrected_deta_sub = jetHadron::ProjectDeta( averagedMixedEventCorrectedSub, selector, "mixing_corrected_deta_sub", true );
  
  // do background subtraction
  jetHadron::SubtractBackgroundDphi( corrected_dphi_lead, selector );
  jetHadron::SubtractBackgroundDphi( corrected_dphi_sub, selector );
  jetHadron::SubtractBackgroundDeta( corrected_deta_lead, selector );
  jetHadron::SubtractBackgroundDeta( corrected_deta_sub, selector );
  
  // normalize with 1/dijets 1/bin width
  jetHadron::Normalize1D( corrected_dphi_lead, nEvents );
  jetHadron::Normalize1D( corrected_dphi_sub, nEvents );
  jetHadron::Normalize1D( corrected_deta_lead, nEvents );
  jetHadron::Normalize1D( corrected_deta_sub, nEvents );
  
  // do final fitting
  std::vector<std::vector<TF1*> > corrected_dphi_lead_fit = jetHadron::FitDphi( corrected_dphi_lead, selector );
  std::vector<std::vector<TF1*> > corrected_dphi_sub_fit = jetHadron::FitDphi( corrected_dphi_sub, selector );
  std::vector<std::vector<TF1*> > corrected_deta_lead_fit = jetHadron::FitDeta( corrected_deta_lead, selector );
  std::vector<std::vector<TF1*> > corrected_deta_sub_fit = jetHadron::FitDeta( corrected_deta_sub, selector );
  
  std::vector<std::vector<double> > corrected_dphi_fit_yield, corrected_dphi_fit_width, corrected_dphi_fit_width_err, corrected_dphi_fit_yield_err;
  std::vector<std::vector<double> > corrected_dphi_sub_fit_yield, corrected_dphi_sub_fit_width, corrected_dphi_sub_fit_width_err, corrected_dphi_sub_fit_yield_err;
  std::vector<std::vector<double> > corrected_deta_fit_yield, corrected_deta_fit_width, corrected_deta_fit_width_err, corrected_deta_fit_yield_err;
  std::vector<std::vector<double> > corrected_deta_sub_fit_yield, corrected_deta_sub_fit_width, corrected_deta_sub_fit_width_err, corrected_deta_sub_fit_yield_err;
  
  jetHadron::ExtractFitVals( corrected_dphi_lead_fit, corrected_dphi_fit_yield, corrected_dphi_fit_width, corrected_dphi_fit_yield_err, corrected_dphi_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_dphi_sub_fit, corrected_dphi_sub_fit_yield, corrected_dphi_sub_fit_width, corrected_dphi_sub_fit_yield_err, corrected_dphi_sub_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_deta_lead_fit, corrected_deta_fit_yield, corrected_dphi_fit_width, corrected_deta_fit_yield_err, corrected_deta_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_deta_sub_fit, corrected_deta_sub_fit_yield, corrected_deta_sub_fit_width, corrected_deta_sub_fit_yield_err, corrected_deta_sub_fit_width_err, selector  );
  
  // now overlay and save
  jetHadron::Print1DHistogramsOverlayedDphiWFit( corrected_dphi_lead, corrected_dphi_lead_fit, outputDirBase+"/corrected_dphi_lead"+analysisNames[0], analysisNames, selector );
  jetHadron::Print1DHistogramsOverlayedDphiWFit( corrected_dphi_sub, corrected_dphi_sub_fit, outputDirBase+"/corrected_dphi_sub"+analysisNames[0], analysisNames, selector );
  jetHadron::Print1DHistogramsOverlayedDetaWFitRestricted( corrected_deta_lead, corrected_deta_lead_fit, outputDirBase+"/corrected_deta_lead"+analysisNames[0], analysisNames, selector );
  jetHadron::Print1DHistogramsOverlayedDetaWFitRestricted( corrected_deta_sub, corrected_deta_sub_fit, outputDirBase+"/corrected_deta_sub"+analysisNames[0], analysisNames, selector );
  
  jetHadron::PrintGraphWithErrors( ptBinCenters, corrected_dphi_fit_yield, zeros, corrected_dphi_fit_yield_err, outputDirBase+"/corrected_dphi_graph", analysisNames, "Trigger Jet Yields", selector, graphPtBinLow, graphPtBinHigh );
  jetHadron::PrintGraphWithErrors( ptBinCenters, corrected_dphi_sub_fit_yield, zeros, corrected_dphi_sub_fit_yield_err, outputDirBase+"/corrected_dphi_sub_graph", analysisNames, "Recoil Jet Yields", selector, graphPtBinLow, graphPtBinHigh );
  jetHadron::PrintGraphWithErrors( ptBinCenters, corrected_deta_fit_yield, zeros, corrected_deta_fit_yield_err, outputDirBase+"/corrected_deta_graph", analysisNames, "Trigger Jet Yields", selector, graphPtBinLow, graphPtBinHigh );
  jetHadron::PrintGraphWithErrors( ptBinCenters, corrected_deta_sub_fit_yield, zeros, corrected_deta_sub_fit_yield_err, outputDirBase+"/corrected_deta_sub_graph", analysisNames, "Recoil Jet Yields", selector, graphPtBinLow, graphPtBinHigh );

  // *******************************
  // get 1D projections for Aj split
  // *******************************
  std::vector<std::vector<TH1F*> > aj_balanced_dphi = jetHadron::ProjectDphi( averagedSignalBalanced, selector, "aj_balanced_", false );
  std::vector<std::vector<TH1F*> > aj_unbalanced_dphi = jetHadron::ProjectDphi( averagedSignalUnbalanced, selector, "aj_unbalanced_", false );
  
  jetHadron::Normalize1DAjSplit( aj_balanced_dphi, nEvents, 1, ajSplitBin );
  jetHadron::Normalize1DAjSplit( aj_unbalanced_dphi, nEvents, ajSplitBin+1, 20 );
  
  // now we need to do background subtraction
  jetHadron::SubtractBackgroundDphi( aj_balanced_dphi, selector );
  jetHadron::SubtractBackgroundDphi( aj_unbalanced_dphi, selector );
  // now do the subtraction
  std::vector<std::vector<TH1F*> > aj_subtracted = jetHadron::Subtract1D( aj_balanced_dphi, aj_unbalanced_dphi, "aj_split" );
  
  // and print it out
  jetHadron::Print1DHistogramsOverlayedDphi( aj_subtracted, outputDirBase+"/aj_subtracted_dif", analysisNames, selector );
  // also print out the overlayed
  for ( int i = 0; i < nFiles; ++i ) {
    std::vector<std::string> tmpVec;
    tmpVec.push_back("balanced");
    tmpVec.push_back("unbalanced");
    jetHadron::Print1DHistogramsOverlayedDphiOther( aj_balanced_dphi[i], aj_unbalanced_dphi[i], outputDirBase+"/aj_dif_dphi"+analysisNames[i], tmpVec[0], tmpVec[1], selector );
  }
  
  
  //******************************************
  // Now we need to test the systematic errors
  //******************************************
  TFile sysIn( "out/added/pp/trg5.6/sys.root", "READ");
  std::vector<TH1F*> deta_sys, deta_sys_sub, dphi_sys, dphi_sys_sub;

  for ( int i = 0; i < corrected_dphi_lead[1].size(); ++i ) {
    std::string tmpdEta = "deta_sys_quad_pt_" + patch::to_string(i);
    std::string tmpdEtaSub = "sub_deta_sys_quad_pt_" + patch::to_string(i);
    std::string tmpdPhi = "dphi_sys_quad_pt_" + patch::to_string(i);
    std::string tmpdPhiSub = "sub_dphi_sys_quad_pt_" + patch::to_string(i);
    
    deta_sys.push_back( (TH1F*) sysIn.Get( tmpdEta.c_str() ) );
    deta_sys_sub.push_back( (TH1F*) sysIn.Get( tmpdEtaSub.c_str() ) );
    dphi_sys.push_back( (TH1F*) sysIn.Get( tmpdPhi.c_str() ) );
    dphi_sys_sub.push_back( (TH1F*) sysIn.Get( tmpdPhiSub.c_str() ) );
    
  }
  
  // and do printouts
  jetHadron::Print1DDPhiHistogramsWithSysErr( corrected_dphi_lead[1], dphi_sys, selector, outputDirBase+"/dphi_sys_lead", -0.8, 0.8  );
  jetHadron::Print1DDPhiHistogramsWithSysErr( corrected_dphi_sub[1], dphi_sys_sub, selector, outputDirBase+"/dphi_sys_sub", -0.8 , 0.8  );
  jetHadron::Print1DDEtaHistogramsWithSysErr( corrected_deta_lead[1], deta_sys, selector, outputDirBase+"/deta_sys_lead", -0.8 , 0.8  );
  jetHadron::Print1DDEtaHistogramsWithSysErr( corrected_deta_sub[1], deta_sys_sub, selector, outputDirBase+"/deta_sys_sub", -0.8 , 0.8  );
  
  // generate some 5% error histograms for each of the histograms we use
  std::vector<std::vector<TH1F*> > dphi_yield_err = jetHadron::BuildYieldError( corrected_dphi_lead, selector, analysisNames, "dphi_lead_yield_err" );
  std::vector<std::vector<TH1F*> > dphi_sub_yield_err = jetHadron::BuildYieldError( corrected_dphi_sub, selector, analysisNames, "dphi_sub_yield_err" );
  std::vector<std::vector<TH1F*> > deta_yield_err = jetHadron::BuildYieldError( corrected_deta_lead, selector, analysisNames, "deta_lead_yield_err" );
  std::vector<std::vector<TH1F*> > deta_sub_yield_err = jetHadron::BuildYieldError( corrected_deta_sub, selector, analysisNames, "deta_sub_yield_err" );
  
  // ******************************************
  // now we're testing yields from bin counting
  // ******************************************
  std::vector<std::vector<double> > dphi_lead_bin_int, dphi_sub_bin_int, deta_lead_bin_int, deta_sub_bin_int;
  std::vector<std::vector<double> > dphi_lead_bin_int_err, dphi_sub_bin_int_err, deta_lead_bin_int_err, deta_sub_bin_int_err;
  
  jetHadron::ExtractIntegral( corrected_dphi_lead, dphi_lead_bin_int, dphi_lead_bin_int_err, selector, -0.6, 0.6 );
  jetHadron::ExtractIntegral( corrected_dphi_sub, dphi_sub_bin_int, dphi_sub_bin_int_err, selector, -0.6, 0.6 );
  jetHadron::ExtractIntegral( corrected_deta_lead, deta_lead_bin_int, deta_lead_bin_int_err, selector, -0.6, 0.6 );
  jetHadron::ExtractIntegral( corrected_deta_sub, deta_sub_bin_int, deta_sub_bin_int_err, selector, -0.6, 0.6 );
  
  std::vector<TGraphErrors*> dphi_yield_graph = jetHadron::MakeGraphs( ptBinCenters, dphi_lead_bin_int, zeros, dphi_lead_bin_int_err, 1, 5, selector, analysisNames, "dphi" );
  std::vector<TGraphErrors*> dphi_sub_yield_graph = jetHadron::MakeGraphs( ptBinCenters, dphi_sub_bin_int, zeros, dphi_sub_bin_int_err, 1, 5, selector, analysisNames, "dphi_sub" );
  std::vector<TGraphErrors*> deta_yield_graph = jetHadron::MakeGraphs( ptBinCenters, deta_lead_bin_int, zeros, deta_lead_bin_int_err, 1, 5, selector, analysisNames, "deta" );
  std::vector<TGraphErrors*> deta_sub_yield_graph = jetHadron::MakeGraphs( ptBinCenters, deta_sub_bin_int, zeros, deta_sub_bin_int_err, 1, 5, selector, analysisNames, "deta_sub" );
  
  // and to do that we need some errors as well! We already have relative systematics, but need the yields for the 5% tracking
  std::vector<std::vector<double> > dphi_lead_sys_bin_int, dphi_sub_sys_bin_int, deta_lead_sys_bin_int, deta_sub_sys_bin_int;
  std::vector<std::vector<double> > dphi_lead_sys_bin_int_err, dphi_sub_sys_bin_int_err, deta_lead_sys_bin_int_err, deta_sub_sys_bin_int_err;
  
  jetHadron::ExtractIntegral( dphi_yield_err, dphi_lead_sys_bin_int, dphi_lead_sys_bin_int_err, selector, -0.6, 0.6 );
  jetHadron::ExtractIntegral( dphi_sub_yield_err, dphi_sub_sys_bin_int, dphi_sub_sys_bin_int_err, selector, -0.6, 0.6 );
  jetHadron::ExtractIntegral( deta_yield_err, deta_lead_sys_bin_int, deta_lead_sys_bin_int_err, selector, -0.6, 0.6 );
  jetHadron::ExtractIntegral( deta_sub_yield_err, deta_sub_sys_bin_int, deta_sub_sys_bin_int_err, selector, -0.6, 0.6 );
  
  std::vector<TGraphErrors*> dphi_yield_sys_graph = jetHadron::MakeGraphs( ptBinCenters, dphi_lead_sys_bin_int, zeros, dphi_lead_sys_bin_int_err, 1, 5, selector, analysisNames, "dphi_sys" );
  std::vector<TGraphErrors*> dphi_sub_yield_sys_graph = jetHadron::MakeGraphs( ptBinCenters, dphi_sub_sys_bin_int, zeros, dphi_sub_sys_bin_int_err, 1, 5, selector, analysisNames, "dphi_sub_sys" );
  std::vector<TGraphErrors*> deta_yield_sys_graph = jetHadron::MakeGraphs( ptBinCenters, deta_lead_sys_bin_int, zeros, deta_lead_sys_bin_int_err, 1, 5, selector, analysisNames, "deta_sys" );
  std::vector<TGraphErrors*> deta_sub_yield_sys_graph = jetHadron::MakeGraphs( ptBinCenters, deta_sub_sys_bin_int, zeros, deta_sub_sys_bin_int_err, 1, 5, selector, analysisNames, "deta_sub_sys" );
  
  
  return 0;
}


