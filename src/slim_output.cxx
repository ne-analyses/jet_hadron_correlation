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
#include "TError.h"

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
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetOptTitle(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing( 1.0 );
  gStyle->SetHatchesLineWidth( 2 );
  gErrorIgnoreLevel = kError;

  
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
  std::vector<TFile*> corrFilesHard;
  std::vector<TFile*> mixFilesHard;
  std::vector<std::string> analysisNames;
  std::vector<std::string> analysisNamesHard;
  
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
      TFile* tmp = new TFile( "out/added/auau/trg6.00/corrv5.root", "READ" );
      TFile* tmpMix = new TFile( "out/added/auau/trg6.00/mix.root", "READ" );
      corrFiles.push_back( tmp );
      mixFiles.push_back( tmpMix );
      
      tmp = new TFile( "out/added/pp/trg6.00/corrv5.root", "READ" );
      tmpMix = new TFile( "out/added/pp/trg6.00/mixv4.root", "READ");
      
      corrFiles.push_back( tmp );
      mixFiles.push_back( tmpMix );
      
      tmp = new TFile( "out/added/ppembedhard/trg6.00/corrv5.root");
      tmpMix = new TFile( "out/added/ppembedhard/trg6.00/mixv4.root");
      
      corrFilesHard.push_back( tmp );
      mixFilesHard.push_back( tmpMix );
      analysisNamesHard.push_back( "ppDijetHard" );

      ajSplitBin = 5;
      analysisNames = defaultCorrNames;
      outputDirBase = "/results/slim_jet_20_10_trig_6";
      
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
  for ( int i = 0; i < 6; ++i ) {
    std::cout<<"pt bin: " << i << " low pt: " << selector.ptBinEdgeLo[i] << " high pt: " << selector.ptBinEdgeHi[i] << std::endl;
  }
  
  // and choose whether to plot the full range or remove the lowest bin
  int graphPtBinLow = 1;
  int graphPtBinHigh = 5;
  if ( includeLowPt ) {
    graphPtBinLow = 0;
  }
  
  // Build our initial histogram holders
  std::vector<TH3F*> nEvents;
  std::vector<TH3F*> nEventsHard;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leadingCorrelationIn;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > subleadingCorrelationIn;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leadingCorrelationInHard;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > subleadingCorrelationInHard;
  // and event mixing
  std::vector<TH3F*> nEventsMixing;
  std::vector<TH3F*> nEventsMixingHard;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leadingMixIn;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > subleadingMixIn;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leadingMixInHard;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > subleadingMixInHard;
  // reading in the histograms
  jetHadron::ReadInFiles( corrFiles, leadingCorrelationIn, subleadingCorrelationIn, nEvents, selector );
  jetHadron::ReadInFiles( corrFilesHard, leadingCorrelationInHard, subleadingCorrelationInHard, nEventsHard, selector, "pp_hard" );
  jetHadron::ReadInFilesMix( mixFiles, leadingMixIn, subleadingMixIn, nEventsMixing, selector );
  jetHadron::ReadInFilesMix( mixFilesHard, leadingMixInHard, subleadingMixInHard, nEventsMixingHard, selector, "pp_hard" );
  
  // Find the pt bin center for future use
  std::vector<TH1F*> ptSpectra;
  std::vector<std::vector<double> > ptBinCenters = jetHadron::FindPtBinCenter( leadingCorrelationIn, ptSpectra, selector );
  
  // building veros vector for graphs later
  std::vector<std::vector<double> > zeros;
  zeros.resize( ptBinCenters.size() );
  for ( int i = 0; i < ptBinCenters.size(); ++i ) {
    zeros[i].resize( ptBinCenters[i].size() );
  }
  
  // Now build the correlation histograms, both the total and the aj split
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > leadingCorrelation;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > subleadingCorrelation;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > leadingCorrelationHard;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > subleadingCorrelationHard;
  
  jetHadron::BuildSingleCorrelation( leadingCorrelationIn, leadingCorrelation, selector, "lead_uncorr" );
  jetHadron::BuildSingleCorrelation( subleadingCorrelationIn, subleadingCorrelation, selector, "sublead_uncorr" );
  jetHadron::BuildSingleCorrelation( leadingCorrelationInHard, leadingCorrelationHard, selector, "lead_uncorr_pp_hard" );
  jetHadron::BuildSingleCorrelation( subleadingCorrelationInHard, subleadingCorrelationHard, selector, "sublead_uncorr_pp_hard" );
  
  // get averaged correlations
  std::vector<std::vector<TH2F*> > averagedSignal = jetHadron::AverageCorrelations( leadingCorrelation, selector, "uncorr_avg" );
  std::vector<std::vector<TH2F*> > averagedSignalSub = jetHadron::AverageCorrelations( subleadingCorrelation, selector, "uncorr_sub_avg" );
  
  // Now build and scale the event mixing histograms.
  // we will use the averaged event mixing for now
  // First average over all centrality/vz/aj, project into pt
  std::vector<std::vector<TH2F*> > leadingMix =  jetHadron::RecombineMixedEvents( leadingMixIn, selector, "avg_mix_" );
  std::vector<std::vector<TH2F*> > subleadingMix = jetHadron::RecombineMixedEvents( subleadingMixIn, selector, "avg_mix_sub" );
  std::vector<std::vector<TH2F*> > leadingMixHard =  jetHadron::RecombineMixedEvents( leadingMixInHard, selector, "avg_mix_hard" );
  std::vector<std::vector<TH2F*> > subleadingMixHard = jetHadron::RecombineMixedEvents( subleadingMixInHard, selector, "avg_mix_sub_hard" );
  
  // And Scale so that MaxBinContent = 1
  jetHadron::ScaleMixedEvents( leadingMix );
  jetHadron::ScaleMixedEvents( subleadingMix );
  jetHadron::ScaleMixedEvents( leadingMixHard );
  jetHadron::ScaleMixedEvents( subleadingMixHard );
  
  // And get the event mixing corrected both averaged or not
  std::vector<std::vector<TH2F*> > averagedMixedEventCorrected = jetHadron::EventMixingCorrection( leadingCorrelation, leadingMix, selector, "leading_avg"  );
  std::vector<std::vector<TH2F*> > averagedMixedEventCorrectedSub = jetHadron::EventMixingCorrection( subleadingCorrelation, subleadingMix, selector, "subleading_avg"  );
  
  std::vector<std::vector<TH2F*> > averagedMixedEventCorrectedHard = jetHadron::EventMixingCorrection( leadingCorrelationHard, leadingMixHard, selector, "leading_avg_hard"  );
  std::vector<std::vector<TH2F*> > averagedMixedEventCorrectedSubHard = jetHadron::EventMixingCorrection( subleadingCorrelationHard, subleadingMixHard, selector, "subleading_avg_hard"  );
  
  
  // ***************************
  // print out the 2d histograms
  // ***************************
  for ( int i = 0; i < nFiles; ++ i ) {
    jetHadron::Print2DHistogramsMixing( leadingMix[i], outputDirBase+"/mixing_"+analysisNames[i], analysisNames[i], selector );
    jetHadron::Print2DHistograms( averagedSignal[i], outputDirBase+"/uncorr_lead_"+analysisNames[i], analysisNames[i], selector );
    jetHadron::Print2DHistogramsEtaRestricted( averagedMixedEventCorrected[i], outputDirBase+"/avg_mix_corrected_lead_"+analysisNames[i], analysisNames[i], selector );
  }
  
  for ( int i = 0; i < nEventsHard.size(); ++i ) {
    jetHadron::Print2DHistogramsMixing( leadingMixHard[i], outputDirBase+"/mixing_"+analysisNamesHard[i], analysisNamesHard[i], selector );
    jetHadron::Print2DHistogramsEtaRestricted( averagedMixedEventCorrectedHard[i], outputDirBase+"/avg_mix_corrected_lead_"+analysisNamesHard[i], analysisNamesHard[i], selector );
  }
  // now turn off the titles
  gStyle->SetOptTitle(0);
  
  
  // ***************************
  // print out the 1d dEta for
  // mixing histograms
  // ***************************

  // first do the 1D dEta for mixing
  std::vector<std::vector<TH1F*> > mixingProjection = jetHadron::ProjectDetaMixing( leadingMix, selector, "mixing_deta" );
  
  // scale by #bins in phi
  for ( int i = 0; i < mixingProjection.size(); ++i ) {
    for ( int j = 0; j < mixingProjection[i].size(); ++j ) {
      
      mixingProjection[i][j]->Scale( 1.0 / leadingMix[i][j]->GetYaxis()->GetNbins() );
      
    }
  }
  
  // and print
  for ( int i = 0; i < mixingProjection.size(); ++i ) {
    jetHadron::Print1DHistogramsDeta( mixingProjection[i], outputDirBase+"/mixing_deta_"+analysisNames[i], analysisNames[i], selector );
  }
  
  __OUT("Starting 1D Correlations")
  
  // *************************************
  // correlations start here!!
  // *************************************
  
  // *************************************
  // Mixing corrected stuff -
  // first Subtracted DPhi
  // *************************************
  // define what "regions" we want the subtraction to be done in
  
  __OUT("now we get near minus far projections for DPhi")
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted = jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrected, selector,  "mixing_corrected_near_far_sub_dphi" );
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted_sub = jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrectedSub, selector,  "mixing_corrected_near_far_sub_dphi_sub"  );
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted_extended = jetHadron::ProjectDphiNearMinusFarExtended( averagedMixedEventCorrected, selector,  "mixing_corrected_near_far_sub_dphi_extended" );
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted_sub_extended = jetHadron::ProjectDphiNearMinusFarExtended( averagedMixedEventCorrectedSub, selector,  "mixing_corrected_near_far_sub_dphi_sub_extended"  );
  
  // and to get the individual near and far histograms
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted_near, corrected_dphi_subtracted_far;
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted_sub_near, corrected_dphi_subtracted_sub_far;
  jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrected, corrected_dphi_subtracted_near, corrected_dphi_subtracted_far, selector,  "mixing_corrected_near_far_sub_dphi" );
  jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrectedSub, corrected_dphi_subtracted_sub_near, corrected_dphi_subtracted_sub_far, selector,  "mixing_corrected_near_far_sub_dphi_sub"  );
  
  __OUT("normalize the subtracted histograms")
  // normalize with 1/dijets 1/bin width
  jetHadron::Normalize1D( corrected_dphi_subtracted, nEvents );
  jetHadron::Normalize1D( corrected_dphi_subtracted_sub, nEvents );
  jetHadron::Normalize1D( corrected_dphi_subtracted_extended, nEvents );
  jetHadron::Normalize1D( corrected_dphi_subtracted_sub_extended, nEvents );
  
  jetHadron::Normalize1D( corrected_dphi_subtracted_near, nEvents );
  jetHadron::Normalize1D( corrected_dphi_subtracted_far, nEvents );
  jetHadron::Normalize1D( corrected_dphi_subtracted_sub_near, nEvents );
  jetHadron::Normalize1D( corrected_dphi_subtracted_sub_far, nEvents );
  
  // make sure the bins are correct...
  jetHadron::FixTheDamnBins( corrected_dphi_subtracted );
  jetHadron::FixTheDamnBins( corrected_dphi_subtracted_sub );
  jetHadron::FixTheDamnBins( corrected_dphi_subtracted_extended );
  jetHadron::FixTheDamnBins( corrected_dphi_subtracted_sub_extended );
  jetHadron::FixTheDamnBins( corrected_dphi_subtracted_near );
  jetHadron::FixTheDamnBins( corrected_dphi_subtracted_far );
  jetHadron::FixTheDamnBins( corrected_dphi_subtracted_sub_near );
  jetHadron::FixTheDamnBins( corrected_dphi_subtracted_sub_far );
  
  
  __OUT("do final fitting and extract fit values")
  // do final fitting
  std::vector<std::vector<TF1*> > corrected_dphi_subtracted_fit = jetHadron::FitDphiRestricted( corrected_dphi_subtracted, selector );
  std::vector<std::vector<TF1*> > corrected_dphi_subtracted_sub_fit = jetHadron::FitDphiRestricted( corrected_dphi_subtracted_sub, selector );
  
  std::vector<std::vector<double> > corrected_dphi_subtracted_fit_yield, corrected_dphi_subtracted_fit_width, corrected_dphi_subtracted_fit_width_err, corrected_dphi_subtracted_fit_yield_err;
  std::vector<std::vector<double> > corrected_dphi_subtracted_sub_fit_yield, corrected_dphi_subtracted_sub_fit_width, corrected_dphi_subtracted_sub_fit_width_err, corrected_dphi_subtracted_sub_fit_yield_err;
  
  jetHadron::ExtractFitVals( corrected_dphi_subtracted_fit, corrected_dphi_subtracted_fit_yield, corrected_dphi_subtracted_fit_width, corrected_dphi_subtracted_fit_yield_err, corrected_dphi_subtracted_fit_width_err, selector  );
  jetHadron::ExtractFitVals( corrected_dphi_subtracted_sub_fit, corrected_dphi_subtracted_sub_fit_yield, corrected_dphi_subtracted_sub_fit_width, corrected_dphi_subtracted_sub_fit_yield_err, corrected_dphi_subtracted_sub_fit_width_err, selector  );
  
  __OUT("print out dphi near minus far subtracted")
  // now overlay and save
  jetHadron::Print1DHistogramsOverlayedDphiWFitRestricted( corrected_dphi_subtracted, corrected_dphi_subtracted_fit, outputDirBase+"/corrected_dphi_subtracted_lead"+analysisNames[0], analysisNames, selector );
  jetHadron::Print1DHistogramsOverlayedDphiWFitRestricted( corrected_dphi_subtracted_sub, corrected_dphi_subtracted_sub_fit, outputDirBase+"/corrected_dphi_subtracted_sub"+analysisNames[0], analysisNames, selector );
  jetHadron::Print1DHistogramsOverlayedDphiWFitRestricted( corrected_dphi_subtracted_extended, corrected_dphi_subtracted_fit, outputDirBase+"/corrected_dphi_subtracted_lead_extended"+analysisNames[0], analysisNames, selector );
  jetHadron::Print1DHistogramsOverlayedDphiWFitRestricted( corrected_dphi_subtracted_sub_extended, corrected_dphi_subtracted_sub_fit, outputDirBase+"/corrected_dphi_subtracted_sub_extended"+analysisNames[0], analysisNames, selector );
  
  __OUT("save near and far projections overlayed")
  // overlay and save near/far
  for ( int i = 0; i < nFiles; ++i ) {
    std::vector<std::string> tmpVec;
    tmpVec.push_back("far");
    tmpVec.push_back("near");
    jetHadron::Print1DHistogramsOverlayedDphiOther( corrected_dphi_subtracted_far[i], corrected_dphi_subtracted_near[i], outputDirBase+"/near_far_overlay_lead_"+analysisNames[i], tmpVec[0], tmpVec[1], selector );
    jetHadron::Print1DHistogramsOverlayedDphiOther( corrected_dphi_subtracted_sub_far[i], corrected_dphi_subtracted_sub_near[i], outputDirBase+"/near_far_overlay_sub_"+analysisNames[i], tmpVec[0], tmpVec[1], selector );
  }
  

  __OUT("now, do projections for deta")
  // Now we will do not subtracted projections
  // and dEta
  std::vector<std::vector<TH1F*> > corrected_deta_lead = jetHadron::ProjectDeta( averagedMixedEventCorrected, selector, "mixing_corrected_deta" );
  std::vector<std::vector<TH1F*> > corrected_deta_sub = jetHadron::ProjectDeta( averagedMixedEventCorrectedSub, selector, "mixing_corrected_deta_sub" );
  std::vector<std::vector<TH1F*> > corrected_deta_lead_extended = jetHadron::ProjectDetaExtended( averagedMixedEventCorrected, selector, "mixing_corrected_deta" );
  std::vector<std::vector<TH1F*> > corrected_deta_sub_extended = jetHadron::ProjectDetaExtended( averagedMixedEventCorrectedSub, selector, "mixing_corrected_deta_sub" );
  
  __OUT("Background subtraction")
  // do background subtraction
  jetHadron::SubtractBackgroundDeta( corrected_deta_lead, selector );
  jetHadron::SubtractBackgroundDeta( corrected_deta_sub, selector );
  jetHadron::SubtractBackgroundDeta( corrected_deta_lead_extended, selector );
  jetHadron::SubtractBackgroundDeta( corrected_deta_sub_extended, selector );
  
  __OUT("Normalization")
  // normalize with 1/dijets 1/bin width
  jetHadron::Normalize1D( corrected_deta_lead, nEvents );
  jetHadron::Normalize1D( corrected_deta_sub, nEvents );
  jetHadron::Normalize1D( corrected_deta_lead_extended, nEvents );
  jetHadron::Normalize1D( corrected_deta_sub_extended, nEvents );
  
  
  // ************************************
  // ************************************
  // we now redo all the fitting above
  // to make the secondary pp with hard
  // embedding included to subtract from
  // the auau
  // ************************************
  // ************************************
  __OUT("performing similar computations on the pp reference used for AuAu subtraction")
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted_hard = jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrectedHard, selector,  "mixing_corrected_near_far_sub_dphi_hard" );
  std::vector<std::vector<TH1F*> > corrected_dphi_subtracted_sub_hard = jetHadron::ProjectDphiNearMinusFar( averagedMixedEventCorrectedSubHard, selector,  "mixing_corrected_near_far_sub_dphi_sub_hard"  );
  std::vector<std::vector<TH1F*> > corrected_deta_lead_hard = jetHadron::ProjectDeta( averagedMixedEventCorrectedHard, selector, "mixing_corrected_deta_hard" );
  std::vector<std::vector<TH1F*> > corrected_deta_sub_hard = jetHadron::ProjectDeta( averagedMixedEventCorrectedSubHard, selector, "mixing_corrected_deta_sub_hard" );
  
  // now do the subtraction for dEta
  jetHadron::SubtractBackgroundDeta( corrected_deta_lead_hard, selector );
  jetHadron::SubtractBackgroundDeta( corrected_deta_sub_hard, selector );
  
  // and normalize
  jetHadron::Normalize1D( corrected_dphi_subtracted_hard, nEventsHard );
  jetHadron::Normalize1D( corrected_dphi_subtracted_sub_hard, nEventsHard );
  jetHadron::Normalize1D( corrected_deta_lead_hard, nEventsHard );
  jetHadron::Normalize1D( corrected_deta_sub_hard, nEventsHard );
  
  // aaaaaand we'll get ahead of ourselves and extract the integrals
  std::vector<std::vector<double> > dphi_lead_bin_int_hard, dphi_sub_bin_int_hard, deta_lead_bin_int_hard, deta_sub_bin_int_hard;
  std::vector<std::vector<double> > dphi_lead_bin_int_hard_err, dphi_sub_bin_int_hard_err, deta_lead_bin_int_hard_err, deta_sub_bin_int_hard_err;
  
  jetHadron::ExtractIntegraldPhi( corrected_dphi_subtracted_hard, dphi_lead_bin_int_hard, dphi_lead_bin_int_hard_err, selector );
  jetHadron::ExtractIntegraldPhi( corrected_dphi_subtracted_sub_hard, dphi_sub_bin_int_hard, dphi_sub_bin_int_hard_err, selector );
  jetHadron::ExtractIntegraldEta( corrected_deta_lead_hard, deta_lead_bin_int_hard, deta_lead_bin_int_hard_err, selector );
  jetHadron::ExtractIntegraldEta( corrected_deta_sub_hard, deta_sub_bin_int_hard, deta_sub_bin_int_hard_err, selector );
  
  // print for my sanity, the histograms overlayed
  jetHadron::PrintPPHardOverlay( corrected_deta_lead_hard[0], corrected_deta_lead[1], outputDirBase + "/DEBUG_PP_OVERLAY_DETA", selector );
  jetHadron::PrintPPHardOverlay( corrected_deta_sub_hard[0], corrected_deta_sub[1], outputDirBase + "/DEBUG_PP_OVERLAY_DETA_SUB", selector );
  jetHadron::PrintPPHardOverlay( corrected_dphi_subtracted_hard[0], corrected_dphi_subtracted[1], outputDirBase + "/DEBUG_PP_OVERLAY_DPHI", selector );
  jetHadron::PrintPPHardOverlay( corrected_dphi_subtracted_sub_hard[0], corrected_dphi_subtracted_sub[1], outputDirBase + "/DEBUG_PP_OVERLAY_DPHI_SUB", selector );
  
  __OUT("reading in the systematic errors")
  //******************************************
  // Now we need to test the systematic errors
  //******************************************
  TFile sysIn( "out/added/pp/trg6/sysv5.root", "READ");
  std::vector<std::vector<TH1F*> > deta_sys, deta_sys_sub, dphi_sys, dphi_sys_sub;
  deta_sys.resize( 1 );
  deta_sys_sub.resize( 1 );
  dphi_sys.resize( 1 );
  dphi_sys_sub.resize( 1 );
  for ( int i = 0; i < corrected_deta_lead[1].size(); ++i ) {
    std::string tmpdEta = "deta_sys_quad_pt_" + patch::to_string(i);
    std::string tmpdEtaSub = "sub_deta_sys_quad_pt_" + patch::to_string(i);
    std::string tmpdPhi = "dphi_sys_quad_pt_" + patch::to_string(i);
    std::string tmpdPhiSub = "sub_dphi_sys_quad_pt_" + patch::to_string(i);
    
    deta_sys[0].push_back( (TH1F*) sysIn.Get( tmpdEta.c_str() ) );
    deta_sys_sub[0].push_back( (TH1F*) sysIn.Get( tmpdEtaSub.c_str() ) );
    dphi_sys[0].push_back( (TH1F*) sysIn.Get( tmpdPhi.c_str() ) );
    dphi_sys_sub[0].push_back( (TH1F*) sysIn.Get( tmpdPhiSub.c_str() ) );
    
  }
  // read in error tree
  TTree* tree = (TTree*) sysIn.Get("yield_err");
  std::vector<std::vector<double> > dphi_yield_sys_rel, dphi_yield_sub_sys_rel, deta_yield_sys_rel, deta_yield_sub_sys_rel;
  dphi_yield_sys_rel.resize(1); dphi_yield_sub_sys_rel.resize(1); deta_yield_sys_rel.resize(1); deta_yield_sub_sys_rel.resize(1);
  double dPhi, dPhi_sub, dEta, dEta_sub;
  tree->SetBranchAddress( "dphi_err", &dPhi );
  tree->SetBranchAddress( "dphi_sub_err", &dPhi_sub );
  tree->SetBranchAddress( "deta_err", &dEta );
  tree->SetBranchAddress( "deta_sub_err", &dEta_sub );
  
  for ( int i = 0; i < selector.nPtBins; ++i ) {
    tree->GetEntry(i);
    dphi_yield_sys_rel[0].push_back(dPhi);
    dphi_yield_sub_sys_rel[0].push_back(dPhi_sub);
    deta_yield_sys_rel[0].push_back(dEta);
    deta_yield_sub_sys_rel[0].push_back(dEta_sub);
  }
  
  // and finally read in the raw correlations from the systematic file,
  // its not used for the analysis, just for diagnostic output &
  // for analysis note
  std::vector<std::vector<TH1D*> > dphi_sys_tow(2), dphi_sys_tow_sub(2), dphi_sys_trk(2), dphi_sys_trk_sub(2);
  std::vector<std::vector<TH1D*> > deta_sys_tow(2), deta_sys_tow_sub(2), deta_sys_trk(2), deta_sys_trk_sub(2);
  
  for ( int i = 0; i < 2; ++i ) {
    for ( int j = 0; j < selector.nPtBins; ++j ) {
      std::string deta_tow_name = "corTowDEta_deta_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
      std::string deta_trk_name = "corTrkDEta_deta_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
      std::string deta_tow_sub_name = "corTowDEtaSub_deta_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
      std::string deta_trk_sub_name = "corTrkDEtaSub_deta_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
      std::string dphi_tow_name = "corTowDPhi_dphi_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
      std::string dphi_trk_name = "corTrkDPhi_dphi_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
      std::string dphi_tow_sub_name = "corTowDPhiSub_dphi_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
      std::string dphi_trk_sub_name = "corTrkDPhiSub_dphi_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
      
      deta_sys_tow[i].push_back( (TH1D*) sysIn.Get( deta_tow_name.c_str() ) );
      deta_sys_trk[i].push_back( (TH1D*) sysIn.Get( deta_trk_name.c_str() ) );
      deta_sys_tow_sub[i].push_back( (TH1D*) sysIn.Get( deta_tow_sub_name.c_str() ) );
      deta_sys_trk_sub[i].push_back( (TH1D*) sysIn.Get( deta_trk_sub_name.c_str() ) );
      dphi_sys_tow[i].push_back( (TH1D*) sysIn.Get( dphi_tow_name.c_str() ) );
      dphi_sys_trk[i].push_back( (TH1D*) sysIn.Get( dphi_trk_name.c_str() ) );
      dphi_sys_tow_sub[i].push_back( (TH1D*) sysIn.Get( dphi_tow_sub_name.c_str() ) );
      dphi_sys_trk_sub[i].push_back( (TH1D*) sysIn.Get( dphi_trk_sub_name.c_str() ) );
      
    }
  }
  
  
  // scale the errors by bin width!
  jetHadron::ScaleErrors( dphi_yield_sys_rel, selector );
  jetHadron::ScaleErrors( dphi_yield_sub_sys_rel, selector );
  jetHadron::ScaleErrors( deta_yield_sys_rel, selector );
  jetHadron::ScaleErrors( deta_yield_sub_sys_rel, selector );
  
  for ( int i = 0; i < dphi_yield_sys_rel[0].size(); ++i ) {
    std::cout<<"bin: "<<i<<std::endl;
    std::cout<<"dphi: "<< dphi_yield_sys_rel[0][i] <<std::endl;
    std::cout<<"dphi_sub: "<< dphi_yield_sub_sys_rel[0][i] << std::endl;
    std::cout<<"deta: "<< deta_yield_sys_rel[0][i]<<std::endl;
    std::cout<<"deta_sub: "<<deta_yield_sub_sys_rel[0][i]<<std::endl;
  }
  
  __OUT("resetting the bin contents for errors")
  // reset sys bin errors to be set centered on pp
  jetHadron::ResetSysBinContent( deta_sys[0], corrected_deta_lead[1], selector );
  jetHadron::ResetSysBinContent( deta_sys_sub[0], corrected_deta_sub[1], selector );
  jetHadron::ResetSysBinContent( dphi_sys[0], corrected_dphi_subtracted[1], selector );
  jetHadron::ResetSysBinContent( dphi_sys_sub[0], corrected_dphi_subtracted_sub[1], selector );
  
  
  __OUT("print them out")
  // and do printouts
  jetHadron::Print1DDPhiHistogramsWithSysErr( corrected_dphi_subtracted[1], dphi_sys[0], selector, outputDirBase+"/dphi_sys_lead", -0.8, 0.8  );
  jetHadron::Print1DDPhiHistogramsWithSysErr( corrected_dphi_subtracted_sub[1], dphi_sys_sub[0], selector, outputDirBase+"/dphi_sys_sub", -0.8 , 0.8  );
  jetHadron::Print1DDEtaHistogramsWithSysErr( corrected_deta_lead[1], deta_sys[0], selector, outputDirBase+"/deta_sys_lead", -0.8 , 0.8  );
  jetHadron::Print1DDEtaHistogramsWithSysErr( corrected_deta_sub[1], deta_sys_sub[0], selector, outputDirBase+"/deta_sys_sub", -0.8 , 0.8  );
  
  // generate some TGraphErrors for those
  __OUT("generate graphs for the systematic errors")
  std::vector<std::vector<double> > dphi_lead_sys_rel_bin_int, dphi_sub_sys_rel_bin_int, deta_lead_sys_rel_bin_int, deta_sub_sys_rel_bin_int;
  std::vector<std::vector<double> > dphi_lead_sys_rel_bin_int_err, dphi_sub_sys_rel_bin_int_err, deta_lead_sys_rel_bin_int_err, deta_sub_sys_rel_bin_int_err;
  
  jetHadron::ExtractIntegraldPhi( dphi_sys, dphi_lead_sys_rel_bin_int, dphi_lead_sys_rel_bin_int_err, selector );
  jetHadron::ExtractIntegraldPhi( dphi_sys_sub, dphi_sub_sys_rel_bin_int, dphi_sub_sys_rel_bin_int_err, selector  );
  jetHadron::ExtractIntegraldEta( deta_sys, deta_lead_sys_rel_bin_int, deta_lead_sys_rel_bin_int_err, selector  );
  jetHadron::ExtractIntegraldEta( deta_sys_sub, deta_sub_sys_rel_bin_int, deta_sub_sys_rel_bin_int_err, selector );
  
  std::vector<TGraphErrors*> dphi_yield_graph_sys_rel = jetHadron::MakeGraphs( ptBinCenters, dphi_lead_sys_rel_bin_int, zeros, dphi_yield_sys_rel, 1, 5, selector, analysisNames, "dphi_sys_rel" );
  std::vector<TGraphErrors*> dphi_sub_yield_graph_sys_rel = jetHadron::MakeGraphs( ptBinCenters, dphi_sub_sys_rel_bin_int, zeros, dphi_yield_sub_sys_rel, 1, 5, selector, analysisNames, "dphi_sub_sys_rel" );
  std::vector<TGraphErrors*> deta_yield_graph_sys_rel = jetHadron::MakeGraphs( ptBinCenters, deta_lead_sys_rel_bin_int, zeros, deta_yield_sys_rel, 1, 5, selector, analysisNames, "deta_sys_rel" );
  std::vector<TGraphErrors*> deta_sub_yield_graph_sys_rel = jetHadron::MakeGraphs( ptBinCenters, deta_sub_sys_rel_bin_int, zeros, deta_yield_sub_sys_rel, 1, 5, selector, analysisNames, "deta_sub_sys_rel" );
  
  __OUT("Generate 5% systematic errors on the yields")
  // generate some 5% error histograms for each of the histograms we use
  std::vector<std::vector<TH1F*> > dphi_yield_err = jetHadron::BuildYieldError( corrected_dphi_subtracted, selector, analysisNames, "dphi_lead_yield_err" );
  std::vector<std::vector<TH1F*> > dphi_sub_yield_err = jetHadron::BuildYieldError( corrected_dphi_subtracted_sub, selector, analysisNames, "dphi_sub_yield_err" );
  std::vector<std::vector<TH1F*> > deta_yield_err = jetHadron::BuildYieldError( corrected_deta_lead, selector, analysisNames, "deta_lead_yield_err" );
  std::vector<std::vector<TH1F*> > deta_sub_yield_err = jetHadron::BuildYieldError( corrected_deta_sub, selector, analysisNames, "deta_sub_yield_err" );
  
  __OUT("now extract the integrals from bin counting")
  // ******************************************
  // now we're testing yields from bin counting
  // ******************************************
  // we'll do the actual one and the one that has the extended range on it
  std::vector<std::vector<double> > dphi_lead_bin_int, dphi_sub_bin_int, deta_lead_bin_int, deta_sub_bin_int;
  std::vector<std::vector<double> > dphi_lead_bin_int_err, dphi_sub_bin_int_err, deta_lead_bin_int_err, deta_sub_bin_int_err;
  std::vector<std::vector<double> > dphi_lead_bin_int_extended, dphi_sub_bin_int_extended, deta_lead_bin_int_extended, deta_sub_bin_int_extended;
  std::vector<std::vector<double> > dphi_lead_bin_int_err_extended, dphi_sub_bin_int_err_extended, deta_lead_bin_int_err_extended, deta_sub_bin_int_err_extended;
  
  // first, the actual results
  jetHadron::ExtractIntegraldPhi( corrected_dphi_subtracted, dphi_lead_bin_int, dphi_lead_bin_int_err, selector );
  jetHadron::ExtractIntegraldPhi( corrected_dphi_subtracted_sub, dphi_sub_bin_int, dphi_sub_bin_int_err, selector );
  jetHadron::ExtractIntegraldEta( corrected_deta_lead, deta_lead_bin_int, deta_lead_bin_int_err, selector );
  jetHadron::ExtractIntegraldEta( corrected_deta_sub, deta_sub_bin_int, deta_sub_bin_int_err, selector );
  
  // now the extended region used for systematics
  jetHadron::ExtractIntegraldPhi( corrected_dphi_subtracted_extended, dphi_lead_bin_int_extended, dphi_lead_bin_int_err_extended, selector );
  jetHadron::ExtractIntegraldPhi( corrected_dphi_subtracted_sub_extended, dphi_sub_bin_int_extended, dphi_sub_bin_int_err_extended, selector );
  jetHadron::ExtractIntegraldEta( corrected_deta_lead_extended, deta_lead_bin_int_extended, deta_lead_bin_int_err_extended, selector );
  jetHadron::ExtractIntegraldEta( corrected_deta_sub_extended, deta_sub_bin_int_extended, deta_sub_bin_int_err_extended, selector );
  
  __OUT("subtracting the normal and extended bin range yields")
  __OUT("using as a systematic error")
  // get difference of yields
  std::vector<std::vector<double> > dphi_subtracted_yield_dif, dphi_sub_subtracted_yield_dif, deta_subtracted_yield_dif, deta_sub_subtracted_yield_dif;
  
  dphi_subtracted_yield_dif.resize( corrected_dphi_subtracted.size() );
  dphi_sub_subtracted_yield_dif.resize( corrected_dphi_subtracted.size() );
  deta_subtracted_yield_dif.resize( corrected_dphi_subtracted.size() );
  deta_sub_subtracted_yield_dif.resize( corrected_dphi_subtracted.size() );
  for ( int i = 0; i < corrected_dphi_subtracted.size(); ++i ) {
    for ( int j = 0; j < corrected_dphi_subtracted[i].size(); ++j ) {
      dphi_subtracted_yield_dif[i].push_back( fabs( dphi_lead_bin_int[i][j] - dphi_lead_bin_int_extended[i][j] ) );
      dphi_sub_subtracted_yield_dif[i].push_back( fabs( dphi_sub_bin_int[i][j] - dphi_sub_bin_int_extended[i][j] ) );
      deta_subtracted_yield_dif[i].push_back( fabs( deta_lead_bin_int[i][j] - deta_lead_bin_int_extended[i][j] ) );
      deta_sub_subtracted_yield_dif[i].push_back( fabs( deta_sub_bin_int[i][j] - deta_sub_bin_int_extended[i][j] ) );
    }
  }
  
  // before subtraction
  std::vector<TGraphErrors*> before_dphi_yield_graph = jetHadron::MakeGraphs( ptBinCenters, dphi_lead_bin_int, zeros, dphi_lead_bin_int_err, 1, 5, selector, analysisNames, "dphi" );
  std::vector<TGraphErrors*> before_dphi_sub_yield_graph = jetHadron::MakeGraphs( ptBinCenters, dphi_sub_bin_int, zeros, dphi_sub_bin_int_err, 1, 5, selector, analysisNames, "dphi_sub" );
  std::vector<TGraphErrors*> before_deta_yield_graph = jetHadron::MakeGraphs( ptBinCenters, deta_lead_bin_int, zeros, deta_lead_bin_int_err, 1, 5, selector, analysisNames, "deta" );
  std::vector<TGraphErrors*> before_deta_sub_yield_graph = jetHadron::MakeGraphs( ptBinCenters, deta_sub_bin_int, zeros, deta_sub_bin_int_err, 1, 5, selector, analysisNames, "deta_sub" );
  
  // just debug for the yield difference between pp with & without hard bkg included
  std::vector<double> dphi_trig_pp_yield_diff{0,0}, dphi_recoil_pp_yield_diff{0,0}, deta_trig_pp_yield_diff{0,0}, deta_recoil_pp_yield_diff{0,0};
  
  __OUT("TESTING THE SUBTRACTION FOR AUAU YIELDS BEING CORRECTED")
  for ( int i = 2; i < dphi_lead_bin_int[0].size(); ++i ) {
    std::cout<<"pt bin: "<<i << std::endl;
    dphi_lead_bin_int[0][i] -= (dphi_lead_bin_int_hard[0][i] - dphi_lead_bin_int[1][i] )  ;
    dphi_trig_pp_yield_diff.push_back( dphi_lead_bin_int_hard[0][i] - dphi_lead_bin_int[1][i] );
    std::cout<<"dphi: "<< dphi_lead_bin_int_hard[0][i] - dphi_lead_bin_int[1][i] << std::endl;
    dphi_sub_bin_int[0][i] -=  ( dphi_sub_bin_int_hard[0][i] - dphi_sub_bin_int[1][i] ) ;
    dphi_recoil_pp_yield_diff.push_back( dphi_sub_bin_int_hard[0][i] - dphi_sub_bin_int[1][i] );
    std::cout<<"dphi sub: "<< dphi_sub_bin_int_hard[0][i] - dphi_sub_bin_int[1][i] << std::endl;
    deta_lead_bin_int[0][i] -= ( deta_lead_bin_int_hard[0][i] - deta_lead_bin_int[1][i] ) ;
    deta_trig_pp_yield_diff.push_back( deta_lead_bin_int_hard[0][i] - deta_lead_bin_int[1][i] );
    std::cout<<"deta: "<< deta_lead_bin_int_hard[0][i] - deta_lead_bin_int[1][i]<<std::endl;
    deta_sub_bin_int[0][i] -=  ( deta_sub_bin_int_hard[0][i] - deta_sub_bin_int[1][i] ) ;
    deta_recoil_pp_yield_diff.push_back( deta_sub_bin_int_hard[0][i] - deta_sub_bin_int[1][i] );
    std::cout<<"deta sub: "<< deta_sub_bin_int_hard[0][i] - deta_sub_bin_int[1][i] << std::endl;
  }
  
  // build the subtracted background graphs first
  // ********************************************
  
  // first step is to get the histogram for random cones
  TH1F* randomConePt = (TH1F*) corrFiles[0]->Get("randomcone");
  randomConePt->Scale( 1.0 / nEvents[0]->GetEntries() );
  randomConePt->RebinX(4);
  
  std::vector<TGraphErrors*> pp_hard_diff_graphs;
  std::vector<std::string> pp_hard_diff_graph_names { "dphi trig", "dphi recoil", "deta trig", "deta recoil" };
  pp_hard_diff_graphs.push_back(jetHadron::MakeGraph( ptBinCenters[1], dphi_trig_pp_yield_diff, zeros[0], zeros[0], 2, 5, selector, "pp_diff", "dphi_pp_diff" ));
  pp_hard_diff_graphs.push_back(jetHadron::MakeGraph( ptBinCenters[1], dphi_recoil_pp_yield_diff, zeros[0], zeros[0], 2, 5, selector, "pp_diff", "dphi_sub_pp_diff" ));
  pp_hard_diff_graphs.push_back(jetHadron::MakeGraph( ptBinCenters[1], deta_trig_pp_yield_diff, zeros[0], zeros[0], 2, 5, selector, "pp_diff", "deta_pp_diff" ));
  pp_hard_diff_graphs.push_back(jetHadron::MakeGraph( ptBinCenters[1], deta_recoil_pp_yield_diff, zeros[0], zeros[0], 2, 5, selector, "pp_diff", "deta_sub_pp_diff" ));
  
  jetHadron::PrintSimpleGraphOverLayWithHistogram( pp_hard_diff_graphs, randomConePt, outputDirBase+"/DEBUG_PP_HARD_DIFF_YIELD", pp_hard_diff_graph_names, "PP_HARD_DIFF" );
  
  // ********************************************

  std::vector<TGraphErrors*> dphi_yield_graph = jetHadron::MakeGraphs( ptBinCenters, dphi_lead_bin_int, zeros, dphi_lead_bin_int_err, 1, 5, selector, analysisNames, "dphi" );
  std::vector<TGraphErrors*> dphi_sub_yield_graph = jetHadron::MakeGraphs( ptBinCenters, dphi_sub_bin_int, zeros, dphi_sub_bin_int_err, 1, 5, selector, analysisNames, "dphi_sub" );
  std::vector<TGraphErrors*> deta_yield_graph = jetHadron::MakeGraphs( ptBinCenters, deta_lead_bin_int, zeros, deta_lead_bin_int_err, 1, 5, selector, analysisNames, "deta" );
  std::vector<TGraphErrors*> deta_sub_yield_graph = jetHadron::MakeGraphs( ptBinCenters, deta_sub_bin_int, zeros, deta_sub_bin_int_err, 1, 5, selector, analysisNames, "deta_sub" );
  
  
  //// annddddddd make some systematic error graphs from the differences
  std::vector<TGraphErrors*> dphi_yield_graph_projection_sys = jetHadron::MakeGraphs( ptBinCenters, dphi_lead_bin_int, zeros, dphi_subtracted_yield_dif, 1, 5, selector, analysisNames, "dphi_proj_sys" );
  std::vector<TGraphErrors*> dphi_sub_yield_graph_projection_sys = jetHadron::MakeGraphs( ptBinCenters, dphi_sub_bin_int, zeros, dphi_sub_subtracted_yield_dif, 1, 5, selector, analysisNames, "dphi_sub_proj_sys" );
  std::vector<TGraphErrors*> deta_yield_graph_projection_sys = jetHadron::MakeGraphs( ptBinCenters, deta_lead_bin_int, zeros, deta_subtracted_yield_dif, 1, 5, selector, analysisNames, "deta_proj_sys" );
  std::vector<TGraphErrors*> deta_sub_yield_graph_projection_sys = jetHadron::MakeGraphs( ptBinCenters, deta_sub_bin_int, zeros, deta_sub_subtracted_yield_dif, 1, 5, selector, analysisNames, "deta_sub_proj_sys" );
  
  
  // and to do that we need some errors as well! We already have relative systematics, but need the yields for the 5% tracking uncertainty
  std::vector<std::vector<double> > dphi_lead_sys_bin_int, dphi_sub_sys_bin_int, deta_lead_sys_bin_int, deta_sub_sys_bin_int;
  std::vector<std::vector<double> > dphi_lead_sys_bin_int_err, dphi_sub_sys_bin_int_err, deta_lead_sys_bin_int_err, deta_sub_sys_bin_int_err;
  
  jetHadron::ExtractIntegraldPhi( dphi_yield_err, dphi_lead_sys_bin_int, dphi_lead_sys_bin_int_err, selector );
  jetHadron::ExtractIntegraldPhi( dphi_sub_yield_err, dphi_sub_sys_bin_int, dphi_sub_sys_bin_int_err, selector );
  jetHadron::ExtractIntegraldEta( deta_yield_err, deta_lead_sys_bin_int, deta_lead_sys_bin_int_err, selector );
  jetHadron::ExtractIntegraldEta( deta_sub_yield_err, deta_sub_sys_bin_int, deta_sub_sys_bin_int_err, selector );

  
  // we are changing from extracting the 5% yield error for the yield plots
  // from the integral of the error histograms, to a raw 5% overall yield on the
  // yield value itself
  dphi_lead_sys_bin_int_err = jetHadron::BuildYieldError( dphi_lead_bin_int, selector );
  dphi_sub_sys_bin_int_err = jetHadron::BuildYieldError( dphi_sub_bin_int, selector );
  deta_lead_sys_bin_int_err = jetHadron::BuildYieldError( deta_lead_bin_int, selector );
  deta_sub_sys_bin_int_err = jetHadron::BuildYieldError( deta_sub_bin_int, selector );
  
  // reset the Au+Au bin centers to the subtracted Au+Au values
  for ( int i = 0; i < dphi_lead_sys_bin_int[0].size(); ++i ) {
    dphi_lead_sys_bin_int[0][i] = dphi_lead_bin_int[0][i];
    dphi_sub_sys_bin_int[0][i] = dphi_sub_bin_int[0][i];
    deta_lead_sys_bin_int[0][i] = deta_lead_bin_int[0][i];
    deta_sub_sys_bin_int[0][i] = deta_sub_bin_int[0][i];
  }
  
  
  std::vector<TGraphErrors*> dphi_yield_sys_graph = jetHadron::MakeGraphs( ptBinCenters, dphi_lead_sys_bin_int, zeros, dphi_lead_sys_bin_int_err, 1, 5, selector, analysisNames, "dphi_sys" );
  std::vector<TGraphErrors*> dphi_sub_yield_sys_graph = jetHadron::MakeGraphs( ptBinCenters, dphi_sub_sys_bin_int, zeros, dphi_sub_sys_bin_int_err, 1, 5, selector, analysisNames, "dphi_sub_sys" );
  std::vector<TGraphErrors*> deta_yield_sys_graph = jetHadron::MakeGraphs( ptBinCenters, deta_lead_sys_bin_int, zeros, deta_lead_sys_bin_int_err, 1, 5, selector, analysisNames, "deta_sys" );
  std::vector<TGraphErrors*> deta_sub_yield_sys_graph = jetHadron::MakeGraphs( ptBinCenters, deta_sub_sys_bin_int, zeros, deta_sub_sys_bin_int_err, 1, 5, selector, analysisNames, "deta_sub_sys" );
  
  
  jetHadron::PrintGraphsWithSystematics( dphi_yield_graph, dphi_yield_sys_graph, dphi_yield_graph_sys_rel, outputDirBase+"/new_trig_yield_dphi", analysisNames, "Trigger Jet Yield #Delta#phi", selector );
  jetHadron::PrintSystematicsOverlayGraph(  dphi_yield_sys_graph[0], nullptr,  outputDirBase+"/new_trig_yield_dphi", "dphi_trigger_auau_err", "p_{T}", "Error", "5% Yield Uncertainty", "" );
  jetHadron::PrintSystematicsOverlayGraph(  dphi_yield_sys_graph[1], dphi_yield_graph_sys_rel[0],  outputDirBase+"/new_trig_yield_dphi", "dphi_trigger_pp_err", "p_{T}", "Error", "5% Yield Uncertainty", "JES uncertainty" );
  
  jetHadron::PrintGraphsWithSystematics( dphi_sub_yield_graph, dphi_sub_yield_sys_graph, dphi_sub_yield_graph_sys_rel, outputDirBase+"/new_recoil_yield_dphi", analysisNames, "Recoil Jet Yield #Delta#phi", selector );
  jetHadron::PrintSystematicsOverlayGraph(  dphi_sub_yield_sys_graph[0], nullptr,  outputDirBase+"/new_recoil_yield_dphi", "dphi_recoil_auau_err", "p_{T}", "Error", "5% Yield Uncertainty", "" );
  jetHadron::PrintSystematicsOverlayGraph(  dphi_sub_yield_sys_graph[1], dphi_sub_yield_graph_sys_rel[0],  outputDirBase+"/new_recoil_yield_dphi", "dphi_recoil_pp_err", "p_{T}", "Error", "5% Yield Uncertainty", "JES uncertainty" );
  
  jetHadron::PrintGraphsWithSystematics( deta_yield_graph, deta_yield_sys_graph, deta_yield_graph_sys_rel, outputDirBase+"/new_trig_yield_deta", analysisNames, "Trigger Jet Yield #Delta#eta", selector );
  jetHadron::PrintSystematicsOverlayGraph(  deta_yield_sys_graph[0], nullptr,  outputDirBase+"/new_trig_yield_deta", "deta_trigger_auau_err", "p_{T}", "Error", "5% Yield Uncertainty", "" );
  jetHadron::PrintSystematicsOverlayGraph(  deta_yield_sys_graph[1], deta_yield_graph_sys_rel[0],  outputDirBase+"/new_trig_yield_deta", "deta_trigger_pp_err", "p_{T}", "Error", "5% Yield Uncertainty", "JES uncertainty" );
  
  jetHadron::PrintGraphsWithSystematics( deta_sub_yield_graph, deta_sub_yield_sys_graph, deta_sub_yield_graph_sys_rel, outputDirBase+"/new_recoil_yield_deta", analysisNames, "Recoil Jet Yield #Delta#eta", selector );
  jetHadron::PrintSystematicsOverlayGraph(  deta_sub_yield_sys_graph[0], nullptr,  outputDirBase+"/new_recoil_yield_deta", "deta_recoil_auau_err", "p_{T}", "Error", "5% Yield Uncertainty", "" );
  jetHadron::PrintSystematicsOverlayGraph(  deta_sub_yield_sys_graph[1], deta_sub_yield_graph_sys_rel[0],  outputDirBase+"/new_recoil_yield_deta", "deta_recoil_pp_err", "p_{T}", "Error", "5% Yield Uncertainty", "JES uncertainty" );
  
  
  std::vector<std::string> phiText, phiTextSub;
  phiText.push_back("|#Delta#eta|<0.45");
  phiTextSub.push_back("|#Delta#eta|<0.45");
  phiText.push_back("trigger jet");
  phiTextSub.push_back("recoil jet");
  std::vector<std::string> etaText, etaTextSub;
  etaText.push_back("|#Delta#phi|<0.71");
  etaTextSub.push_back("|#Delta#phi|<0.71");
  etaText.push_back("trigger jet");
  etaTextSub.push_back("recoil jet");
  
  // check errors on yields
  jetHadron::Print1DDPhiHistogramsWithSysErr( corrected_dphi_subtracted, dphi_yield_err, selector, outputDirBase+"/dphi_yield_err_lead", -0.8, 0.8  );
  jetHadron::Print1DDPhiHistogramsWithSysErr( corrected_dphi_subtracted_sub, dphi_sub_yield_err, selector, outputDirBase+"/dphi_yield_err_sub", -0.8 , 0.8  );
  jetHadron::Print1DDEtaHistogramsWithSysErr( corrected_deta_lead, deta_yield_err, selector, outputDirBase+"/deta_yield_err_lead", -0.8 , 0.8  );
  jetHadron::Print1DDEtaHistogramsWithSysErr( corrected_deta_sub, deta_sub_yield_err, selector, outputDirBase+"/deta_yield_err_sub", -0.8 , 0.8  );
  
  
  // and plot the full histograms
  std::vector<TH1F*> tmp_vec;
  
  jetHadron::FixTheDamnBins( corrected_dphi_subtracted );
  jetHadron::Print1DDPhiHistogramsWithSysErrFull( corrected_dphi_subtracted, dphi_yield_err, dphi_sys[0], selector, outputDirBase+"/new_dphi_yield_err_lead", phiText, -0.8, 0.8  );
  jetHadron::PrintSystematicsOverlay( dphi_yield_err[0], tmp_vec, outputDirBase+"/new_dphi_yield_err_lead_auau", selector, "#Delta#phi", "Error", "5% Yield Error", "", -0.8, 0.8 );
  jetHadron::PrintSystematicsOverlay( dphi_yield_err[1], dphi_sys[0], outputDirBase+"/new_dphi_yield_err_lead_pp", selector, "#Delta#phi", "Error", "5% Yield Error", "Relative Tracking Uncertainty", -0.8, 0.8 );
  
  jetHadron::FixTheDamnBins( corrected_dphi_subtracted_sub );
  jetHadron::Print1DDPhiHistogramsWithSysErrFull( corrected_dphi_subtracted_sub, dphi_sub_yield_err, dphi_sys_sub[0], selector, outputDirBase+"/new_dphi_yield_err_sub", phiTextSub, -0.8 , 0.8  );
  jetHadron::PrintSystematicsOverlay( dphi_sub_yield_err[0], tmp_vec, outputDirBase+"/new_dphi_yield_err_sub_auau", selector, "#Delta#phi", "Error", "5% Yield Error", "", -0.8, 0.8 );
  jetHadron::PrintSystematicsOverlay( dphi_sub_yield_err[1], dphi_sys_sub[0], outputDirBase+"/new_dphi_yield_err_sub_pp", selector, "#Delta#phi", "Error", "5% Yield Error", "Relative Tracking Uncertainty", -0.8, 0.8 );
  
  jetHadron::FixTheDamnBins( corrected_deta_lead );
  jetHadron::Print1DDEtaHistogramsWithSysErrFull( corrected_deta_lead, deta_yield_err, deta_sys[0], selector, outputDirBase+"/new_deta_yield_err_lead", etaText, -0.8 , 0.8  );
  jetHadron::PrintSystematicsOverlay( deta_yield_err[0], tmp_vec, outputDirBase+"/new_deta_yield_err_lead_auau", selector, "#Delta#eta", "Error", "5% Yield Error", "", -0.8, 0.8 );
  jetHadron::PrintSystematicsOverlay( deta_yield_err[1], deta_sys[0], outputDirBase+"/new_deta_yield_err_lead_pp", selector, "#Delta#eta", "Error", "5% Yield Error", "Relative Tracking Uncertainty", -0.8, 0.8 );
  
  jetHadron::FixTheDamnBins( corrected_deta_sub );
  jetHadron::Print1DDEtaHistogramsWithSysErrFull( corrected_deta_sub, deta_sub_yield_err, deta_sys_sub[0], selector, outputDirBase+"/new_deta_yield_err_sub", etaTextSub, -0.8 , 0.8  );
  jetHadron::PrintSystematicsOverlay( deta_sub_yield_err[0], tmp_vec, outputDirBase+"/new_deta_yield_err_sub_auau", selector, "#Delta#eta", "Error", "5% Yield Error", "", -0.8, 0.8 );
  jetHadron::PrintSystematicsOverlay( deta_sub_yield_err[1], deta_sys_sub[0], outputDirBase+"/new_deta_yield_err_sub_pp", selector, "#Delta#eta", "Error", "5% Yield Error", "Relative Tracking Uncertainty", -0.8, 0.8 );
  
  std::vector<std::string> overlayText, overlayTextSub;
  overlayText.push_back("trigger jet");
  overlayTextSub.push_back("recoil jet");

  return 0;
}