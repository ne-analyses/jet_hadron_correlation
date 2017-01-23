// Basic functions used in the
// jet-hadron correlation output
// workflow histogram adding/dividing,
// drawing, etc

// Nick Elsey

// include the constants from the
// correlation analysis
#include "corrParameters.hh"

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
#include "TLegend.h"
#include "TGraphErrors.h"

#ifndef OUTPUTFUNCTIONS_HH
#define OUTPUTFUNCTIONS_HH

namespace patch {
  template < typename T > std::string to_string( const T& n );
}

// first!
// used to clean up these vectors of large amounts of histograms
// mainly for my sanity but also to reduce chance of errors...
template <class T>
void ClearHistograms( std::vector<std::vector<T*> >& histograms )  {
  for ( int i = 0; i < histograms.size(); ++i ) {
    for ( int j = 0; j < histograms[0].size(); ++j ) {
      
      if ( histograms[i][j] )
        delete histograms[i][j];
      
    }
  }
  histograms.resize(0);
}

template <class T>
void ClearHistograms( std::vector<std::vector<std::vector<T*> > >& histograms ) {
  for ( int i = 0; i < histograms.size(); ++i ) {
    for ( int j = 0; j < histograms[i].size(); ++j ) {
      for ( int k = 0; k < histograms[i][j].size(); ++k ) {
        if ( histograms[i][j][k] )
          delete histograms[i][j][k];
        
      }
    }
  }
  histograms.resize(0);
}

template <class T>
void ClearHistograms( std::vector<std::vector<std::vector<std::vector<T*> > > >& histograms ) {
  for ( int i = 0; i < histograms.size(); ++i ) {
    for ( int j = 0; j < histograms[i].size(); ++j ) {
      for ( int k = 0; k < histograms[i][j].size(); ++k ) {
        for ( int l = 0; l < histograms[i][j][k].size(); ++l ) {
          if ( histograms[i][j][k][l] )
            delete histograms[i][j][k][l];
          
        }
      }
    }
  }
  histograms.resize(0);
}


namespace jetHadron {
  
  
  // Used for reading in histograms, makes things a little simpler
  // Default settings are for all Aj bins, 0-20% centrality, and all Vz bins
  // also builds the proper variables for histogram pt projections
  struct binSelector {
    
    // used to select which histograms to read in
    unsigned centLow = 6;
    unsigned centHigh = 8;
    unsigned vzLow = 0;
    unsigned vzHigh = binsVz-1;
    unsigned ajLow = 0;
    unsigned ajHigh = binsAj-1;
    
    // histogram bins
    std::string analysisStrings[2] = { "Au+Au 0-20%", "P+P" };
    
    // pt with low bin
    double ptBinEdgeLo[6] = { 0.5, 1.0, 2.0, 3.0, 4.0, 6.0 };
    double ptBinEdgeHi[6] = { 1.0, 2.0, 3.0, 4.0, 6.0, 10.0 };
    double ptBinWidth = ( ptHighEdge - ptLowEdge ) / binsPt;
    std::string ptBinString[6] = { "0.5 < p_{T} < 1.0", "1.0 < p_{T} < 2.0", "2.0 < p_{T} < 3.0", "3.0 < p_{T} < 4.0", "4.0 < p_{T} < 6.0", "6.0 < p_{T} < 10.0" };
    std::string ptBinStringMix[3] = { "0.5 < p_{T} < 1.0", "1.0 < p_{T} < 2.0", "2.0 < p_{T}" };
    
    const int nPtBins = 6;
    
    // without low bin
    //double ptBinEdgeLo[4] = { 1.0, 2.0, 3.0, 4.0 };
    //double ptBinEdgeHi[4] = { 2.0, 3.0, 4.0, 6.0 };
    //double ptBinWidth = ( ptHighEdge - ptLowEdge ) / binsPt;
    //std::string ptBinString[4] = { "1.0 < p_{T} < 2.0", "2.0 < p_{T} < 3.0", "3.0 < p_{T} < 4.0", "4.0 < p_{T} < 6.0" };
    //const int nPtBins = 4;
    
    double GetPtBinWidth( int i ) { return ptBinEdgeHi[i] - ptBinEdgeLo[i]; }
    
    double ptBinLowEdge( int i ) {if (i < 6 && i >= 0 ) return ( ptBinEdgeLo[i]/ptBinWidth ) + 1; else __ERR("bad pt bin index") return 0;  }
    double ptBinHighEdge( int i ) {if (i < 6 && i >= 0 ) return ( ptBinEdgeHi[i]/ptBinWidth ); else __ERR("bad pt bin index") return 0;  }
    
    // these can be used to help with deta and dphi binning
    unsigned bindEta = 22;
    double dEtaLow = -2.0;
    double dEtaHigh = 2.0;
    double dEtaAcceptanceLow = -1.6;
    double dEtaAcceptanceHigh = 1.6;
    unsigned bindPhi = 22;
    double dPhiLow = -pi/2.0;
    double dPhiHigh = 3.0*pi/2.0;
    
    // if not using the defaults, this can set them all
    // from one of the initial histograms
    void SetHistogramBins( TH2F* h );
    
    // setting ranges for the histogram
    // if the radius is NOT 0.4, this needs to be used...
    void ChangeRadius( double R = 0.4);
    
    // for fitting dEta, the fit range
    const double eta_fit_low_edge = -1.2;
    const double eta_fit_high_edge = -1.0*eta_fit_low_edge;
    
    // and for dPhi without near - far correction
    const double phi_fit_low_edge = - pi / 2.0;
    const double phi_fit_high_edge = 3.0 * pi / 2.0;
    
    // and for dphi with the correction
    const double phi_corrected_fit_low_edge = - 1.2;
    const double phi_corrected_fit_high_edge = 1.2;
    
    // used to select the near side correlation for dEta
    const double eta_projection_phi_bound_low = -0.6;
    const double eta_projection_phi_bound_high = -1.0*eta_projection_phi_bound_low;
    
    // used to select the near side correlation for dEta
    // with a wider range for our systematic errors
    const double eta_projection_phi_bound_low_extended = -0.9;
    const double eta_projection_phi_bound_high_extended = -1.0*eta_projection_phi_bound_low_extended;
    
    // used when projecting the 2d to get  dPhi
    const double phi_projection_eta_bound_low = -0.45;
    const double phi_projection_eta_bound_high = -1.0*phi_projection_eta_bound_low;
    
    // used when projecting the 2d to get  dPhi
    // with a wider range for our systematic errors
    const double phi_projection_eta_bound_low_extended = -0.5;
    const double phi_projection_eta_bound_high_extended = -1.0*phi_projection_eta_bound_low_extended;
    
    // ranges used for integration
    // corresponds to the projection range in the other dimension
    // when taking 1D projections
    const double phi_projection_integral_range_low = eta_projection_phi_bound_low;
    const double phi_projection_integral_range_high = eta_projection_phi_bound_high;;
    
    const double eta_projection_integral_range_low = phi_projection_eta_bound_low;
    const double eta_projection_integral_range_high = phi_projection_eta_bound_high;
    
    const double phi_projection_subtraction_regions[4] = { -0.99, phi_projection_eta_bound_low, phi_projection_eta_bound_high, 0.99 };
    const double phi_projection_subtraction_regions_extended[4] = { -0.99, phi_projection_eta_bound_low_extended, phi_projection_eta_bound_high_extended, 0.99 };
    
    // used when doing near - far subtraction for flow in dPhi
    double restricted_near_phi_projection_eta_bound_low() { return phi_projection_eta_bound_low / 2.0; }
    double restricted_near_phi_projection_eta_bound_high() { return phi_projection_eta_bound_high / 2.0; }
    
    // This takes the entire acceptance instead of the event mixing accepted region
    double near_phi_projection_eta_bound_low() { return dEtaAcceptanceLow / 2.0; }
    double near_phi_projection_eta_bound_high() { return dEtaAcceptanceHigh / 2.0; }
    
    
  };
  
  // Function used to read in histograms from
  // the files passed in - it returns the correlations,
  // and the number of events, and selects using the centralities,
  // vz bin range, and aj ranges passed in via binSelector
  int ReadInFiles(std::vector<TFile*>& filesIn, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& leadingCorrelations, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& subLeadingCorrelations, std::vector<TH3F*>& nEvents, binSelector selector, std::string uniqueID = "" );
  int ReadInFilesMix(std::vector<TFile*>& filesIn, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& leadingMix, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& subLeadingMix, std::vector<TH3F*>& nEvents, binSelector selector, std::string uniqueID = "" );
  
  
  // Function used to find the weighted center
  // for each pt bin for each file - vector<vector<double> >
  // and also creates pt spectra for each file
  std::vector<std::vector<double> > FindPtBinCenter( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& correlations, std::vector<TH1F*>& ptSpectra, binSelector selector, std::string uniqueID = "" );
  
  // Functions to project out the Aj dependence -
  // can either produce a single, Aj independent bin
  // or splits on an ajbin
  void BuildSingleCorrelation( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& correlations, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& reducedCorrelations, binSelector selector, std::string uniqueID = "" );
  void BuildAjSplitCorrelation( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& correlations, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& reducedCorrelationsHigh, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& reducedCorrelationsLow, binSelector selector, int ajBinSplit, std::string uniqueID = "" );
  
  // Averages over all vz and centralities
  // to show uncorrected signals
  std::vector<std::vector<TH2F*> > AverageCorrelations( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& correlations, binSelector selector, std::string uniqueID = "" );
  
  // Used to recombine Aj and split in pt
  // to give 2D projections we can turn use
  // to correct the correlations
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > BuildMixedEvents( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& mixedEvents, binSelector selector, std::string uniqueID = "" );
  
  // Used to average the mixed event data to help
  // with the lower statistics
  std::vector<std::vector<TH2F*> > RecombineMixedEvents( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& mixedEvents, binSelector selector, std::string uniqueID = "" );
  
  // Used to partially recomine over Vz bins
  // but leaves centrality untouched
  std::vector<std::vector<std::vector<TH2F*> > > PartialRecombineMixedEvents( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& mixedEvents, binSelector selector, std::string uniqueID = "" );
  
  // Used to normalize mixed event histograms so
  // that the maximum bin content = 1
  // version for both the independent mixed events and the weighed averages
  void ScaleMixedEvents( std::vector<std::vector<TH2F*> >& mixedEvents );
  void ScaleMixedEvents( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& mixedEvents );
  
  // Used to perform the mixed event division
  // And add up everything into a 2D structure
  // only keeping differntial in file and Pt
  // Has a version for both the averaged and non
  // averaged event mixing
  std::vector<std::vector<TH2F*> > EventMixingCorrection( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& correlations, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& mixedEvents, binSelector selector, std::string uniqueID = "" );
  std::vector<std::vector<TH2F*> > EventMixingCorrection( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& correlations, std::vector<std::vector<std::vector<TH2F*> > >& mixedEvents, binSelector selector, std::string uniqueID = "" );
  std::vector<std::vector<TH2F*> > EventMixingCorrection( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& correlations, std::vector<std::vector<TH2F*> >& mixedEvents, binSelector selector, std::string uniqueID = "" );
  
  
  // Used to extract 1D projections from
  // the 2D histograms - allows for setting
  // ranges for the projection ( e.g. projecting
  // only the near side of the dPhi range in a dEta
  // projection ) via binSelector
  std::vector<std::vector<TH1F*> > ProjectDphi( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, std::string uniqueID = "" );
  std::vector<std::vector<TH1F*> > ProjectDphiNearMinusFar( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, std::string uniqueID = "" );
  void ProjectDphiNearMinusFar( std::vector<std::vector<TH2F*> >& correlation2d, std::vector<std::vector<TH1F*> >& near, std::vector<std::vector<TH1F*> >& far, binSelector selector, std::string uniqueID = "" );
  
  std::vector<std::vector<TH1F*> > ProjectDeta( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, std::string uniqueID = "" );
  
  // extended range for our systematic errors
  std::vector<std::vector<TH1F*> > ProjectDetaExtended( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, std::string uniqueID = "" );
  std::vector<std::vector<TH1F*> > ProjectDphiNearMinusFarExtended( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, std::string uniqueID = "" );
  
  // Normalizes 1D histograms based on what
  // type of analysis they are
  void Normalize1D( std::vector<std::vector<TH1F*> >& histograms, std::vector<TH3F*>& nEvents );
  void Normalize1DAjSplit( std::vector<std::vector<TH1F*> >& histograms, std::vector<TH3F*>& nEvents, int ajBinLow, int ajBinHigh );
  
  // Used to subtract one set of 1D histograms
  // from another - does not do background sub
  // or anything like that
  std::vector<std::vector<TH1F*> > Subtract1D( std::vector<std::vector<TH1F*> >& base, std::vector<std::vector<TH1F*> >& subtraction, std::string = "" );
  
  // Used to subtract background from each histogram
  void SubtractBackgroundDeta( std::vector<std::vector<TH1F*> >& histograms, binSelector selector );
  void SubtractBackgroundDphi( std::vector<std::vector<TH1F*> >& histograms, binSelector selector );
  
  // Used to fit each histogram
  std::vector<std::vector<TF1*> > FitDeta( std::vector<std::vector<TH1F*> >& histograms, binSelector selector, std::string uniqueID = "" );
  std::vector<std::vector<TF1*> > FitDphi( std::vector<std::vector<TH1F*> >& histograms, binSelector selector, std::string uniqueID = "" );
  std::vector<std::vector<TF1*> > FitDphiRestricted( std::vector<std::vector<TH1F*> >& histograms, binSelector selector, std::string uniqueID = "" );
  
  // Extracts the yield and errors from the fits
  // only extracts for near side so can be used
  // for both dphi and deta
  void ExtractFitVals( std::vector<std::vector<TF1*> >& fits, std::vector<std::vector<double> >& yields, std::vector<std::vector<double> >& widths, std::vector<std::vector<double> >& normError, std::vector<std::vector<double> >& widthError, binSelector selector );
  
  // Used to get the integrals of the
  // histograms, and errors on the integrals
  void ExtractIntegraldPhi( std::vector<std::vector<TH1F*> >& histograms, std::vector<std::vector<double> >& integrals, std::vector<std::vector<double> >& errors, binSelector selector );
  void ExtractIntegraldEta( std::vector<std::vector<TH1F*> >& histograms, std::vector<std::vector<double> >& integrals, std::vector<std::vector<double> >& errors, binSelector selector );
  
  // testing function to fix histogram
  // bins not being drawn if the content is small
  void FixTheDamnBins( std::vector<std::vector<TH1F*> >& histograms );
  void FixTheDamnBins( std::vector<TH1F*>& histograms );
  
  std::vector<TGraphErrors*> MakeGraphs( std::vector<std::vector<double> >& x, std::vector<std::vector<double> >& y, std::vector<std::vector<double> >& x_err, std::vector<std::vector<double> >& y_err, int ptBinLow, int ptBinHigh, binSelector selector, std::vector<std::string> analysisName, std::string uniqueID = "" );
  
  // ***************************************
  // these are used for building uncertainty
  // bands for the pp data
  // ****************************************
  
  // Used to create systematic uncertainty bands
  // from varying tower energy scale / TPC tracking variables
  std::vector<TH1F*> BuildSystematicHistogram( std::vector<TH1F*>& upper, std::vector<TH1F*>& lower, binSelector selector, std::string uniqueID = "" );
  
  // Used to add systematic errors in quadrature
  std::vector<TH1F*> AddInQuadrature( std::vector<TH1F*> upper, std::vector<TH1F*> lower, binSelector selector, std::string uniqueID = "" );
  
  // and using only numbers
  std::vector<double> AddInQuadrature( std::vector<double> upper, std::vector<double> lower );
  
  
  // used to make 5% errors on yields due to tracking
  std::vector<std::vector<TH1F*> > BuildYieldError( std::vector<std::vector<TH1F*> > histograms, binSelector selector, std::vector<std::string> analysisName, std::string uniqueID = "" );
  // used to make full 5% errors on yields
  std::vector<std::vector<double> > BuildYieldError ( std::vector<std::vector<double> > &yields, binSelector selector );
  
   // testing function to reset the bin contents to those of the histogram
  void ResetSysBinContent( std::vector<TH1F*>& errors, std::vector<TH1F*>& histograms, binSelector selector );
  
  // and used to scale errors by the pt bin width
  void ScaleErrors( std::vector<std::vector<double> > errors, binSelector selector );
  
  // used to extract only the yields, dont need the errors
  std::vector<std::vector<double> > OnlyYieldsEta( std::vector<std::vector<TH1F*> >& histograms, binSelector selector );
  std::vector<std::vector<double> > OnlyYieldsPhi( std::vector<std::vector<TH1F*> >& histograms, binSelector selector );
  
  // used to get the difference between two sets of values... specific use
  std::vector<double> GetDifference( std::vector<std::vector<double> >& yields );
  
  // Methods for Printing out and Saving
  // ***********************************
  
  void FindGood1DUserRange( std::vector<TH1F*> histograms, double& max, double& min );
  void FindGood1DUserRange( std::vector<TH1F*> histograms, double& max, double& min, double xMax, double xMin );
  
  // Used to print out and save the 2D prots ( correlations, mixed events )
  void Print2DHistograms( std::vector<TH2F*>& histograms, std::string outputDir, std::string analysisName, binSelector selector );
  void Print2DHistogramsMixing( std::vector<TH2F*>& histograms, std::string outputDir, std::string analysisName, binSelector selector );
  void Print2DHistogramsEtaRestricted( std::vector<TH2F*>& histograms, std::string outputDir, std::string analysisName, binSelector selector );
  
  void Print1DHistogramsDphi( std::vector<TH1F*>& histograms, std::string outputDir, std::string analysisName, binSelector selector );
  void Print1DHistogramsOverlayedDphi( std::vector<std::vector<TH1F*> >& histograms, std::string outputDir, std::vector<std::string> analysisName, binSelector selector );
  void Print1DHistogramsOverlayedDphiWFit( std::vector<std::vector<TH1F*> >& histograms, std::vector<std::vector<TF1*> >& fits, std::string outputDir, std::vector<std::string> analysisName, binSelector selector );
  void Print1DHistogramsOverlayedDphiWFitRestricted( std::vector<std::vector<TH1F*> >& histograms, std::vector<std::vector<TF1*> >& fits, std::string outputDir, std::vector<std::string> analysisName, binSelector selector );
  
  void Print1DHistogramsDeta( std::vector<TH1F*>& histograms, std::string outputDir, std::string analysisName, binSelector selector );
  void Print1DHistogramsOverlayedDeta( std::vector<std::vector<TH1F*> >& histograms, std::string outputDir, std::vector<std::string> analysisName, binSelector selector );
  void Print1DHistogramsOverlayedDetaWFit( std::vector<std::vector<TH1F*> >& histograms, std::vector<std::vector<TF1*> >& fits, std::string outputDir, std::vector<std::string> analysisName, binSelector selector );
  void Print1DHistogramsOverlayedDetaWFitRestricted( std::vector<std::vector<TH1F*> >& histograms, std::vector<std::vector<TF1*> >& fits, std::string outputDir, std::vector<std::string> analysisName, binSelector selector );
  
  void Print1DHistogramsOverlayedDphiOther( std::vector<TH1F*>& histograms, std::vector<TH1F*>& histograms2, std::string outputDir, std::string analysisName1, std::string analysisName2, binSelector selector );
  
  void PrintGraphWithErrors( std::vector<std::vector<double> > x, std::vector<std::vector<double> > y, std::vector<std::vector<double> > x_err, std::vector<std::vector<double> > y_err, std::string outputDir, std::vector<std::string> analysisName, std::string title, binSelector selector, const int pt_min, const int pt_max );
  
  // printing with errors
  void Print1DDPhiHistogramsWithSysErr( std::vector<TH1F*>& histograms, std::vector<TH1F*>& errors, binSelector selector, std::string outputDir, double rangeLow, double rangeHigh  );
  void Print1DDEtaHistogramsWithSysErr( std::vector<TH1F*>& histograms, std::vector<TH1F*>& errors, binSelector selector, std::string outputDir, double rangeLow, double rangeHigh  );
  void Print1DDPhiHistogramsWithSysErr( std::vector<std::vector<TH1F*> >& histograms, std::vector<std::vector<TH1F*> >& errors, binSelector selector, std::string outputDir, double rangeLow, double rangeHigh  );
  void Print1DDEtaHistogramsWithSysErr( std::vector<std::vector<TH1F*> >& histograms, std::vector<std::vector<TH1F*> >& errors, binSelector selector, std::string outputDir, double rangeLow, double rangeHigh  );
  
  // and the plotting with full set of errors
  void Print1DDPhiHistogramsWithSysErrFull( std::vector<std::vector<TH1F*> >& histograms, std::vector<std::vector<TH1F*> >& errors, std::vector<TH1F*>& errors2, binSelector selector, std::string outputDir, std::vector<std::string> text, double rangeLow, double rangeHigh  );
  void Print1DDEtaHistogramsWithSysErrFull( std::vector<std::vector<TH1F*> >& histograms, std::vector<std::vector<TH1F*> >& errors, std::vector<TH1F*>& errors2, binSelector selector, std::string outputDir, std::vector<std::string> text, double rangeLow, double rangeHigh  );
  
  // and a function to print out the near/far dphi correlations
  void PrintNearFarDPhiCorrelations( std::vector<TH1F*> hist1, std::vector<TH1F*> hist2, binSelector selector, std::string outputDir, std::vector<std::string> text, double rangeLow, double rangeHigh );
  
  // printing some graphs with some systematic errors as well
  void PrintGraphsWithSystematics( std::vector<TGraphErrors*>& graphs, std::vector<TGraphErrors*>& sys1, std::vector<TGraphErrors*> sys2, std::string outputDir, std::vector<std::string> analysisName, std::string title, binSelector selector );
  void PrintGraphsWithSystematics( std::vector<TGraphErrors*>& graphs, std::vector<TGraphErrors*>& sys1, std::vector<TGraphErrors*> sys2, std::vector<TGraphErrors*> sys3, std::string outputDir, std::vector<std::string> analysisName, std::string title, binSelector selector );
  
} // end namespace

#endif
