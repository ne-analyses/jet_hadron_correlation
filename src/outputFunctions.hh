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
    
    double ptBinLowEdge( int i ) {if (i < 5 && i >= 0 ) return ( ptBinEdgeLo[i]/ptBinWidth ) + 1; else __ERR("bad pt bin index") return 0;  }
    double ptBinHighEdge( int i ) {if (i < 5 && i >= 0 ) return ( ptBinEdgeHi[i]/ptBinWidth ); else __ERR("bad pt bin index") return 0;  }
    
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
    
    // used to select the near size correlation for dEta
    double eta_projection_phi_bound_low = -1.0;
    double eta_projection_phi_bound_high = 1.0;
    
    // used when projecting the mixed events to get  dPhi
    double phi_projection_eta_bound_low = -1.2;
    double phi_projection_eta_bound_high = 1.2;
    
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
  int ReadInFiles(std::vector<TFile*>& filesIn, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& leadingCorrelations, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& subLeadingCorrelations, std::vector<TH3F*>& nEvents, binSelector selector );
  int ReadInFilesMix(std::vector<TFile*>& filesIn, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& leadingMix, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& subLeadingMix, std::vector<TH3F*>& nEvents, binSelector selector );
  
  
  // Function used to find the weighted center
  // for each pt bin for each file - vector<vector<double> >
  // and also creates pt spectra for each file
  std::vector<std::vector<double> > FindPtBinCenter( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& correlations, std::vector<TH1F*>& ptSpectra, binSelector selector );
  
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
  // projection )
  std::vector<std::vector<TH1F*> > ProjectDphi( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, std::string uniqueID = "", bool restrictDeta = false );
  std::vector<std::vector<TH1F*> > ProjectDphiNearMinusFar( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, double* edges, std::string uniqueID = "", bool restrictDeta = false);
  void ProjectDphiNearMinusFar( std::vector<std::vector<TH2F*> >& correlation2d, std::vector<std::vector<TH1F*> >& near, std::vector<std::vector<TH1F*> >& far, binSelector selector, double* edges, std::string uniqueID = "", bool restrictDeta = false);
  std::vector<std::vector<TH1F*> > ProjectDeta( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, std::string uniqueID = "", bool restrictDphi = false );
  
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
  void ExtractIntegral( std::vector<std::vector<TH1F*> > histograms, std::vector<std::vector<double> > integrals, std::vector<std::vector<double> > errors, binSelector selector, double lowEdge, double highEdge );
  
  
  // Methods for Printing out and Saving
  // ***********************************
  
  void FindGood1DUserRange( std::vector<TH1F*> histograms, double& max, double& min );
  
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
  
} // end namespace

#endif
