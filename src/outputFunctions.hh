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

#ifndef OUTPUTFUNCTIONS_HH
#define OUTPUTFUNCTIONS_HH

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
    
    // pt
    double ptBinEdgeLo[5] = { 0.5, 1.0, 2.0, 3.0, 4.0 };
    double ptBinEdgeHi[5] = { 1.0, 2.0, 3.0, 4.0, 6.0 };
    double ptBinWidth = ( ptHighEdge - ptLowEdge ) / binsPt;
    const int nPtBins = 5;
    
    double ptBinLowEdge( int i ) {if (i < 5 && i >= 0 ) return ( ptBinEdgeLo[i]/ptBinWidth ) + 1; else __ERR("bad pt bin index") return 0;  }
    double ptBinHighEdge( int i ) {if (i < 5 && i >= 0 ) return ( ptBinEdgeHi[i]/ptBinWidth ); else __ERR("bad pt bin index") return 0;  }
    
    // these can be used to help with deta and dphi binning
    unsigned bindEta = 22;
    double dEtaLow = -2.0;
    double dEtaHigh = 2.0;
    unsigned bindPhi = 22;
    double dPhiLow = -pi/2.0;
    double dPhiHigh = 3.0*pi/2.0;
    
    void SetHistogramBins( TH3F* h );
    
    
  };
  
  // Function used to read in histograms from
  // the files passed in - it returns the correlations,
  // and the number of events, and selects using the centralities,
  // vz bin range, and aj ranges passed in via binSelector
  void ReadInFiles(std::vector<TFile*>& filesIn, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& leadingCorrelations, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& subLeadingCorrelations, std::vector<TH3F*>& nEvents, binSelector selector );
  
  
  // Function used to find the weighted center
  // for each pt bin for each file - vector<vector<double> >
  // and also creates pt spectra for each file
  std::vector<std::vector<double> > FindPtBinCenter( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& correlations, std::vector<TH1F*>& ptSpectra, binSelector selector );
  
} // end namespace

#endif
