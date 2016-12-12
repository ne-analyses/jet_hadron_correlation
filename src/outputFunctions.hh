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
    unsigned ajLow = 0;
    unsigned ajHigh = binsAj-1;
    unsigned centLow = 6;
    unsigned centHigh = 8;
    unsigned vzLow = 0;
    unsigned vzHigh = binsVz-1;
    
    // histogram bins
    double ptBinWidth = ( ptHighEdge - ptLowEdge ) / binsPt;
    const int nPtBins = 5;
    double histogramPtLo[nPtBins] = { (0.5/ptBinWidth)+1, (1.0/ptBinWidth)+1, (2.0/ptBinWidth)+1, (3.0/ptBinWidth)+1, (4.0/ptBinWidth)+1 };
    double histogramPtHi[nPtBins] = { (1.0/ptBinWidth), (2.0/ptBinWidth), (3.0/ptBinWidth), (4.0/ptBinWidth), (6.0/ptBinWidth) };
    

  };
  
  // Function used to read in histograms from
  // the files passed in - it returns the correlations,
  // and the number of events, and selects using the centralities,
  // vz bin range, and aj ranges passed in via binSelector
  void ReadInFiles(std::vector<TFile*>& filesIn, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& leadingCorrelations, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& subLeadingCorrelations, std::vector<TH3F*>& nEvents, binSelector selector );
  
  
} // end namespace

#endif
