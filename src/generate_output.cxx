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

// Because the grid is ridiculous and doesnt
// Have std::to_string
namespace patch
{
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}

// list all input files as arguments -
// 0 = corr1
// 1 = mix1
// 2 = analysis1 identifying string
// 3 = corr2
// 4 = mix2
// .......

int main( int argc, const char** argv) {
  
  // Pt bins
  const int nPtBins = 5;
  double ptBinLo[nPtBins] = { 3, 5, 9, 13, 17 };
  double ptBinHi[nPtBins] = { 4, 8, 12, 16, 24 };
  std::string ptBinString[nPtBins] = { "0.5-1.0", "1.0-2.0", "2.0-3.0", "3.0-4.0", "4.0-6.0" };
  
  // analysis names
  std::vector<std::string> defaultCorrNames;
  defaultCorrNames.resize(4);
  defaultCorrNames[0] = "Dijet";
  defaultCorrNames[1] = "10 < Jet < 15";
  defaultCorrNames[2] = "15 < Jet < 20";
  defaultCorrNames[3] = "Jet > 20";
  
  
  // First check to make sure we're located properly
  std::string currentDirectory = corrAnalysis::getPWD( );
  
  // If we arent in the analysis directory, exit
  if ( !(corrAnalysis::HasEnding ( currentDirectory, "jet_hadron_corr" ) || corrAnalysis::HasEnding ( currentDirectory, "jet_hadron_correlation" )) ) {
    std::cerr << "Error: Need to be in jet_hadron_corr directory" << std::endl;
    return -1;
  }
  
  
  // files and naming
  TFile** corrFiles;
  TFile** mixFiles;
  std::vector<std::string> analysisNames;
  
  switch ( argc ) {
    case 1: { // Default case
      __OUT( "Using Default Settings" )
      corrFiles = new TFile*[4];
      mixFiles = new TFile*[4];
      
      // default files
      corrFiles[0] = new TFile( "out/tmp/dijet_corr.root", "READ" );
      corrFiles[1] = new TFile( "out/tmp/jet10_corr.root", "READ" );
      corrFiles[2] = new TFile( "out/tmp/jet15_corr.root", "READ" );
      corrFiles[3] = new TFile( "out/tmp/jet20_corr.root", "READ" );
      mixFiles[0] = new TFile( "out/tmp/dijet_mix.root", "READ" );
      mixFiles[1] = new TFile( "out/tmp/jet10_mix.root", "READ" );
      mixFiles[2] = new TFile( "out/tmp/jet15_mix.root", "READ" );
      mixFiles[3] = new TFile( "out/tmp/jet20_mix.root", "READ" );
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
      
      corrFiles = new TFile*[ nCorrFiles ];
      mixFiles = new TFile*[ nCorrFiles ];
      
      for ( int i = 0; i < nCorrFiles; ++i ) {
        corrFiles[i] = new TFile( arguments[ 3*i ].c_str(), "READ" );
        mixFiles[i] = new TFile( arguments[ (3*i)+1 ].c_str(), "READ" );
        analysisNames[i] = arguments[ (3*i)+2 ];
      }
    }
  }
  
  int nFiles = analysisNames.size();
  
  // Load in the histograms
  TH2D* nEvents[ nFiles ];
  TH1D* hVz[ nFiles ];
  TH3D* corrHist[ nFiles ];
  TH3D* mixHist[ nFiles ];
  std::vector<std::vector<std::vector<TH3D*> > > corrCentVz;
  std::vector<std::vector<std::vector<TH3D*> > > subCentVz;
  std::vector<std::vector<std::vector<TH3D*> > > mixCentVz;
  std::vector<std::vector<std::vector<TH3D*> > > mixSubCentVz;
  corrCentVz.resize( nFiles );
  subCentVz.resize( nFiles );
  mixCentVz.resize( nFiles );
  mixSubCentVz.resize( nFiles );
  for ( int i = 0; i < nFiles; ++i ) {
    corrCentVz[i].resize( corrAnalysis::binsCentrality );
    subCentVz[i].resize( corrAnalysis::binsCentrality );
    mixCentVz[i].resize( corrAnalysis::binsCentrality );
    mixSubCentVz[i].resize( corrAnalysis::binsCentrality );
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      corrCentVz[i][j].resize( corrAnalysis::binsVz );
      subCentVz[i][j].resize( corrAnalysis::binsVz );
      mixCentVz[i][j].resize( corrAnalysis::binsVz );
      mixSubCentVz[i][j].resize( corrAnalysis::binsVz );
    }
  }
  
  for ( int i = 0; i < nFiles; ++i ) {
    
    std::string neventsBaseName = "nevents_"; neventsBaseName += analysisNames[i];
    std::string hvzBaseName = "hvz_"; hvzBaseName += analysisNames[i];
    std::string corrhistBaseName = "corrHist_"; corrhistBaseName += analysisNames[i];
    std::string mixhistBaseName = "mixHist_"; mixhistBaseName += analysisNames[i];

    nEvents[i] = (TH2D*) corrFiles[i]->Get( "nevents" );
    nEvents[i]->SetName( neventsBaseName.c_str() );
    hVz[i] = (TH1D*) corrFiles[i]->Get( "vzdist" );
    hVz[i]->SetName( hvzBaseName.c_str() );
    corrHist[i] = (TH3D*) corrFiles[i]->Get( "leadjetcorr" );
    corrHist[i]->SetName( corrhistBaseName.c_str() );
    mixHist[i] = (TH3D*) mixFiles[i]->Get( "leadjetcorr" );
    mixHist[i]->SetName( mixhistBaseName.c_str() );
    
    // pull in the cent/vz diffentiated histograms
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j )
      for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {

        // make the initial name
        std::string corrDifInitName = "lead_cent_"; corrDifInitName += patch::to_string(j);
        std::string subDifInitName = "sub_cent_"; subDifInitName += patch::to_string(j);
        std::string mixDifInitName = "mix_lead_cent_"; mixDifInitName += patch::to_string(j);
        std::string mixSubDifInitName = "mix_sub_cent_";
        mixSubDifInitName += patch::to_string(j);
        
        corrDifInitName += "_vz_"; corrDifInitName += patch::to_string(k);
        subDifInitName += "_vz_"; subDifInitName += patch::to_string(k);
        mixDifInitName += "_vz_"; mixDifInitName += patch::to_string(k);
        mixSubDifInitName += "_vz_"; mixSubDifInitName += patch::to_string(k);
        
        
        // make the new histogram name
        std::string corrDifBaseName = "corr_file_"; corrDifBaseName += patch::to_string(i);
        std::string subDifBaseName = "sub_file_"; subDifBaseName += patch::to_string(i);
        std::string mixDifBaseName = "mix_file_"; mixDifBaseName += patch::to_string(i);
        std::string mixSubDifBaseName = "mix_file_"; mixSubDifBaseName += patch::to_string(i);
        
        corrDifBaseName += "_cent_"; corrDifBaseName += patch::to_string(j);
        corrDifBaseName += "_vz_"; corrDifBaseName += patch::to_string(k);
        
        subDifBaseName += "_cent_"; subDifBaseName += patch::to_string(j);
        subDifBaseName += "_vz_"; subDifBaseName += patch::to_string(k);
        
        mixDifBaseName += "_cent_"; mixDifBaseName += patch::to_string(j);
        mixDifBaseName += "_vz_"; mixDifBaseName += patch::to_string(k);
        
        mixSubDifBaseName += "_cent_"; mixSubDifBaseName += patch::to_string(j);
        mixSubDifBaseName += "_vz_"; mixSubDifBaseName += patch::to_string(k);
        
        // get the histograms
        corrCentVz[i][j][k] = (TH3D*) corrFiles[i]->Get( corrDifInitName.c_str() );
        corrCentVz[i][j][k]->SetName( corrDifBaseName.c_str() );
        subCentVz[i][j][k] = (TH3D*) corrFiles[i]->Get( subDifInitName.c_str() );
        subCentVz[i][j][k]->SetName( subDifBaseName.c_str() );
        mixCentVz[i][j][k] = (TH3D*) mixFiles[i]->Get( mixDifInitName.c_str() );
        mixCentVz[i][j][k]->SetName( mixDifBaseName.c_str() );
        mixSubCentVz[i][j][k] = (TH3D*) mixFiles[i]->Get( mixSubDifInitName.c_str() );
        mixSubCentVz[i][j][k]->SetName( mixSubDifBaseName.c_str() );
      }
  }
  
  // setup for 2d projections along pt axis
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > corrCentVzPt;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > subCentVzPt;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > mixCentVzPt;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > mixSubCentVzPt;
  corrCentVzPt.resize( nFiles );
  subCentVzPt.resize( nFiles );
  mixCentVzPt.resize( nFiles );
  mixSubCentVzPt.resize( nFiles );
  for ( int i = 0; i < nFiles; ++i ) {
    corrCentVzPt[i].resize( corrAnalysis::binsCentrality );
    subCentVzPt[i].resize( corrAnalysis::binsCentrality );
    mixCentVzPt[i].resize( corrAnalysis::binsCentrality );
    mixSubCentVzPt[i].resize( corrAnalysis::binsCentrality );
    
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      corrCentVzPt[i][j].resize( corrAnalysis::binsVz );
      subCentVzPt[i][j].resize( corrAnalysis::binsVz );
      mixCentVzPt[i][j].resize( corrAnalysis::binsVz );
      mixSubCentVzPt[i][j].resize( corrAnalysis::binsVz );
      
      for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
        corrCentVzPt[i][j][k].resize( nPtBins );
        subCentVzPt[i][j][k].resize( nPtBins );
        mixCentVzPt[i][j][k].resize( nPtBins );
        mixSubCentVzPt[i][j][k].resize( nPtBins );
      }
    }
  }
  
  // now get the pt projections
  for ( int i = 0; i < nFiles; ++i ) {
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      for ( int k = 0; k < corrAnalysis::binsVz; ++ k ) {
        for ( int l = 0; l < nPtBins; ++l ) {
          
          corrCentVz[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          subCentVz[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          mixCentVz[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          mixSubCentVz[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          
          corrCentVzPt[i][j][k][l] = (TH2D*) ((TH2D*) corrCentVz[i][j][k]->Project3D( "YX" ))->Clone();
          subCentVzPt[i][j][k][l] = (TH2D*) ((TH2D*) subCentVz[i][j][k]->Project3D("YX"))->Clone();
          mixCentVzPt[i][j][k][l] = (TH2D*) ((TH2D*) mixCentVz[i][j][k]->Project3D( "YX" ))->Clone();
          mixSubCentVzPt[i][j][k][l] = (TH2D*) ((TH2D*) mixSubCentVz[i][j][k]->Project3D("YX"))->Clone();
          
        }
      }
    }
  }
  
  // scale the mixing histograms
  for ( int i = 0; i < nFiles; ++i ) {
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
        for ( int l = 0; l < nPtBins; ++l ) {
          mixCentVzPt[i][j][k][l]->Scale( 1.0/mixCentVzPt[i][j][k][l]->GetMaximum() );
          mixSubCentVzPt[i][j][k][l]->Scale( 1.0/mixSubCentVzPt[i][j][k][l]->GetMaximum() );
        }
      }
    }
  }
  
  // TESTING PAST HERE
  // going back to using each bin independently not averaging over vz and cent
  // ( that is found in generate_output_avg.cxx )
  
  // make the container for the recombined histograms
  std::vector<std::vector<TH2D*> > recombinedCorr;
  std::vector<std::vector<TH2D*> > recombinedSub;
  std::vector<std::vector<TH2D*> > recombinedPre;
  std::vector<std::vector<TH2D*> > recombinedSubPre;
  recombinedCorr.resize( nFiles );
  recombinedSub.resize( nFiles );
  recombinedPre.resize( nFiles );
  recombinedSubPre.resize( nFiles );
  
  for (int i = 0; i < nFiles; ++ i ) {
    
    recombinedCorr[i].resize( nPtBins );
    recombinedSub[i].resize( nPtBins );
    recombinedPre[i].resize( nPtBins );
    recombinedSubPre[i].resize( nPtBins );
    
    for ( int l = 0; l < nPtBins; ++l ) {
      
      std::string corrName = analysisNames[i] + " " + ptBinString[l];
      std::string subName = analysisNames[i] + "_sub " + ptBinString[l];
      std::string preName = "pre_" + analysisNames[i] + " " + ptBinString[l];
      std::string subPreName = "pre_" + analysisNames[i] + "_sub " + ptBinString[l];
      
      recombinedCorr[i][l] = new TH2D( corrName.c_str(), corrName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      recombinedSub[i][l] = new TH2D( subName.c_str(), subName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      recombinedPre[i][l] = new TH2D( preName.c_str(), preName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      recombinedSubPre[i][l] = new TH2D( subPreName.c_str(), subPreName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      if ( l <= 2 ) {
      
        for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
          for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
            if ( mixCentVzPt[i][j][k][l]->GetEntries() != 0 && corrCentVzPt[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedPre[i][l]->Add( corrCentVzPt[i][j][k][l] );
              
              corrCentVzPt[i][j][k][l]->Divide( mixCentVzPt[i][j][k][l] );
            
              recombinedCorr[i][l]->Add( corrCentVzPt[i][j][k][l] );
            }
            if ( subCentVzPt[i][j][k][l]->GetEntries() != 0 && mixSubCentVzPt[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedSubPre[i][l]->Add( subCentVzPt[i][j][k][l] );
              
              subCentVzPt[i][j][k][l]->Divide( mixSubCentVzPt[i][j][k][l] );
              
              recombinedSub[i][l]->Add( subCentVzPt[i][j][k][l] );
              
            }
          }
        }
      }
      
      else {
        for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
          for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
            if ( mixCentVzPt[i][j][k][2]->GetEntries() != 0 && corrCentVzPt[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedPre[i][l]->Add( corrCentVzPt[i][j][k][l] );
              
              corrCentVzPt[i][j][k][l]->Divide( mixCentVzPt[i][j][k][2] );
              
              recombinedCorr[i][l]->Add( corrCentVzPt[i][j][k][l] );
            }
            if ( subCentVzPt[i][j][k][l]->GetEntries() != 0 && mixSubCentVzPt[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedSubPre[i][l]->Add( subCentVzPt[i][j][k][l] );
              
              subCentVzPt[i][j][k][l]->Divide( mixSubCentVzPt[i][j][k][2] );
              
              recombinedSub[i][l]->Add( subCentVzPt[i][j][k][l] );
              
            }
          }
        }
      }
    }
  }
  
  // get the reduced eta and phi ranges for projections
  double etaMax = 1.4;
  double etaMin = -1.4;
  double etaNearMin = -0.7;
  double etaNearMax = 0.7;
  double phiMin = -corrAnalysis::pi/2.0;
  double phiMaxClose = corrAnalysis::pi/2.0;
  double phiMaxFar = 3.0*corrAnalysis::pi/2.0;
  double phiMax = 3.0*corrAnalysis::pi/2.0;
  
  
  // going to get the 1D projections
  std::vector<std::vector<TH1D*> > dPhiLead;
  std::vector<std::vector<TH1D*> > dPhiLeadNear;
  std::vector<std::vector<TH1D*> > dPhiLeadFar;
  std::vector<std::vector<TH1D*> > dEtaLead;
  std::vector<std::vector<TH1D*> > dPhiSub;
  std::vector<std::vector<TH1D*> > dPhiSubNear;
  std::vector<std::vector<TH1D*> > dPhiSubFar;
  std::vector<std::vector<TH1D*> > dEtaSub;
  
  dPhiLead.resize( nFiles );
  dPhiLeadNear.resize( nFiles );
  dPhiLeadFar.resize( nFiles );
  dEtaLead.resize( nFiles );
  dPhiSub.resize( nFiles );
  dPhiSubNear.resize( nFiles );
  dPhiSubFar.resize( nFiles );
  dEtaSub.resize( nFiles );
  
  for ( int i = 0; i < nFiles; ++i ) {
    dPhiLead[i].resize( nPtBins );
    dPhiLeadNear[i].resize( nPtBins );
    dPhiLeadFar[i].resize( nPtBins );
    dEtaLead[i].resize( nPtBins );
    dPhiSub[i].resize( nPtBins );
    dPhiSubNear[i].resize( nPtBins );
    dPhiSubFar[i].resize( nPtBins );
    dEtaSub[i].resize( nPtBins );
    
    for ( int j = 0; j < nPtBins; ++j ) {
      // first restrict the eta range
      recombinedCorr[i][j]->GetXaxis()->SetRangeUser( etaMin, etaMax  );
      recombinedSub[i][j]->GetXaxis()->SetRangeUser( etaMin, etaMax );
      
      dPhiLead[i][j] = (TH1D*) ((TH1D*) recombinedCorr[i][j]->ProjectionY())->Clone();
      dPhiSub[i][j] = (TH1D*) ((TH1D*) recombinedSub[i][j]->ProjectionY())->Clone();
      
      recombinedCorr[i][j]->GetYaxis()->SetRangeUser( phiMin, phiMaxClose );
      recombinedSub[i][j]->GetYaxis()->SetRangeUser( phiMin, phiMaxClose );
      
      dEtaLead[i][j] = (TH1D*) ((TH1D*) recombinedCorr[i][j]->ProjectionX())->Clone();
      dEtaSub[i][j] = (TH1D*) ((TH1D*) recombinedSub[i][j]->ProjectionX())->Clone();
      
      recombinedCorr[i][j]->GetYaxis()->SetRangeUser( phiMin, phiMaxFar );
      recombinedSub[i][j]->GetYaxis()->SetRangeUser( phiMin, phiMaxFar );
      
      // now get dphi in "near" and "far" eta ranges
      recombinedCorr[i][j]->GetXaxis()->SetRangeUser( etaNearMin, etaNearMax  );
      recombinedSub[i][j]->GetXaxis()->SetRangeUser( etaNearMin, etaNearMax );
      
      dPhiLeadNear[i][j] = (TH1D*) ((TH1D*) recombinedCorr[i][j]->ProjectionY())->Clone();
      dPhiSubNear[i][j] = (TH1D*) ((TH1D*) recombinedSub[i][j]->ProjectionY())->Clone();
      
      recombinedCorr[i][j]->GetXaxis()->SetRangeUser( etaNearMax, etaMax  );
      recombinedSub[i][j]->GetXaxis()->SetRangeUser( etaNearMax, etaMax );
      
      dPhiLeadFar[i][j] = (TH1D*) ((TH1D*) recombinedCorr[i][j]->ProjectionY())->Clone();
      dPhiSubFar[i][j] = (TH1D*) ((TH1D*) recombinedSub[i][j]->ProjectionY())->Clone();
      
      recombinedCorr[i][j]->GetXaxis()->SetRangeUser( etaMin, etaNearMin  );
      recombinedSub[i][j]->GetXaxis()->SetRangeUser( etaMin, etaNearMin );
      
      dPhiLeadFar[i][j]->Add( recombinedCorr[i][j]->ProjectionY() );
      dPhiSubFar[i][j]->Add( recombinedSub[i][j]->ProjectionY() );
      
      // Now do the subtraction
      dPhiLeadNear[i][j]->Add( dPhiLeadFar[i][j], -1 );
      dPhiSubNear[i][j]->Add( dPhiSubFar[i][j], -1 );
      
    }
  }
  
  // Now to do some fitting and subtract the background
  // define the fits
  // ---------------
  std::string phiForm = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)+[4]*exp(-0.5*((x-[5])/[6])**2)";
  std::string etaForm = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)";
  
  // do a first, temporary fit to remove background
  for ( int i = 0; i < nFiles; ++i ) {
    for ( int j = 0; j < nPtBins; ++j ) {
      std::string dPhiLeadName = "fit_lead_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubName = "fit_sub_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiLeadNameDif = "fit_lead_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubNameDif = "fit_sub_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaLeadName = "fit_lead_eta_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaSubName = "fit_sub_eta_" + patch::to_string(i) + patch::to_string(j);
      
      TF1* leadPhiInitFit = new TF1( dPhiLeadName.c_str(), phiForm.c_str(), phiMin, phiMax );
      leadPhiInitFit->FixParameter( 2, 0 );
      leadPhiInitFit->SetParameter( 5, corrAnalysis::pi );
      leadPhiInitFit->SetParameter( 3, 0.1 );
      leadPhiInitFit->SetParameter( 6, 0.1 );
      
      TF1* leadPhiDifInitFit = new TF1( dPhiLeadNameDif.c_str(), phiForm.c_str(), phiMin, phiMax );
      leadPhiDifInitFit->FixParameter( 2, 0 );
      leadPhiDifInitFit->SetParameter( 5, corrAnalysis::pi );
      leadPhiDifInitFit->SetParameter( 3, 0.1 );
      leadPhiDifInitFit->SetParameter( 6, 0.1 );
      
      TF1* subPhiInitFit = new TF1( dPhiSubName.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiInitFit->FixParameter( 2, 0 );
      subPhiInitFit->SetParameter( 5, corrAnalysis::pi );
      subPhiInitFit->SetParameter( 3, 0.1 );
      subPhiInitFit->SetParameter( 6, 0.1 );
      
      TF1* subPhiDifInitFit = new TF1( dPhiSubNameDif.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiDifInitFit->FixParameter( 2, 0 );
      subPhiDifInitFit->SetParameter( 5, corrAnalysis::pi );
      subPhiDifInitFit->SetParameter( 3, 0.1 );
      subPhiDifInitFit->SetParameter( 6, 0.1 );
      
      TF1* leadEtaInitFit = new TF1( dEtaLeadName.c_str(), etaForm.c_str(), etaMin, etaMax );
      leadEtaInitFit->FixParameter( 2, 0 );
      leadEtaInitFit->SetParameter( 3, 0.1 );
      
      TF1* subEtaInitFit = new TF1( dEtaSubName.c_str(), etaForm.c_str(), etaMin, etaMax );
      subEtaInitFit->FixParameter( 2, 0 );
      subEtaInitFit->SetParameter( 3, 0.1 );
      
      dPhiLead[i][j]->Fit( dPhiLeadName.c_str(), "R" );
      dPhiSub[i][j]->Fit( dPhiSubName.c_str(), "R" );
      dPhiLeadNear[i][j]->Fit( dPhiLeadNameDif.c_str(), "R" );
      dPhiSubNear[i][j]->Fit( dPhiSubNameDif.c_str(), "R" );
      dEtaLead[i][j]->Fit( dEtaLeadName.c_str(), "R" );
      dEtaSub[i][j]->Fit( dEtaSubName.c_str(), "R" );
      
      // Now to subtract the constants
      TF1* subConst = new TF1( "subConst", "[0]", phiMin, phiMax);
      TF1* subConstEta = new TF1("subConstEta", "[0]", etaMin, etaMax);
      subConst->SetParameter( 0, leadPhiInitFit->GetParameter(0) );
      dPhiLead[i][j]->Add( subConst, -1 );
      subConst->SetParameter( 0, leadPhiDifInitFit->GetParameter(0));
      dPhiLeadNear[i][j]->Add( subConst, -1 );
      subConstEta->SetParameter( 0, leadEtaInitFit->GetParameter(0));
      dEtaLead[i][j]->Add( subConstEta, -1 );
      
      subConst->SetParameter( 0, subPhiInitFit->GetParameter(0) );
      dPhiSub[i][j]->Add( subConst, -1 );
      subConst->SetParameter( 0, subPhiDifInitFit->GetParameter(0));
      dPhiSubNear[i][j]->Add( subConst, -1 );
      subConstEta->SetParameter( 0, subEtaInitFit->GetParameter(0));
      dEtaSub[i][j]->Add( subConstEta, -1 );
      
      // Now reset the constants to zero...
      leadPhiInitFit->SetParameter( 0, 0 );
      leadPhiDifInitFit->SetParameter( 0, 0 );
      leadEtaInitFit->SetParameter( 0, 0 );
      subPhiInitFit->SetParameter( 0, 0 );
      subPhiDifInitFit->SetParameter( 0, 0 );
      subEtaInitFit->SetParameter( 0, 0 );
      
      
    }
  }
  
  TCanvas c1;
  dPhiLead[0][0]->Draw();
  c1.SaveAs("test.pdf");
  
  
  return 0;
}

  // find the bins by looping over the axes
  
  //  // first look for eta
  //  int etaMinBin, etaMaxBin, phiMinCloseBin, phiMinFarBin, phiMaxCloseBin, phiMaxFarBin;
  //
  //  for ( int i = 1; i <= recombinedCorr[0][0]->GetXaxis()->GetNbins(); ++i ) {
  //    if ( recombinedCorr[0][0]->GetXaxis()->GetBinLowEdge(i) >= etaMin && recombinedCorr[0][0]->GetXaxis()->GetBinLowEdge(i-1) < etaMin  )
  //      etaMinBin = i;
  //    if ( recombinedCorr[0][0]->GetXaxis()->GetBinUpEdge(i) > etaMax && recombinedCorr[0][0]->GetXaxis()->GetBinUpEdge(i-1) <= etaMax  )
  //      etaMaxBin = i;
  //  }
  //
  //  // we can manually set the low and high for phi
  //  phiMinCloseBin = 1;
  //  phiMaxFarBin = recombinedCorr[0][0]->GetXaxis()->GetNbins();
  //
  //  for ( int i = 1; i <= recombinedCorr[0][0]->GetYaxis()->GetNbins(); ++i ) {
  //    if ( recombinedCorr[0][0]->GetYaxis()->GetBinLowEdge(i) >= phiMaxClose && recombinedCorr[0][0]->GetYaxis()->GetBinLowEdge(i-1) < phiMaxClose  ) {
  //      phiMinFarBin = i;
  //      phiMaxCloseBin = i-1;
  //
  //    }
  //
  //  }
  
  //  // test output
  //  TCanvas c1;
  //  for ( int i = 0; i < nFiles; ++i ) {
  //    for ( int j = 0; j < nPtBins; ++j ) {
  //       std::string preCorrNameOut = "tmp/pre_" + analysisNames[i]; preCorrNameOut += ptBinString[j]; preCorrNameOut += ".pdf";
  //      std::string corrNameOut = "tmp/" + analysisNames[i]; corrNameOut += ptBinString[j]; corrNameOut += ".pdf";
  //      std::string mixNameOut = "tmp/" + analysisNames[i]; mixNameOut += ptBinString[j]; mixNameOut += " Mix.pdf";
  //      std::string projYNameOut = "tmp/" + analysisNames[i]; projYNameOut += ptBinString[j]; projYNameOut += "projectY.pdf";
  //      std::string projXNameOut = "tmp/" + analysisNames[i]; projXNameOut += ptBinString[j]; projXNameOut += "projectX.pdf";
  //      std::string preProjYNameOut = "tmp/pre_" + analysisNames[i]; preProjYNameOut += ptBinString[j]; preProjYNameOut += "projectY.pdf";
  //      std::string preProjXNameOut = "tmp/pre_" + analysisNames[i]; preProjXNameOut += ptBinString[j]; preProjXNameOut += "projectX.pdf";
  //
  //      recombinedPre[i][j]->Draw( "surf1" );
  //      c1.SaveAs( preCorrNameOut.c_str() );
  //      recombinedPre[i][j]->ProjectionY()->Draw();
  //      c1.SaveAs( preProjYNameOut.c_str() );
  //      recombinedPre[i][j]->ProjectionX()->Draw();
  //      c1.SaveAs( preProjXNameOut.c_str() );
  //
  //
  //      recombinedCorr[i][j]->Draw( "surf1" );
  //      c1.SaveAs(corrNameOut.c_str() );
  //      recombinedCorr[i][j]->ProjectionY()->Draw();
  //      c1.SaveAs( projYNameOut.c_str() );
  //      recombinedCorr[i][j]->ProjectionX()->Draw();
  //      c1.SaveAs( projXNameOut.c_str() );
  //
  //
  //      std::string postProjYNameOut = "tmp/post_" + analysisNames[i]; postProjYNameOut += ptBinString[j]; postProjYNameOut += "projectY.pdf";
  //
  //      recombinedCorr[i][j]->GetXaxis()->SetRangeUser( etaMin, etaMax );
  //
  //      recombinedCorr[i][j]->ProjectionY()->Draw();
  //      c1.SaveAs( postProjYNameOut.c_str() );
  //
  //      std::string postProjXNameOutNear = "tmp/post_" + analysisNames[i]; postProjXNameOutNear += ptBinString[j]; postProjXNameOutNear += "projectXNear.pdf";
  //      std::string postProjXNameOutFar = "tmp/post_" + analysisNames[i]; postProjXNameOutFar += ptBinString[j]; postProjXNameOutFar += "projectXFar.pdf";
  //
  //      recombinedCorr[i][j]->GetYaxis()->SetRange( phiMinCloseBin, phiMaxCloseBin );
  //      recombinedCorr[i][j]->ProjectionX()->Draw();
  //      c1.SaveAs( postProjXNameOutNear.c_str() );
  //
  //      recombinedCorr[i][j]->GetYaxis()->SetRange( phiMinFarBin, phiMaxFarBin );
  //      recombinedCorr[i][j]->ProjectionX()->Draw();
  //      c1.SaveAs( postProjXNameOutFar.c_str() );
  //    }
  //  }
  //  
