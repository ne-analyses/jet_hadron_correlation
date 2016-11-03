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
  
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(true);
  
  // Pt bins
  const int nPtBins = 5;
  const int startPtBin = 1;
  double ptBinLo[nPtBins] = { 3, 5, 9, 13, 17 };
  double ptBinHi[nPtBins] = { 4, 8, 12, 16, 24 };
  std::string ptBinString[nPtBins] = { "0.5 < p_{T} < 1.0", "1.0 < p_{T} < 2.0", "2.0 < p_{T} < 3.0", "3.0 < p_{T} < 4.0", "4.0 < p_{T} < 6.0" };
  double ptBinWidth[nPtBins];
  for ( int i = 0; i < nPtBins; ++i ) {
    ptBinWidth[i] = ( ptBinHi[i] - ptBinLo[i] + 1 ) * 0.25;
  }
  
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
  
  // Load in the histograms and get the pt spectra
  TH2D* nEvents[ nFiles ];
  TH1D* hVz[ nFiles ];
  TH3D* corrHist[ nFiles ];
  TH3D* mixHist[ nFiles ];
  std::vector<std::vector<std::vector<TH3D*> > > corrCentVz;
  std::vector<std::vector<std::vector<TH3D*> > > subCentVz;
  std::vector<std::vector<std::vector<TH3D*> > > mixCentVz;
  std::vector<std::vector<std::vector<TH3D*> > > mixSubCentVz;
  std::vector<TH1D*> recombinedPtLead;
  std::vector<TH1D*> recombinedPtSub;
  corrCentVz.resize( nFiles );
  subCentVz.resize( nFiles );
  mixCentVz.resize( nFiles );
  mixSubCentVz.resize( nFiles );
  recombinedPtLead.resize( nFiles );
  recombinedPtSub.resize( nFiles );
  for ( int i = 0; i < nFiles; ++i ) {
    corrCentVz[i].resize( corrAnalysis::binsCentrality );
    subCentVz[i].resize( corrAnalysis::binsCentrality );
    mixCentVz[i].resize( corrAnalysis::binsCentrality );
    mixSubCentVz[i].resize( corrAnalysis::binsCentrality );
    
    std::string ptLeadName = analysisNames[i] + "_pt_lead";
    std::string ptSubName = analysisNames[i] + "_pt_sub";
    
    recombinedPtLead[i] = new TH1D( ptLeadName.c_str(), "p_{T} Spectrum Trigger Jet", corrAnalysis::binsPt, corrAnalysis::ptLowEdge, corrAnalysis::ptHighEdge );
    recombinedPtSub[i] = new TH1D( ptSubName.c_str(), "p_{T} Spectrum Recoil Jet", corrAnalysis::binsPt, corrAnalysis::ptLowEdge, corrAnalysis::ptHighEdge );
    
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      corrCentVz[i][j].resize( corrAnalysis::binsVz );
      subCentVz[i][j].resize( corrAnalysis::binsVz );
      mixCentVz[i][j].resize( corrAnalysis::binsVz );
      mixSubCentVz[i][j].resize( corrAnalysis::binsVz );
    }
  }
  
  for ( int i = 0; i < nFiles; ++i ) {
    std::cout<<"loading histograms - file: "<< i <<std::endl;
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
        
        // now we can get the pt spectrum as well
        recombinedPtLead[i]->Add( (TH1D*) corrCentVz[i][j][k]->Project3D("Z") );
        recombinedPtSub[i]->Add( (TH1D*) subCentVz[i][j][k]->Project3D("Z") );
      }
  }
  
  std::cout<<"loaded all histograms"<<std::endl;
  
  // Printing out eta and phi bins
  std::cout<<"ETA BINS: "<<mixSubCentVz[0][0][0]->GetXaxis()->GetNbins()<<std::endl;
  for ( int i = 1; i <= mixSubCentVz[0][0][0]->GetXaxis()->GetNbins(); ++i ) {
    std::cout<<"bin: "<<i<<" low: "<< mixSubCentVz[0][0][0]->GetXaxis()->GetBinLowEdge(i)<<" upper: "<<mixSubCentVz[0][0][0]->GetXaxis()->GetBinUpEdge(i)<<std::endl;
  }
  
  std::cout<<"PHI BINS: "<<mixSubCentVz[0][0][0]->GetYaxis()->GetNbins()<<std::endl;
  for ( int i = 1; i <= mixSubCentVz[0][0][0]->GetYaxis()->GetNbins(); ++i ) {
    std::cout<<"bin: "<<i<<" low: "<< mixSubCentVz[0][0][0]->GetYaxis()->GetBinLowEdge(i)<<" upper: "<<mixSubCentVz[0][0][0]->GetYaxis()->GetBinUpEdge(i)<<std::endl;
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
  double etaMax = 1.2;
  double etaMin = -1.2;
  double etaNearMin = etaMin/2.0;
  double etaNearMax = etaMax/2.0;
  int etaMinBin = 5;
  int etaMaxBin = etaMinBin + 3;
  int etaNearMinBin = etaMaxBin + 1;
  int etaNearMaxBin = etaNearMinBin + 7;
  int etaFarMinBin = etaNearMaxBin + 1;
  int etaFarMaxBin = etaFarMinBin + 3;
  double phiMin = -corrAnalysis::pi/2.0;
  double phiMinClose = -0.6;
  double phiMaxClose = 0.6;
  double phiMaxFar = 3.0*corrAnalysis::pi/2.0;
  double phiMax = 3.0*corrAnalysis::pi/2.0;
  double phiDifMax = corrAnalysis::pi/2.0;
  
  
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
      
      // save the 2D histograms
      std::string leadOutName = "tmp/lead2d_" + analysisNames[i] +"_pt_" +patch::to_string(j) + ".pdf";
      std::string subOutName = "tmp/sub2d_" + analysisNames[i] +"_pt_"+ patch::to_string(j) + ".pdf";
      
      TCanvas c1;
      recombinedCorr[i][j]->Draw("surf1");
      c1.SaveAs( leadOutName.c_str() );
      recombinedSub[i][j]->Draw("surf1");
      c1.SaveAs( subOutName.c_str() );
      
      
      dPhiLead[i][j] = (TH1D*) ((TH1D*) recombinedCorr[i][j]->ProjectionY())->Clone();
      dPhiSub[i][j] = (TH1D*) ((TH1D*) recombinedSub[i][j]->ProjectionY())->Clone();
      
      recombinedCorr[i][j]->GetYaxis()->SetRangeUser( phiMinClose, phiMaxClose );
      recombinedSub[i][j]->GetYaxis()->SetRangeUser( phiMinClose, phiMaxClose );
      
      dEtaLead[i][j] = (TH1D*) ((TH1D*) recombinedCorr[i][j]->ProjectionX())->Clone();
      dEtaSub[i][j] = (TH1D*) ((TH1D*) recombinedSub[i][j]->ProjectionX())->Clone();
      
      recombinedCorr[i][j]->GetYaxis()->SetRangeUser( phiMin, phiMaxFar );
      recombinedSub[i][j]->GetYaxis()->SetRangeUser( phiMin, phiMaxFar );
      
      // now get dphi in "near" and "far" eta ranges
      recombinedCorr[i][j]->GetXaxis()->SetRange( etaNearMinBin, etaNearMaxBin  );
      recombinedSub[i][j]->GetXaxis()->SetRange( etaNearMinBin, etaNearMaxBin );
      
      dPhiLeadNear[i][j] = (TH1D*) ((TH1D*) recombinedCorr[i][j]->ProjectionY())->Clone();
      dPhiSubNear[i][j] = (TH1D*) ((TH1D*) recombinedSub[i][j]->ProjectionY())->Clone();
      
      recombinedCorr[i][j]->GetXaxis()->SetRange( etaFarMinBin, etaFarMaxBin  );
      recombinedSub[i][j]->GetXaxis()->SetRange( etaFarMinBin, etaFarMaxBin );
      
      dPhiLeadFar[i][j] = (TH1D*) ((TH1D*) recombinedCorr[i][j]->ProjectionY())->Clone();
      dPhiSubFar[i][j] = (TH1D*) ((TH1D*) recombinedSub[i][j]->ProjectionY())->Clone();
      
      recombinedCorr[i][j]->GetXaxis()->SetRange( etaMinBin, etaMaxBin  );
      recombinedSub[i][j]->GetXaxis()->SetRange( etaMinBin, etaMaxBin );
      
      dPhiLeadFar[i][j]->Add( recombinedCorr[i][j]->ProjectionY() );
      dPhiSubFar[i][j]->Add( recombinedSub[i][j]->ProjectionY() );
    }
  }
  
  // now  first overlay and output,
  // then subtract near from far eta regions
  for ( int i = 0; i < nFiles; ++i )
    for ( int j = 0; j < nPtBins; ++j ) {
      std::string leadPhiName = "tmp/lead_phi_near_far_"+analysisNames[i]+"_pt_"+patch::to_string(j)+".pdf";
      std::string subPhiName = "tmp/sub_phi_near_far_"+analysisNames[i]+"_pt_"+patch::to_string(j)+".pdf";
      
      TCanvas c1;
      dPhiLeadNear[i][j]->SetLineColor(kBlack);
      dPhiLeadNear[i][j]->SetMarkerStyle(29);
      dPhiLeadNear[i][j]->SetMarkerSize(3);
      dPhiLeadNear[i][j]->SetMarkerColor(kBlack);
      dPhiLeadNear[i][j]->Draw();
      dPhiLeadFar[i][j]->SetLineColor(kRed);
      dPhiLeadFar[i][j]->SetMarkerStyle(29);
      dPhiLeadFar[i][j]->SetMarkerSize(3);
      dPhiLeadFar[i][j]->SetMarkerColor(kRed);
      dPhiLeadFar[i][j]->Draw("SAME");
      c1.SaveAs( leadPhiName.c_str() );
      dPhiSubNear[i][j]->SetLineColor(kBlack);
      dPhiSubNear[i][j]->SetMarkerStyle(29);
      dPhiSubNear[i][j]->SetMarkerSize(3);
      dPhiSubNear[i][j]->SetMarkerColor(kBlack);
      dPhiSubNear[i][j]->Draw();
      dPhiSubFar[i][j]->SetLineColor(kRed);
      dPhiSubFar[i][j]->SetMarkerStyle(29);
      dPhiSubFar[i][j]->SetMarkerSize(3);
      dPhiSubFar[i][j]->SetMarkerColor(kRed);
      dPhiSubFar[i][j]->Draw("SAME");
      c1.SaveAs( subPhiName.c_str() );
      
      // Now do the subtraction
      dPhiLeadNear[i][j]->Add( dPhiLeadFar[i][j], -1 );
      dPhiSubNear[i][j]->Add( dPhiSubFar[i][j], -1 );
    }

  
  // Now to do some fitting and subtract the background
  // define the fits
  // ---------------
  std::string phiForm = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)+[4]*exp(-0.5*((x-[5])/[6])**2)";
  std::string etaForm = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)";
  std::string phiDifForm = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)";
  
  // do a first, temporary fit to remove background
  for ( int i = 0; i < nFiles; ++i ) {
    for ( int j = 0; j < nPtBins; ++j ) {
      std::string dPhiLeadName = "tmp_fit_lead_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubName = "tmp_fit_sub_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiLeadNameDif = "tmp_fit_lead_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubNameDif = "tmp_fit_sub_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaLeadName = "tmp_fit_lead_eta_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaSubName = "tmp_fit_sub_eta_" + patch::to_string(i) + patch::to_string(j);
      
      TF1* leadPhiInitFit = new TF1( dPhiLeadName.c_str(), phiForm.c_str(), phiMin, phiDifMax );
      leadPhiInitFit->FixParameter( 2, 0 );
      leadPhiInitFit->FixParameter( 5, corrAnalysis::pi );
      leadPhiInitFit->SetParameter( 3, 0.2 );
      leadPhiInitFit->SetParameter( 6, 0.2 );
      
      TF1* leadPhiDifInitFit = new TF1( dPhiLeadNameDif.c_str(), phiDifForm.c_str(), phiMin, phiMax );
      leadPhiDifInitFit->FixParameter( 2, 0 );
      leadPhiDifInitFit->SetParameter( 3, 0.2 );
      
      TF1* subPhiInitFit = new TF1( dPhiSubName.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiInitFit->FixParameter( 2, 0 );
      subPhiInitFit->FixParameter( 5, corrAnalysis::pi );
      subPhiInitFit->SetParameter( 3, 0.2 );
      subPhiInitFit->SetParameter( 6, 0.2 );
      
      TF1* subPhiDifInitFit = new TF1( dPhiSubNameDif.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      subPhiDifInitFit->FixParameter( 2, 0 );
      subPhiDifInitFit->SetParameter( 3, 0.2 );
      
      TF1* leadEtaInitFit = new TF1( dEtaLeadName.c_str(), etaForm.c_str(), etaMin, etaMax );
      leadEtaInitFit->FixParameter( 2, 0 );
      leadEtaInitFit->SetParameter( 3, 0.2 );
      
      TF1* subEtaInitFit = new TF1( dEtaSubName.c_str(), etaForm.c_str(), etaMin, etaMax );
      subEtaInitFit->FixParameter( 2, 0 );
      subEtaInitFit->SetParameter( 3, 0.2 );
      
      dPhiLead[i][j]->Fit( dPhiLeadName.c_str(), "RM" );
      dPhiSub[i][j]->Fit( dPhiSubName.c_str(), "RM" );
      dPhiLeadNear[i][j]->Fit( dPhiLeadNameDif.c_str(), "RM" );
      dPhiSubNear[i][j]->Fit( dPhiSubNameDif.c_str(), "RM" );
      dEtaLead[i][j]->Fit( dEtaLeadName.c_str(), "RM" );
      dEtaSub[i][j]->Fit( dEtaSubName.c_str(), "RM" );
      
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
      
      // now scale the histograms
      dPhiLead[i][j]->Scale( 1.0 / dPhiLead[i][j]->GetXaxis()->GetBinWidth(1) );
      dPhiLead[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
      dPhiLeadNear[i][j]->Scale( 1.0 / dPhiLeadNear[i][j]->GetXaxis()->GetBinWidth(1) );
      dPhiLeadNear[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
      dEtaLead[i][j]->Scale( 1.0 / dEtaLead[i][j]->GetXaxis()->GetBinWidth(1) );
      dEtaLead[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
      dPhiSub[i][j]->Scale( 1.0 / dPhiSub[i][j]->GetXaxis()->GetBinWidth(1) );
      dPhiSub[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
      dPhiSubNear[i][j]->Scale( 1.0 / dPhiSubNear[i][j]->GetXaxis()->GetBinWidth(1) );
      dPhiSubNear[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
      dEtaSub[i][j]->Scale( 1.0 / dEtaSub[i][j]->GetXaxis()->GetBinWidth(1) );
      dEtaSub[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
    }
  }
  
  // final fitting
  std::vector<std::vector<TF1*> > leadPhiFit;
  leadPhiFit.resize( nFiles );
  std::vector<std::vector<TF1*> > leadPhiDifFit;
  leadPhiDifFit.resize( nFiles );
  std::vector<std::vector<TF1*> > leadEtaFit;
  leadEtaFit.resize( nFiles );
  std::vector<std::vector<TF1*> > subPhiFit;
  subPhiFit.resize( nFiles );
  std::vector<std::vector<TF1*> > subPhiDifFit;
  subPhiDifFit.resize( nFiles );
  std::vector<std::vector<TF1*> > subEtaFit;
  subEtaFit.resize( nFiles );
  
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiFit[i].resize( nPtBins );
    leadPhiDifFit[i].resize( nPtBins );
    leadEtaFit[i].resize( nPtBins );
    subPhiFit[i].resize( nPtBins );
    subPhiDifFit[i].resize( nPtBins );
    subEtaFit[i].resize( nPtBins );
    
    for ( int j = 0; j < nPtBins; ++j ) {
      std::string dPhiLeadName = "fit_lead_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubName = "fit_sub_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiLeadNameDif = "fit_lead_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubNameDif = "fit_sub_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaLeadName = "fit_lead_eta_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaSubName = "fit_sub_eta_" + patch::to_string(i) + patch::to_string(j);
    
      leadPhiFit[i][j] = new TF1( dPhiLeadName.c_str(), phiForm.c_str(), phiMin, phiMax );
      leadPhiFit[i][j]->FixParameter( 2, 0 );
      leadPhiFit[i][j]->FixParameter( 5, corrAnalysis::pi );
      leadPhiFit[i][j]->SetParameter( 3, 0.2 );
      leadPhiFit[i][j]->SetParameter( 6, 0.2 );
      leadPhiFit[i][j]->SetLineColor( i + 1 );
      
      leadPhiDifFit[i][j] = new TF1( dPhiLeadNameDif.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      leadPhiDifFit[i][j]->FixParameter( 2, 0 );
      leadPhiDifFit[i][j]->SetParameter( 3, 0.2 );
      leadPhiDifFit[i][j]->SetLineColor( i + 1 );
      
      subPhiFit[i][j] = new TF1( dPhiSubName.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiFit[i][j]->FixParameter( 2, 0 );
      subPhiFit[i][j]->FixParameter( 5, corrAnalysis::pi );
      subPhiFit[i][j]->SetParameter( 3, 0.2 );
      subPhiFit[i][j]->SetParameter( 6, 0.2 );
      subPhiFit[i][j]->SetLineColor( i + 1 );
      
      subPhiDifFit[i][j] = new TF1( dPhiSubNameDif.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      subPhiDifFit[i][j]->FixParameter( 2, 0 );
      subPhiDifFit[i][j]->SetParameter( 3, 0.2 );
      subPhiDifFit[i][j]->SetLineColor( i + 1 );
      
      leadEtaFit[i][j] = new TF1( dEtaLeadName.c_str(), etaForm.c_str(), etaMin, etaMax );
      leadEtaFit[i][j]->FixParameter( 2, 0 );
      leadEtaFit[i][j]->SetParameter( 3, 0.2 );
      leadEtaFit[i][j]->SetLineColor( i + 1 );
      
      subEtaFit[i][j] = new TF1( dEtaSubName.c_str(), etaForm.c_str(), etaMin, etaMax );
      subEtaFit[i][j]->FixParameter( 2, 0 );
      subEtaFit[i][j]->SetParameter( 3, 0.2 );
      subEtaFit[i][j]->SetLineColor( i + 1 );

      // Now set same colors and fit
      dPhiLead[i][j]->SetLineColor( i + 1 );
      dPhiLead[i][j]->Fit( dPhiLeadName.c_str(), "RM" );
      dPhiSub[i][j]->SetLineColor( i + 1 );
      dPhiSub[i][j]->Fit( dPhiSubName.c_str(), "RM" );
      dPhiLeadNear[i][j]->SetLineColor( i + 1 );
      dPhiLeadNear[i][j]->Fit( dPhiLeadNameDif.c_str(), "RM" );
      dPhiSubNear[i][j]->SetLineColor( i + 1 );
      dPhiSubNear[i][j]->Fit( dPhiSubNameDif.c_str(), "RM" );
      dEtaLead[i][j]->SetLineColor( i + 1 );
      dEtaLead[i][j]->Fit( dEtaLeadName.c_str(), "RM" );
      dEtaSub[i][j]->SetLineColor( i + 1 );
      dEtaSub[i][j]->Fit( dEtaSubName.c_str(), "RM" );
      
    }
  }
  
  
  // Now start making output
  std::string outBase = "tmp/";
  std::string leadPhiOutBase = outBase + "leadphi_pt";
  std::string leadPhiDifOutBase = outBase + "leadphidif_pt";
  std::string leadEtaOutBase = outBase + "leadeta_pt";
  std::string subPhiOutBase = outBase + "subphi_pt";
  std::string subPhiDifOutBase = outBase + "subphidif_pt";
  std::string subEtaOutBase = outBase + "subeta_pt";
  std::string outExt = ".pdf";
  
  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string leadPhiOut = leadPhiOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      if ( j == 0 ) {
        std::string outTitle = "Trigger Jet #Delta#phi " + ptBinString[i];
        dPhiLead[j][i]->SetTitle( outTitle.c_str() );
        dPhiLead[j][i]->SetLineColor(j+1);
        dPhiLead[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiLead[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiLead[j][i]->SetMarkerStyle(29);
        dPhiLead[j][i]->SetMarkerSize(3);
        dPhiLead[j][i]->SetMarkerColor(j+1);
        dPhiLead[j][i]->Draw();
      }
      else {
        dPhiLead[j][i]->Draw("same");
      }
    }
    c1.SaveAs( leadPhiOut.c_str() );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string leadPhiDifOut = leadPhiDifOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      if ( j == 0 ) {
        std::string outTitle = "Trigger Jet #Delta#eta subtracted #Delta#phi " + ptBinString[i];
        dPhiLeadNear[j][i]->SetTitle( outTitle.c_str() );
        dPhiLeadNear[j][i]->GetXaxis()->SetRangeUser( -corrAnalysis::pi/2.0, corrAnalysis::pi/2.0 );
        dPhiLeadNear[j][i]->SetLineColor(j+1);
        dPhiLeadNear[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiLeadNear[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiLeadNear[j][i]->SetMarkerStyle(29);
        dPhiLeadNear[j][i]->SetMarkerSize(3);
        dPhiLeadNear[j][i]->SetMarkerColor(j+1);
        dPhiLeadNear[j][i]->Draw();
      }
      else {
        dPhiLeadNear[j][i]->Draw("same");
      }
    }
    c1.SaveAs( leadPhiDifOut.c_str() );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string leadEtaOut = leadEtaOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      if ( j == 0 ) {
        std::string outTitle = "Trigger Jet #Delta#eta " + ptBinString[i];
        dEtaLead[j][i]->SetTitle( outTitle.c_str() );
        dEtaLead[j][i]->SetLineColor(j+1);
        dEtaLead[j][i]->GetXaxis()->SetTitle("#Delta#eta");
        dEtaLead[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#eta");
        dEtaLead[j][i]->SetMarkerStyle(29);
        dEtaLead[j][i]->SetMarkerSize(3);
        dEtaLead[j][i]->SetMarkerColor(j+1);
        dEtaLead[j][i]->Draw();
      }
      else {
        dEtaLead[j][i]->Draw("same");
      }
    }
    c1.SaveAs( leadEtaOut.c_str() );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string subPhiOut = subPhiOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      if ( j == 0 ) {
        std::string outTitle = "Recoil Jet #Delta#phi " + ptBinString[i];
        dPhiSub[j][i]->SetTitle( outTitle.c_str() );
        dPhiSub[j][i]->SetLineColor(j+1);
        dPhiSub[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiSub[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiSub[j][i]->SetMarkerStyle(29);
        dPhiSub[j][i]->SetMarkerSize(3);
        dPhiSub[j][i]->SetMarkerColor(j+1);
        dPhiSub[j][i]->Draw();
      }
      else {
        dPhiSub[j][i]->Draw("same");
      }
    }
    c1.SaveAs( subPhiOut.c_str() );
  }

  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string subPhiDifOut = subPhiDifOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      if ( j == 0 ) {
        std::string outTitle = "Recoil Jet #Delta#eta subtracted #Delta#phi " + ptBinString[i];
        dPhiSubNear[j][i]->SetTitle( outTitle.c_str() );
        dPhiSubNear[j][i]->GetXaxis()->SetRangeUser( -corrAnalysis::pi/2.0, corrAnalysis::pi/2.0 );
        dPhiSubNear[j][i]->SetLineColor(j+1);
        dPhiSubNear[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiSubNear[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiSubNear[j][i]->SetMarkerStyle(29);
        dPhiSubNear[j][i]->SetMarkerSize(3);
        dPhiSubNear[j][i]->SetMarkerColor(j+1);
        dPhiSubNear[j][i]->Draw();
      }
      else {
        dPhiSubNear[j][i]->Draw("same");
      }
    }
    c1.SaveAs( subPhiDifOut.c_str() );
  }

  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string subEtaOut = subEtaOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      if ( j == 0 ) {
        std::string outTitle = "Recoil Jet #Delta#eta " + ptBinString[i];
        dEtaSub[j][i]->SetTitle( outTitle.c_str() );
        dEtaSub[j][i]->SetLineColor(j+1);
        dEtaSub[j][i]->GetXaxis()->SetTitle("#Delta#eta");
        dEtaSub[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#eta");
        dEtaSub[j][i]->SetMarkerStyle(29);
        dEtaSub[j][i]->SetMarkerSize(3);
        dEtaSub[j][i]->SetMarkerColor(j+1);
        dEtaSub[j][i]->Draw();
      }
      else {
        dEtaSub[j][i]->Draw("same");
      }
    }
    c1.SaveAs( subEtaOut.c_str() );
  }
  
  // now to get yields
  std::vector<std::vector<double> > leadPhiYield( nFiles );
  std::vector<std::vector<double> > leadPhiWidth( nFiles );
  std::vector<std::vector<double> > leadPhiError( nFiles );
  std::vector<std::vector<double> > leadPhiWidthError( nFiles );
  std::vector<std::vector<double> > leadPhiDifYield( nFiles );
  std::vector<std::vector<double> > leadPhiDifWidth( nFiles );
  std::vector<std::vector<double> > leadPhiDifError( nFiles );
  std::vector<std::vector<double> > leadPhiDifWidthError( nFiles );
  std::vector<std::vector<double> > leadEtaYield( nFiles );
  std::vector<std::vector<double> > leadEtaWidth( nFiles );
  std::vector<std::vector<double> > leadEtaError( nFiles );
  std::vector<std::vector<double> > leadEtaWidthError( nFiles );
  std::vector<std::vector<double> > subPhiYield( nFiles );
  std::vector<std::vector<double> > subPhiWidth( nFiles );
  std::vector<std::vector<double> > subPhiError( nFiles );
  std::vector<std::vector<double> > subPhiWidthError( nFiles );
  std::vector<std::vector<double> > subPhiDifYield( nFiles );
  std::vector<std::vector<double> > subPhiDifWidth( nFiles );
  std::vector<std::vector<double> > subPhiDifError( nFiles );
  std::vector<std::vector<double> > subPhiDifWidthError( nFiles );
  std::vector<std::vector<double> > subEtaYield( nFiles );
  std::vector<std::vector<double> > subEtaWidth( nFiles );
  std::vector<std::vector<double> > subEtaError( nFiles );
  std::vector<std::vector<double> > subEtaWidthError( nFiles );
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiYield[i].resize( nPtBins );
    leadPhiWidth[i].resize( nPtBins );
    leadPhiError[i].resize( nPtBins );
    leadPhiWidthError[i].resize( nPtBins );
    leadPhiDifYield[i].resize( nPtBins );
    leadPhiDifWidth[i].resize( nPtBins );
    leadPhiDifError[i].resize( nPtBins );
    leadPhiDifWidthError[i].resize( nPtBins );
    leadEtaYield[i].resize( nPtBins );
    leadEtaWidth[i].resize( nPtBins );
    leadEtaError[i].resize( nPtBins );
    leadEtaWidthError[i].resize( nPtBins );
    subPhiYield[i].resize( nPtBins );
    subPhiWidth[i].resize( nPtBins );
    subPhiError[i].resize( nPtBins );
    subPhiWidthError[i].resize( nPtBins );
    subPhiDifYield[i].resize( nPtBins );
    subPhiDifWidth[i].resize( nPtBins );
    subPhiDifError[i].resize( nPtBins );
    subPhiDifWidthError[i].resize( nPtBins );
    subEtaYield[i].resize( nPtBins );
    subEtaWidth[i].resize( nPtBins );
    subEtaError[i].resize( nPtBins );
    subEtaWidthError[i].resize( nPtBins );
    for ( int j = 0; j < nPtBins; ++j ) {
      leadPhiYield[i][j] = leadPhiFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(leadPhiFit[i][j]->GetParameter(3))/ptBinWidth[j];
      leadPhiError[i][j] = leadPhiFit[i][j]->GetParError(1);
      leadPhiWidth[i][j] = fabs(leadPhiFit[i][j]->GetParameter(3));
      leadPhiWidthError[i][j] = leadPhiFit[i][j]->GetParError(3);
      leadPhiDifYield[i][j] = leadPhiDifFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(leadPhiDifFit[i][j]->GetParameter(3))/ptBinWidth[j];
      leadPhiDifError[i][j] = leadPhiDifFit[i][j]->GetParError(1);
      leadPhiDifWidth[i][j] = fabs(leadPhiDifFit[i][j]->GetParameter(3));
      leadPhiDifWidthError[i][j] = leadPhiDifFit[i][j]->GetParError(3);
      leadEtaYield[i][j] = leadEtaFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(leadEtaFit[i][j]->GetParameter(3))/ptBinWidth[j];
      leadEtaError[i][j] = leadEtaFit[i][j]->GetParError(1);
      leadEtaWidth[i][j] = fabs(leadEtaFit[i][j]->GetParameter(3));
      leadEtaWidthError[i][j] = leadEtaFit[i][j]->GetParError(3);
      subPhiYield[i][j] = subPhiFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(subPhiFit[i][j]->GetParameter(3))/ptBinWidth[j];
      subPhiError[i][j] = subPhiFit[i][j]->GetParError(1);
      subPhiWidth[i][j] = fabs(subPhiFit[i][j]->GetParameter(3));
      subPhiWidthError[i][j] = subPhiFit[i][j]->GetParError(3);
      subPhiDifYield[i][j] = subPhiDifFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(subPhiDifFit[i][j]->GetParameter(3))/ptBinWidth[j];
      subPhiDifError[i][j] = subPhiDifFit[i][j]->GetParError(1);
      subPhiDifWidth[i][j] = fabs(subPhiDifFit[i][j]->GetParameter(3));
      subPhiDifWidthError[i][j] = subPhiDifFit[i][j]->GetParError(3);
      subEtaYield[i][j] = subEtaFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(subEtaFit[i][j]->GetParameter(3))/ptBinWidth[j];
      subEtaError[i][j] = subEtaFit[i][j]->GetParError(1);
      subEtaWidth[i][j] = fabs(subEtaFit[i][j]->GetParameter(3));
      subEtaWidthError[i][j] = subEtaFit[i][j]->GetParError(3);
    }
  }
  
  double ptBins[5] = { 0.75, 1.5, 2.5, 3.5, 5 };
  
  std::vector<TGraphErrors*> leadPhiGraph( nFiles );
  std::vector<TGraphErrors*> leadPhiWidthGraph( nFiles );
  std::vector<TGraphErrors*> leadPhiDifGraph( nFiles );
  std::vector<TGraphErrors*> leadPhiDifWidthGraph( nFiles );
  std::vector<TGraphErrors*> leadEtaGraph( nFiles );
  std::vector<TGraphErrors*> leadEtaWidthGraph( nFiles );
  std::vector<TGraphErrors*> subPhiGraph( nFiles );
  std::vector<TGraphErrors*> subPhiWidthGraph( nFiles );
  std::vector<TGraphErrors*> subPhiDifGraph( nFiles );
  std::vector<TGraphErrors*> subPhiDifWidthGraph( nFiles );
  std::vector<TGraphErrors*> subEtaGraph( nFiles );
  std::vector<TGraphErrors*> subEtaWidthGraph( nFiles );
  
  for ( int i = 0; i < nFiles; ++i ) {
    double leadPhiTmp[nPtBins-startPtBin];
    double leadPhiErr[nPtBins-startPtBin];
    double leadPhiWidthTmp[nPtBins-startPtBin];
    double leadPhiWidthErrTmp[nPtBins-startPtBin];
    double leadPhiDifTmp[nPtBins-startPtBin];
    double leadPhiDifErr[nPtBins-startPtBin];
    double leadPhiDifWidthTmp[nPtBins-startPtBin];
    double leadPhiDifWidthErrTmp[nPtBins-startPtBin];
    double leadEtaTmp[nPtBins-startPtBin];
    double leadEtaErr[nPtBins-startPtBin];
    double leadEtaWidthTmp[nPtBins-startPtBin];
    double leadEtaWidthErrTmp[nPtBins-startPtBin];
    double subPhiTmp[nPtBins-startPtBin];
    double subPhiErr[nPtBins-startPtBin];
    double subPhiWidthTmp[nPtBins-startPtBin];
    double subPhiWidthErrTmp[nPtBins-startPtBin];
    double subPhiDifTmp[nPtBins-startPtBin];
    double subPhiDifErr[nPtBins-startPtBin];
    double subPhiDifWidthTmp[nPtBins-startPtBin];
    double subPhiDifWidthErrTmp[nPtBins-startPtBin];
    double subEtaTmp[nPtBins-startPtBin];
    double subEtaErr[nPtBins-startPtBin];
    double subEtaWidthTmp[nPtBins-startPtBin];
    double subEtaWidthErrTmp[nPtBins-startPtBin];
    
    double errX[nPtBins - startPtBin];
    double tmpPtBin[nPtBins - startPtBin];
    for ( int j = 0; j < nPtBins-startPtBin; ++j ) {
      errX[j] = 0;
      tmpPtBin[j] = ptBins[j+startPtBin];
      leadPhiTmp[j] = leadPhiYield[i][j+startPtBin];
      leadPhiErr[j] = leadPhiError[i][j+startPtBin];
      leadPhiWidthTmp[j] = leadPhiWidth[i][j+startPtBin];
      leadPhiWidthErrTmp[j] = leadPhiWidthError[i][j+startPtBin];
      leadPhiDifTmp[j] = leadPhiDifYield[i][j+startPtBin];
      leadPhiDifErr[j] = leadPhiDifError[i][j+startPtBin];
      leadPhiDifWidthTmp[j] = leadPhiDifWidth[i][j+startPtBin];
      leadPhiDifWidthErrTmp[j] = leadPhiDifWidthError[i][j+startPtBin];
      leadEtaTmp[j] = leadEtaYield[i][j+startPtBin];
      leadEtaErr[j] = leadEtaError[i][j+startPtBin];
      leadEtaWidthTmp[j] = leadEtaWidth[i][j+startPtBin];
      leadEtaWidthErrTmp[j] = leadEtaWidthError[i][j+startPtBin];
      subPhiTmp[j] = subPhiYield[i][j+startPtBin];
      subPhiErr[j] = subPhiError[i][j+startPtBin];
      subPhiWidthTmp[j] = subPhiWidth[i][j+startPtBin];
      subPhiWidthErrTmp[j] = subPhiWidthError[i][j+startPtBin];
      subPhiDifTmp[j] = subPhiDifYield[i][j+startPtBin];
      subPhiDifErr[j] = subPhiDifError[i][j+startPtBin];
      subPhiDifWidthTmp[j] = subPhiDifWidth[i][j+startPtBin];
      subPhiDifWidthErrTmp[j] = subPhiDifWidthError[i][j+startPtBin];
      subEtaTmp[j] = subEtaYield[i][j+startPtBin];
      subEtaErr[j] = subEtaError[i][j+startPtBin];
      subEtaWidthTmp[j] = subEtaWidth[i][j+startPtBin];
      subEtaWidthErrTmp[j] = subEtaWidthError[i][j+startPtBin];
    }
    
    leadPhiGraph[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadPhiTmp, errX, leadPhiErr );
    leadPhiWidthGraph[i] = new TGraphErrors( nPtBins - startPtBin, tmpPtBin, leadPhiWidthTmp, errX, leadPhiWidthErrTmp );
    leadPhiDifGraph[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadPhiDifTmp, errX, leadPhiDifErr );
    leadPhiDifWidthGraph[i] = new TGraphErrors( nPtBins - startPtBin, tmpPtBin, leadPhiDifWidthTmp, errX, leadPhiDifWidthErrTmp );
    leadEtaGraph[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadEtaTmp, errX, leadEtaErr );
    leadEtaWidthGraph[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadEtaWidthTmp, errX, leadEtaWidthErrTmp );
    subPhiGraph[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subPhiTmp, errX, subPhiErr );
    subPhiWidthGraph[i] = new TGraphErrors(nPtBins-startPtBin, tmpPtBin, subPhiWidthTmp, errX, subPhiWidthErrTmp );
    subPhiDifGraph[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subPhiDifTmp, errX, subPhiDifErr );
    subPhiDifWidthGraph[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subPhiDifWidthTmp, errX, subPhiDifWidthErrTmp );
    subEtaGraph[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subEtaTmp, errX, subEtaErr );
    subEtaWidthGraph[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subEtaWidthTmp, errX, subEtaWidthErrTmp );
    
  }

  // print out the yields
  
  TCanvas* c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiGraph[i]->SetLineColor(i+1);
    leadPhiGraph[i]->SetMarkerStyle(29);
    leadPhiGraph[i]->SetMarkerSize(3);
    leadPhiGraph[i]->SetMarkerColor(i+1);
    leadPhiGraph[i]->SetTitle("Trigger Jet - #Delta#phi Fit Yield");
    leadPhiGraph[i]->GetXaxis()->SetTitle("p_{T}");
    leadPhiGraph[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    leadPhiGraph[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0)
      leadPhiGraph[i]->Draw();
    else
      leadPhiGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/leadphiyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiDifGraph[i]->SetLineColor(i+1);
    leadPhiDifGraph[i]->SetMarkerStyle(29);
    leadPhiDifGraph[i]->SetMarkerSize(3);
    leadPhiDifGraph[i]->SetMarkerColor(i+1);
    leadPhiDifGraph[i]->SetTitle("Trigger Jet - #eta Subtracted #Delta#phi Fit Yield");
    leadPhiDifGraph[i]->GetXaxis()->SetTitle("p_{T}");
    leadPhiDifGraph[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    leadPhiDifGraph[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0)
      leadPhiDifGraph[i]->Draw();
    else
      leadPhiDifGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/leadphidifyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadEtaGraph[i]->SetLineColor(i+1);
    leadEtaGraph[i]->SetMarkerStyle(29);
    leadEtaGraph[i]->SetMarkerSize(3);
    leadEtaGraph[i]->SetMarkerColor(i+1);
    leadEtaGraph[i]->SetTitle("Trigger Jet - #Delta#eta Fit Yield");
    leadEtaGraph[i]->GetXaxis()->SetTitle("p_{T}");
    leadEtaGraph[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    leadEtaGraph[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0)
      leadEtaGraph[i]->Draw();
    else
      leadEtaGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/leadetayield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiGraph[i]->SetLineColor(i+1);
    subPhiGraph[i]->SetMarkerStyle(29);
    subPhiGraph[i]->SetMarkerSize(3);
    subPhiGraph[i]->SetMarkerColor(i+1);
    subPhiGraph[i]->SetTitle("Recoil Jet - #Delta#phi Fit Yield");
    subPhiGraph[i]->GetXaxis()->SetTitle("p_{T}");
    subPhiGraph[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    subPhiGraph[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0)
      subPhiGraph[i]->Draw();
    else
      subPhiGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/subphiyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiDifGraph[i]->SetLineColor(i+1);
    subPhiDifGraph[i]->SetMarkerStyle(29);
    subPhiDifGraph[i]->SetMarkerSize(3);
    subPhiDifGraph[i]->SetMarkerColor(i+1);
    subPhiDifGraph[i]->SetTitle("Recoil Jet - #eta Subtracted #Delta#phi Fit Yield");
    subPhiDifGraph[i]->GetXaxis()->SetTitle("p_{T}");
    subPhiDifGraph[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    subPhiDifGraph[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0)
      subPhiDifGraph[i]->Draw();
    else
      subPhiDifGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/subphidifyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subEtaGraph[i]->SetLineColor(i+1);
    subEtaGraph[i]->SetMarkerStyle(29);
    subEtaGraph[i]->SetMarkerSize(3);
    subEtaGraph[i]->SetMarkerColor(i+1);
    subEtaGraph[i]->SetTitle("Recoil Jet - #Delta#eta Fit Yield");
    subEtaGraph[i]->GetXaxis()->SetTitle("p_{T}");
    subEtaGraph[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    subEtaGraph[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0) {
      subEtaGraph[i]->Draw();
    }
    else {
      subEtaGraph[i]->Draw("SAME");
    }
  }
  c1->SaveAs("tmp/subetayield.pdf");
  
  // now print out widths
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiWidthGraph[i]->SetLineColor(i+1);
    leadPhiWidthGraph[i]->SetMarkerStyle(29);
    leadPhiWidthGraph[i]->SetMarkerSize(3);
    leadPhiWidthGraph[i]->SetMarkerColor(i+1);
    leadPhiWidthGraph[i]->SetTitle("Trigger Jet #Delta#phi - Gaussian Widths");
    leadPhiWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    leadPhiWidthGraph[i]->GetYaxis()->SetTitle("Width");
    leadPhiWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 4);
    if ( i == 0)
      leadPhiWidthGraph[i]->Draw();
    else
      leadPhiWidthGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/leadphiwidth.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiDifWidthGraph[i]->SetLineColor(i+1);
    leadPhiDifWidthGraph[i]->SetMarkerStyle(29);
    leadPhiDifWidthGraph[i]->SetMarkerSize(3);
    leadPhiDifWidthGraph[i]->SetMarkerColor(i+1);
    leadPhiDifWidthGraph[i]->SetTitle("Trigger Jet #Delta#eta Subtracted #Delta#phi - Gaussian Widths");
    leadPhiDifWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    leadPhiDifWidthGraph[i]->GetYaxis()->SetTitle("Width");
    leadPhiDifWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 4);
    if ( i == 0)
      leadPhiDifWidthGraph[i]->Draw();
    else
      leadPhiDifWidthGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/leadphidifwidth.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadEtaWidthGraph[i]->SetLineColor(i+1);
    leadEtaWidthGraph[i]->SetMarkerStyle(29);
    leadEtaWidthGraph[i]->SetMarkerSize(3);
    leadEtaWidthGraph[i]->SetMarkerColor(i+1);
    leadEtaWidthGraph[i]->SetTitle("Trigger Jet #Delta#eta - Gaussian Widths");
    leadEtaWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    leadEtaWidthGraph[i]->GetYaxis()->SetTitle("Width");
    leadEtaWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 4);
    if ( i == 0)
      leadEtaWidthGraph[i]->Draw();
    else
      leadEtaWidthGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/leadetawidth.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiWidthGraph[i]->SetLineColor(i+1);
    subPhiWidthGraph[i]->SetMarkerStyle(29);
    subPhiWidthGraph[i]->SetMarkerSize(3);
    subPhiWidthGraph[i]->SetMarkerColor(i+1);
    subPhiWidthGraph[i]->SetTitle("Recoil Jet #Delta#phi - Gaussian Widths");
    subPhiWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    subPhiWidthGraph[i]->GetYaxis()->SetTitle("Width");
    subPhiWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 4);
    if ( i == 0)
      subPhiWidthGraph[i]->Draw();
    else
      subPhiWidthGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/subphiwidth.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiDifWidthGraph[i]->SetLineColor(i+1);
    subPhiDifWidthGraph[i]->SetMarkerStyle(29);
    subPhiDifWidthGraph[i]->SetMarkerSize(3);
    subPhiDifWidthGraph[i]->SetMarkerColor(i+1);
    subPhiDifWidthGraph[i]->SetTitle("Recoil Jet #Delta#eta Subtracted #Delta#phi - Gaussian Widths");
    subPhiDifWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    subPhiDifWidthGraph[i]->GetYaxis()->SetTitle("Width");
    subPhiDifWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 4);
    if ( i == 0)
      subPhiDifWidthGraph[i]->Draw();
    else
      subPhiDifWidthGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/subphidifwidth.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subEtaWidthGraph[i]->SetLineColor(i+1);
    subEtaWidthGraph[i]->SetMarkerStyle(29);
    subEtaWidthGraph[i]->SetMarkerSize(3);
    subEtaWidthGraph[i]->SetMarkerColor(i+1);
    subEtaWidthGraph[i]->SetTitle("Recoil Jet #Delta#eta - Gaussian Widths");
    subEtaWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    subEtaWidthGraph[i]->GetYaxis()->SetTitle("Width");
    subEtaWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 4);
    if ( i == 0)
      subEtaWidthGraph[i]->Draw();
    else
      subEtaWidthGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/subetawidth.pdf");

  
  // now put leading and subleading yields on their own plots
  for ( int i = 0; i < nFiles; ++i ) {
    c1 = new TCanvas();
    leadPhiGraph[i]->SetLineColor(1);
    leadPhiGraph[i]->SetMarkerStyle(29);
    leadPhiGraph[i]->SetMarkerSize(2);
    leadPhiGraph[i]->SetMarkerColor(1);
    leadPhiDifGraph[i]->SetLineColor(2);
    leadPhiDifGraph[i]->SetMarkerStyle(20);
    leadPhiDifGraph[i]->SetMarkerSize(2);
    leadPhiDifGraph[i]->SetMarkerColor(2);
    leadEtaGraph[i]->SetLineColor(3);
    leadEtaGraph[i]->SetMarkerStyle(21);
    leadEtaGraph[i]->SetMarkerSize(2);
    leadEtaGraph[i]->SetMarkerColor(3);
    
    leadPhiGraph[i]->Draw();
    leadPhiDifGraph[i]->Draw("SAME");
    leadEtaGraph[i]->Draw("SAME");
    
    TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);
    leg->SetHeader("Trigger Jet Yields"); // option "C" allows to center the header
    leg->AddEntry(leadPhiGraph[i],"#Delta#phi","lep");
    leg->AddEntry(leadPhiDifGraph[i],"#Delta#phi #Delta#eta subtracted","lep");
    leg->AddEntry(leadEtaGraph[i],"#Delta#eta","lep");
    leg->Draw();
    
    std::string graphOutName = "tmp/graph_out_"+analysisNames[i]+"_lead.pdf";
    
    c1->SaveAs( graphOutName.c_str() );
    
    c1 = new TCanvas();
    subPhiGraph[i]->SetLineColor(1);
    subPhiGraph[i]->SetMarkerStyle(29);
    subPhiGraph[i]->SetMarkerSize(2);
    subPhiGraph[i]->SetMarkerColor(1);
    subPhiDifGraph[i]->SetLineColor(2);
    subPhiDifGraph[i]->SetMarkerStyle(20);
    subPhiDifGraph[i]->SetMarkerSize(2);
    subPhiDifGraph[i]->SetMarkerColor(2);
    subEtaGraph[i]->SetLineColor(3);
    subEtaGraph[i]->SetMarkerStyle(21);
    subEtaGraph[i]->SetMarkerSize(2);
    subEtaGraph[i]->SetMarkerColor(3);
    
    subPhiGraph[i]->Draw();
    subPhiDifGraph[i]->Draw("SAME");
    subEtaGraph[i]->Draw("SAME");
    
    leg = new TLegend(0.65,0.7,0.9,0.9);
    leg->SetHeader("Recoil Jet Yields"); // option "C" allows to center the header
    leg->AddEntry(subPhiGraph[i],"#Delta#phi","lep");
    leg->AddEntry(subPhiDifGraph[i],"#Delta#phi #Delta#eta subtracted","lep");
    leg->AddEntry(subEtaGraph[i],"#Delta#eta","lep");
    leg->Draw();

    graphOutName = "tmp/graph_out_"+analysisNames[i]+"_sub.pdf";

    c1->SaveAs( graphOutName.c_str() );
    
    // joern asked for dphi without subtraction overlayed
    graphOutName = "tmp/dphi_trig_recoil_yield_"+analysisNames[i]+".pdf";
    c1 = new TCanvas();
    leadPhiGraph[i]->SetTitle("");
    leadPhiGraph[i]->SetLineColor(1);
    leadPhiGraph[i]->SetMarkerColor(1);
    subPhiGraph[i]->SetLineColor(7);
    subPhiGraph[i]->SetMarkerColor(7);
    subPhiGraph[i]->SetMarkerStyle(20);
    
    leadPhiGraph[i]->Draw();
    subPhiGraph[i]->Draw("SAME");
    
    leg = new TLegend(0.65,0.7,0.9,0.9);
    leg->SetHeader("#Delta#phi Jet Yields");
    leg->AddEntry(leadPhiGraph[i],"Trigger","lep");
    leg->AddEntry(subPhiGraph[i],"Recoil","lep");
    leg->Draw();
    
    c1->SaveAs( graphOutName.c_str() );
  }
  
  // testing pt averaging
  
  double total;
  double numberOfEntries;
  for ( int i = 0; i < nPtBins; ++i ) {
    total = 0;
    numberOfEntries = 0;
    for ( int j = ptBinLo[i]; j <= ptBinHi[i]; ++j ) {
      numberOfEntries += recombinedPtLead[0]->GetBinContent(i);
      total += recombinedPtLead[0]->GetBinContent(j)*recombinedPtLead[0]->GetBinCenter(j);
      std::cout<<"j: "<<j<<std::endl;
      std::cout<<"content: "<< recombinedPtLead[0]->GetBinContent(i) << std::endl;
      std::cout<<"bin center: "<< recombinedPtLead[0]->GetBinCenter(j) << std::endl;
    }
    std::cout<<"total: "<< total << std::endl;
    std::cout<<"bin: "<<i<<" average: "<< total/numberOfEntries<<std::endl;
  }

  
  return 0;
}
