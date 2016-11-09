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
  
  // define if the analyses are split or not
  bool ajSplit[nFiles];
  
  // Load in the histograms
  TH3D* nEvents[ nFiles ];
  TH1D* hVz[ nFiles ];
  TH3D* corrHist[ nFiles ];
  TH3D* mixHist[ nFiles ];
  std::vector<std::vector<std::vector<TH3D*> > > corrCentVzLarge;
  std::vector<std::vector<std::vector<TH3D*> > > corrCentVzSmall;
  std::vector<std::vector<std::vector<TH3D*> > > mixCentVz;
  std::vector<std::vector<std::vector<TH3D*> > > subCentVzLarge;
  std::vector<std::vector<std::vector<TH3D*> > > subCentVzSmall;
  std::vector<std::vector<std::vector<TH3D*> > > mixSubCentVz;
  std::vector<TH1D*> recombinedPtLead;
  std::vector<TH1D*> recombinedPtSub;
  corrCentVzLarge.resize( nFiles );
  corrCentVzSmall.resize( nFiles );
  mixCentVz.resize( nFiles );
  subCentVzLarge.resize( nFiles );
  subCentVzSmall.resize( nFiles );
  mixSubCentVz.resize( nFiles );
  recombinedPtLead.resize( nFiles );
  recombinedPtSub.resize( nFiles );
  
  // while resizing the 3dim arrays we can also set up the pt spectra histograms
  for ( int i = 0; i < nFiles; ++i ) {
    corrCentVzLarge[i].resize( corrAnalysis::binsCentrality );
    corrCentVzSmall[i].resize( corrAnalysis::binsCentrality );
    mixCentVz[i].resize( corrAnalysis::binsCentrality );
    subCentVzLarge[i].resize( corrAnalysis::binsCentrality );
    subCentVzSmall[i].resize( corrAnalysis::binsCentrality );
    mixSubCentVz[i].resize( corrAnalysis::binsCentrality );
    
    std::string ptLeadName = analysisNames[i] + "_pt_lead";
    std::string ptSubName = analysisNames[i] + "_pt_sub";
    
    recombinedPtLead[i] = new TH1D( ptLeadName.c_str(), "p_{T} Spectrum Trigger Jet", corrAnalysis::binsPt, corrAnalysis::ptLowEdge, corrAnalysis::ptHighEdge );
    recombinedPtSub[i] = new TH1D( ptSubName.c_str(), "p_{T} Spectrum Recoil Jet", corrAnalysis::binsPt, corrAnalysis::ptLowEdge, corrAnalysis::ptHighEdge );

    
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      corrCentVzLarge[i][j].resize( corrAnalysis::binsVz );
      corrCentVzSmall[i][j].resize( corrAnalysis::binsVz );
      mixCentVz[i][j].resize( corrAnalysis::binsVz );
      subCentVzLarge[i][j].resize( corrAnalysis::binsVz );
      subCentVzSmall[i][j].resize( corrAnalysis::binsVz );
      mixSubCentVz[i][j].resize( corrAnalysis::binsVz );
    }
  }
  
  // while loading the 3dim vectors we will also get the projections for pt so that the axis resizing in the next step doesnt affect it
  for ( int i = 0; i < nFiles; ++i ) {
    std::cout<<"loading histograms - file: "<< i <<std::endl;
    std::string neventsBaseName = "nevents_"; neventsBaseName += analysisNames[i];
    std::string hvzBaseName = "hvz_"; hvzBaseName += analysisNames[i];
    std::string corrhistBaseName = "corrHist_"; corrhistBaseName += analysisNames[i];
    std::string mixhistBaseName = "mixHist_"; mixhistBaseName += analysisNames[i];

    nEvents[i] = (TH3D*) corrFiles[i]->Get( "nevents" );
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
        // assuming it is a non-split file
        std::string corrDifInitName = "lead_cent_"; corrDifInitName += patch::to_string(j);
        std::string mixDifInitName = "mix_lead_cent_"; mixDifInitName += patch::to_string(j);
        std::string subDifInitName = "sub_cent_"; subDifInitName += patch::to_string(j);
        std::string mixSubDifInitName = "mix_sub_cent_";
        mixSubDifInitName += patch::to_string(j);
        
        corrDifInitName += "_vz_"; corrDifInitName += patch::to_string(k);
        mixDifInitName += "_vz_"; mixDifInitName += patch::to_string(k);
        subDifInitName += "_vz_"; subDifInitName += patch::to_string(k);
        mixSubDifInitName += "_vz_"; mixSubDifInitName += patch::to_string(k);
        
        // make the new histogram name
        std::string corrDifBaseName = "corr_file_"; corrDifBaseName += patch::to_string(i);
        std::string mixDifBaseName = "mix_file_"; mixDifBaseName += patch::to_string(i);
        std::string subDifBaseName = "sub_file_"; subDifBaseName += patch::to_string(i);
        std::string mixSubDifBaseName = "mix_file_"; mixSubDifBaseName += patch::to_string(i);
        
        corrDifBaseName += "_cent_"; corrDifBaseName += patch::to_string(j);
        corrDifBaseName += "_vz_"; corrDifBaseName += patch::to_string(k);
        
        mixDifBaseName += "_cent_"; mixDifBaseName += patch::to_string(j);
        mixDifBaseName += "_vz_"; mixDifBaseName += patch::to_string(k);
        
        subDifBaseName += "_cent_"; subDifBaseName += patch::to_string(j);
        subDifBaseName += "_vz_"; subDifBaseName += patch::to_string(k);
        
        mixSubDifBaseName += "_cent_"; mixSubDifBaseName += patch::to_string(j);
        mixSubDifBaseName += "_vz_"; mixSubDifBaseName += patch::to_string(k);

        // get the histograms
        if ( corrFiles[i]->Get( corrDifInitName.c_str() ) ) {
          corrCentVzLarge[i][j][k] = (TH3D*) corrFiles[i]->Get( corrDifInitName.c_str() );
          corrCentVzLarge[i][j][k]->SetName( corrDifBaseName.c_str() );
          ajSplit[i] = false;
        }
        else {
          std::string smallName = "small_" + corrDifInitName;
          std::string largeName = "large_" + corrDifInitName;
          std::string smallDifName = "small_" + corrDifBaseName;
          std::string largeDifName = "large_" + corrDifBaseName;

          corrCentVzLarge[i][j][k] = (TH3D*) corrFiles[i]->Get( largeName.c_str() );
          corrCentVzLarge[i][j][k]->SetName( largeDifName.c_str() );
          corrCentVzSmall[i][j][k] = (TH3D*) corrFiles[i]->Get( smallName.c_str() );
          corrCentVzSmall[i][j][k]->SetName( smallDifName.c_str() );
          ajSplit[i] = true;
          
        }

        mixCentVz[i][j][k] = (TH3D*) mixFiles[i]->Get( mixDifInitName.c_str() );
        mixCentVz[i][j][k]->SetName( mixDifBaseName.c_str() );
        if ( corrFiles[i]->Get( subDifInitName.c_str() ) ) {
          subCentVzLarge[i][j][k] = (TH3D*) corrFiles[i]->Get( subDifInitName.c_str() );
          subCentVzLarge[i][j][k]->SetName( subDifBaseName.c_str() );
        }
        else {
          std::string smallName = "small_" + subDifInitName;
          std::string largeName = "large_" + subDifInitName;
          std::string smallDifName = "small_" + subDifBaseName;
          std::string largeDifName = "large_" + subDifBaseName;
          
          subCentVzLarge[i][j][k] = (TH3D*) corrFiles[i]->Get( largeName.c_str() );
          subCentVzLarge[i][j][k]->SetName( largeDifName.c_str() );
          subCentVzSmall[i][j][k] = (TH3D*) corrFiles[i]->Get( smallName.c_str() );
          subCentVzSmall[i][j][k]->SetName( smallDifName.c_str() );
        }
        mixSubCentVz[i][j][k] = (TH3D*) mixFiles[i]->Get( mixSubDifInitName.c_str() );
        mixSubCentVz[i][j][k]->SetName( mixSubDifBaseName.c_str() );
        
        // now we can get the pt spectrum as well
        recombinedPtLead[i]->Add( (TH1D*) corrCentVzLarge[i][j][k]->Project3D("Z") );
        recombinedPtSub[i]->Add( (TH1D*) subCentVzLarge[i][j][k]->Project3D("Z") );
        
        if ( corrCentVzSmall[i][j][k] )
          recombinedPtLead[i]->Add( (TH1D*) corrCentVzSmall[i][j][k]->Project3D("Z") );
        if ( subCentVzSmall[i][j][k] )
          recombinedPtSub[i]->Add( (TH1D*) subCentVzSmall[i][j][k]->Project3D("Z") );
        
      }
  }
  
  std::cout<<"finished loading all histograms"<<std::endl;
  
  // setup for 2d projections along pt axis
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > corrCentVzPtLarge;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > corrCentVzPtSmall;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > mixCentVzPt;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > subCentVzPtLarge;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > subCentVzPtSmall;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > mixSubCentVzPt;
  corrCentVzPtLarge.resize( nFiles );
  corrCentVzPtSmall.resize( nFiles );
  mixCentVzPt.resize( nFiles );
  subCentVzPtLarge.resize( nFiles );
  subCentVzPtSmall.resize( nFiles );
  mixSubCentVzPt.resize( nFiles );
  for ( int i = 0; i < nFiles; ++i ) {
    corrCentVzPtLarge[i].resize( corrAnalysis::binsCentrality );
    corrCentVzPtSmall[i].resize( corrAnalysis::binsCentrality );
    mixCentVzPt[i].resize( corrAnalysis::binsCentrality );
    subCentVzPtLarge[i].resize( corrAnalysis::binsCentrality );
    subCentVzPtSmall[i].resize( corrAnalysis::binsCentrality );
    mixSubCentVzPt[i].resize( corrAnalysis::binsCentrality );
    
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      corrCentVzPtLarge[i][j].resize( corrAnalysis::binsVz );
      corrCentVzPtSmall[i][j].resize( corrAnalysis::binsVz );
      mixCentVzPt[i][j].resize( corrAnalysis::binsVz );
      subCentVzPtLarge[i][j].resize( corrAnalysis::binsVz );
      subCentVzPtSmall[i][j].resize( corrAnalysis::binsVz );
      mixSubCentVzPt[i][j].resize( corrAnalysis::binsVz );
      for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
        corrCentVzPtLarge[i][j][k].resize( nPtBins );
        corrCentVzPtSmall[i][j][k].resize( nPtBins );
        mixCentVzPt[i][j][k].resize( nPtBins );
        subCentVzPtLarge[i][j][k].resize( nPtBins );
        subCentVzPtSmall[i][j][k].resize( nPtBins );
        mixSubCentVzPt[i][j][k].resize( nPtBins );
      }
    }
  }

  // now get the pt projections
  for ( int i = 0; i < nFiles; ++i ) {
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      for ( int k = 0; k < corrAnalysis::binsVz; ++ k ) {
        for ( int l = 0; l < nPtBins; ++l ) {
          
          corrCentVzLarge[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          if ( corrCentVzSmall[i][j][k] )
            corrCentVzSmall[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          mixCentVz[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          subCentVzLarge[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          if ( subCentVzSmall[i][j][k] )
            subCentVzSmall[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          mixSubCentVz[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          
          corrCentVzPtLarge[i][j][k][l] = (TH2D*) ((TH2D*) corrCentVzLarge[i][j][k]->Project3D( "YX" ))->Clone();
          mixCentVzPt[i][j][k][l] = (TH2D*) ((TH2D*) mixCentVz[i][j][k]->Project3D( "YX" ))->Clone();
          subCentVzPtLarge[i][j][k][l] = (TH2D*) ((TH2D*) subCentVzLarge[i][j][k]->Project3D("YX"))->Clone();
          mixSubCentVzPt[i][j][k][l] = (TH2D*) ((TH2D*) mixSubCentVz[i][j][k]->Project3D("YX"))->Clone();
          
          if ( corrCentVzSmall[i][j][k] )
            corrCentVzPtSmall[i][j][k][l] = (TH2D*) ((TH2D*) corrCentVzSmall[i][j][k]->Project3D( "YX" ))->Clone();
          if ( subCentVzSmall[i][j][k] )
            subCentVzPtSmall[i][j][k][l] = (TH2D*) ((TH2D*) subCentVzSmall[i][j][k]->Project3D("YX"))->Clone();
        }
      }
    }
  }

  // TESTING PAST HERE
  // Averaging the event mixing over vz/cent
  // NEEDS TO BE UPDATED FOR UPDATED PT BINS
  
  // now to average over all vz/centrality bins
  // (weighted in vz but not centrality)
  std::vector<std::vector<TH2D*> > weightedMix;
  std::vector<std::vector<TH2D*> > weightedSub;
  weightedMix.resize( nFiles );
  weightedSub.resize( nFiles );
  for ( int i = 0; i < nFiles; ++i ) {
    
    weightedMix[i].resize( nPtBins );
    weightedSub[i].resize( nPtBins );
    
    for ( int l = 0; l < nPtBins; ++l ) {
      std::string weightedMixName = "ave_mix_file_"; weightedMixName += patch::to_string( i );
      weightedMixName += "_ptBin_"; weightedMixName += patch::to_string( l );
      std::string weightedSubName = "sub_mix_file_"; weightedSubName += patch::to_string( i );
      weightedSubName += "_ptBin_"; weightedSubName += patch::to_string( l );
      // create new histogram, add all appropriate vz/cent bins
      weightedMix[i][l] = new TH2D( weightedMixName.c_str(), weightedMixName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      weightedSub[i][l] = new TH2D( weightedSubName.c_str(), weightedSubName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
        for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
          
          if ( l <= 2 ) {
            weightedMix[i][l]->Add( mixCentVzPt[i][j][k][l] );
            weightedSub[i][l]->Add( mixSubCentVzPt[i][j][k][l] );
          }
          else {
            weightedMix[i][2]->Add( mixCentVzPt[i][j][k][l] );
            weightedSub[i][2]->Add( mixSubCentVzPt[i][j][k][l] );
          }
        }
      }
    }
  }
  
  // scale the weighted mixing histograms
  for ( int i = 0; i < nFiles; ++i ) {
    for ( int j = 0; j <= 2; ++j ) {
      weightedMix[i][j]->Scale( 1.0/weightedMix[i][j]->GetMaximum() );
      weightedSub[i][j]->Scale( 1.0/weightedSub[i][j]->GetMaximum() );
    }
  }
  
  // make the container for the recombined histograms
  std::vector<std::vector<TH2D*> > recombinedCorrLarge;
  std::vector<std::vector<TH2D*> > recombinedCorrSmall;
  std::vector<std::vector<TH2D*> > recombinedPreLarge;
  std::vector<std::vector<TH2D*> > recombinedPreSmall;
  std::vector<std::vector<TH2D*> > recombinedSubLarge;
  std::vector<std::vector<TH2D*> > recombinedSubSmall;
  std::vector<std::vector<TH2D*> > recombinedSubPreLarge;
  std::vector<std::vector<TH2D*> > recombinedSubPreSmall;
  
  recombinedCorrLarge.resize( nFiles );
  recombinedCorrSmall.resize( nFiles );
  recombinedPreLarge.resize( nFiles );
  recombinedPreSmall.resize( nFiles );
  recombinedSubLarge.resize( nFiles );
  recombinedSubSmall.resize( nFiles );
  recombinedSubPreLarge.resize( nFiles );
  recombinedSubPreSmall.resize( nFiles );
  
  for (int i = 0; i < nFiles; ++ i ) {

    recombinedCorrLarge[i].resize( nPtBins );
    recombinedCorrSmall[i].resize( nPtBins );
    recombinedPreLarge[i].resize( nPtBins );
    recombinedPreSmall[i].resize( nPtBins );
    recombinedSubLarge[i].resize( nPtBins );
    recombinedSubSmall[i].resize( nPtBins );
    recombinedSubPreLarge[i].resize( nPtBins );
    recombinedSubPreSmall[i].resize( nPtBins );

    
    for ( int l = 0; l < nPtBins; ++l ) {
      
      std::string corrName = analysisNames[i] + " " + ptBinString[l];
      std::string corrSmallName = "";
      std::string preName = "pre_" + analysisNames[i] + " " + ptBinString[l];
      std::string preSmallName = "";
      std::string subName = analysisNames[i] + "_sub " + ptBinString[l];
      std::string subSmallName = "";
      std::string subPreName = "pre_" + analysisNames[i] + "_sub " + ptBinString[l];
      std::string subPreSmallName = "";
      
      if ( ajSplit[i] ) {
        corrName = "large " + analysisNames[i] + " " + ptBinString[l];
        corrSmallName = "small " + analysisNames[i] + " " + ptBinString[l];
        subName = "large " + analysisNames[i] + "_sub " + ptBinString[l];
        subSmallName = "small " + analysisNames[i] + "_sub " + ptBinString[l];
        preName = "large pre_" + analysisNames[i] + " " + ptBinString[l];
        preSmallName = "small pre_" + analysisNames[i] + " " + ptBinString[l];
        subPreName = "large pre_" + analysisNames[i] + "_sub " + ptBinString[l];
        subPreSmallName = "small pre_" + analysisNames[i] + "_sub " + ptBinString[l];
      }
      
      recombinedCorrLarge[i][l] = new TH2D( corrName.c_str(), corrName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      recombinedPreLarge[i][l] = new TH2D( preName.c_str(), preName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      recombinedSubLarge[i][l] = new TH2D( subName.c_str(), subName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      recombinedSubPreLarge[i][l] = new TH2D( subPreName.c_str(), subPreName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      if ( ajSplit[i] ) {
        recombinedCorrSmall[i][l] = new TH2D( corrSmallName.c_str(), corrSmallName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
        
        recombinedSubSmall[i][l] = new TH2D( subSmallName.c_str(), subSmallName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
        
        recombinedPreSmall[i][l] = new TH2D( preSmallName.c_str(), preSmallName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
        
        recombinedSubPreSmall[i][l] = new TH2D( subPreSmallName.c_str(), subPreSmallName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
        
        
      }
      
      if ( l <= 2 ) {
      
        for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
          for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
            if ( weightedMix[i][l]->GetEntries() != 0 && corrCentVzPtLarge[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedPreLarge[i][l]->Add( corrCentVzPtLarge[i][j][k][l] );
              
              corrCentVzPtLarge[i][j][k][l]->Divide( weightedMix[i][l] );
            
              recombinedCorrLarge[i][l]->Add( corrCentVzPtLarge[i][j][k][l] );
            }
            if ( ajSplit[i] && weightedMix[i][l]->GetEntries() != 0 && corrCentVzPtSmall[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedPreSmall[i][l]->Add( corrCentVzPtSmall[i][j][k][l] );
              
              corrCentVzPtSmall[i][j][k][l]->Divide( weightedMix[i][l] );
              
              recombinedCorrSmall[i][l]->Add( corrCentVzPtSmall[i][j][k][l] );
            }
            if ( subCentVzPtLarge[i][j][k][l]->GetEntries() != 0 && weightedSub[i][l]->GetEntries() != 0 ) {
              
              recombinedSubPreLarge[i][l]->Add( subCentVzPtLarge[i][j][k][l] );
              
              subCentVzPtLarge[i][j][k][l]->Divide( weightedSub[i][l] );
              
              recombinedSubLarge[i][l]->Add( subCentVzPtLarge[i][j][k][l] );
              
            }
            if ( ajSplit[i] && weightedMix[i][l]->GetEntries() != 0 && subCentVzPtSmall[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedSubPreSmall[i][l]->Add( subCentVzPtSmall[i][j][k][l] );
              
              subCentVzPtSmall[i][j][k][l]->Divide( weightedMix[i][l] );
              
              recombinedSubSmall[i][l]->Add( subCentVzPtSmall[i][j][k][l] );
            }
          }
        }
      }
      
      else {
        for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
          for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
            if ( weightedMix[i][2]->GetEntries() != 0 && corrCentVzPtLarge[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedPreLarge[i][l]->Add( corrCentVzPtLarge[i][j][k][l] );
              
              corrCentVzPtLarge[i][j][k][l]->Divide( weightedMix[i][2] );
              
              recombinedCorrLarge[i][l]->Add( corrCentVzPtLarge[i][j][k][l] );
            }
            if ( ajSplit[i] && weightedMix[i][2]->GetEntries() != 0 && corrCentVzPtSmall[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedPreSmall[i][l]->Add( corrCentVzPtSmall[i][j][k][l] );
              
              corrCentVzPtSmall[i][j][k][l]->Divide( weightedMix[i][2] );
              
              recombinedCorrSmall[i][l]->Add( corrCentVzPtSmall[i][j][k][l] );
            }
            if ( weightedSub[i][2]->GetEntries() != 0 && subCentVzPtLarge[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedSubPreLarge[i][l]->Add( subCentVzPtLarge[i][j][k][l] );
              
              subCentVzPtLarge[i][j][k][l]->Divide( weightedSub[i][2] );
              
              recombinedSubLarge[i][l]->Add( subCentVzPtLarge[i][j][k][l] );
            }
            if ( ajSplit[i] && weightedMix[i][2]->GetEntries() != 0 && subCentVzPtSmall[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedSubPreSmall[i][l]->Add( subCentVzPtSmall[i][j][k][l] );
              
              subCentVzPtSmall[i][j][k][l]->Divide( weightedMix[i][2] );
              
              recombinedSubSmall[i][l]->Add( subCentVzPtSmall[i][j][k][l] );
            }

          }
        }
      }
    }
  }
  
  // Testing Aj Small - Aj Large
  // ***************************
  std::vector<std::vector<TH1D*> > triggerSmallAjdPhi( nFiles );
  std::vector<std::vector<TH1D*> > triggerLargeAjdPhi( nFiles );
  std::vector<std::vector<TH1D*> > recoilSmallAjdPhi( nFiles );
  std::vector<std::vector<TH1D*> > recoilLargeAjdPhi( nFiles );
  // fits
  std::vector<std::vector<TF1*> > triggerSmallAjdPhiFit( nFiles );
  std::vector<std::vector<TF1*> > triggerLargeAjdPhiFit( nFiles );
  std::vector<std::vector<TF1*> > recoilSmallAjdPhiFit( nFiles );
  std::vector<std::vector<TF1*> > recoilLargeAjdPhiFit( nFiles );
  
  // subtracted histograms
  std::vector<std::vector<TH1D*> > triggerSubtracted( nFiles );
  std::vector<std::vector<TH1D*> > recoilSubtracted( nFiles );
  
  std::vector<std::vector<TF1*> > triggerSubFit( nFiles );
  std::vector<std::vector<TF1*> > recoilSubFit( nFiles );
  
  // and we'll need some arrays for the widths and yields
  std::vector<std::vector<double> > triggerLargeYieldNear( nFiles );
  std::vector<std::vector<double> > triggerLargeYieldNearErr( nFiles );
  std::vector<std::vector<double> > triggerLargeWidthNear( nFiles );
  std::vector<std::vector<double> > triggerLargeWidthNearErr( nFiles );
  std::vector<std::vector<double> > triggerSmallYieldNear( nFiles );
  std::vector<std::vector<double> > triggerSmallYieldNearErr( nFiles );
  std::vector<std::vector<double> > triggerSmallWidthNear( nFiles );
  std::vector<std::vector<double> > triggerSmallWidthNearErr( nFiles );
  std::vector<std::vector<double> > recoilLargeYieldNear( nFiles );
  std::vector<std::vector<double> > recoilLargeYieldNearErr( nFiles );
  std::vector<std::vector<double> > recoilLargeWidthNear( nFiles );
  std::vector<std::vector<double> > recoilLargeWidthNearErr( nFiles );
  std::vector<std::vector<double> > recoilSmallYieldNear( nFiles );
  std::vector<std::vector<double> > recoilSmallYieldNearErr( nFiles );
  std::vector<std::vector<double> > recoilSmallWidthNear( nFiles );
  std::vector<std::vector<double> > recoilSmallWidthNearErr( nFiles );
  
  std::vector<std::vector<double> > triggerLargeYieldAway( nFiles );
  std::vector<std::vector<double> > triggerLargeYieldAwayErr( nFiles );
  std::vector<std::vector<double> > triggerLargeWidthAway( nFiles );
  std::vector<std::vector<double> > triggerLargeWidthAwayErr( nFiles );
  std::vector<std::vector<double> > triggerSmallYieldAway( nFiles );
  std::vector<std::vector<double> > triggerSmallYieldAwayErr( nFiles );
  std::vector<std::vector<double> > triggerSmallWidthAway( nFiles );
  std::vector<std::vector<double> > triggerSmallWidthAwayErr( nFiles );
  std::vector<std::vector<double> > recoilLargeYieldAway( nFiles );
  std::vector<std::vector<double> > recoilLargeYieldAwayErr( nFiles );
  std::vector<std::vector<double> > recoilLargeWidthAway( nFiles );
  std::vector<std::vector<double> > recoilLargeWidthAwayErr( nFiles );
  std::vector<std::vector<double> > recoilSmallYieldAway( nFiles );
  std::vector<std::vector<double> > recoilSmallYieldAwayErr( nFiles );
  std::vector<std::vector<double> > recoilSmallWidthAway( nFiles );
  std::vector<std::vector<double> > recoilSmallWidthAwayErr( nFiles );
  
  std::vector<std::vector<double> > triggerSubYieldNear( nFiles );
  std::vector<std::vector<double> > triggerSubYieldNearErr( nFiles );
  std::vector<std::vector<double> > triggerSubWidthNear( nFiles );
  std::vector<std::vector<double> > triggerSubWidthNearErr( nFiles );
  
  std::vector<std::vector<double> > triggerSubYieldAway( nFiles );
  std::vector<std::vector<double> > triggerSubYieldAwayErr( nFiles );
  std::vector<std::vector<double> > triggerSubWidthAway( nFiles );
  std::vector<std::vector<double> > triggerSubWidthAwayErr( nFiles );
  
  std::vector<std::vector<double> > recoilSubYieldNear( nFiles );
  std::vector<std::vector<double> > recoilSubYieldNearErr( nFiles );
  std::vector<std::vector<double> > recoilSubWidthNear( nFiles );
  std::vector<std::vector<double> > recoilSubWidthNearErr( nFiles );

  std::vector<std::vector<double> > recoilSubYieldAway( nFiles );
  std::vector<std::vector<double> > recoilSubYieldAwayErr( nFiles );
  std::vector<std::vector<double> > recoilSubWidthAway( nFiles );
  std::vector<std::vector<double> > recoilSubWidthAwayErr( nFiles );
  
  for ( int i = 0; i < nFiles; ++i ) {
    
    triggerSmallAjdPhi[i].resize( nPtBins );
    triggerLargeAjdPhi[i].resize( nPtBins );
    recoilSmallAjdPhi[i].resize( nPtBins );
    recoilLargeAjdPhi[i].resize( nPtBins );
    
    triggerSmallAjdPhiFit[i].resize( nPtBins );
    triggerLargeAjdPhiFit[i].resize( nPtBins );
    recoilSmallAjdPhiFit[i].resize( nPtBins );
    recoilLargeAjdPhiFit[i].resize( nPtBins );
    
    triggerSubtracted[i].resize( nPtBins );
    recoilSubtracted[i].resize( nPtBins );
    
    triggerSubFit[i].resize( nPtBins );
    recoilSubFit[i].resize( nPtBins );
    
    triggerLargeYieldNear[i].resize( nPtBins );
    triggerLargeYieldNearErr[i].resize( nPtBins );
    triggerLargeWidthNear[i].resize( nPtBins );
    triggerLargeWidthNearErr[i].resize( nPtBins );
    triggerSmallYieldNear[i].resize( nPtBins );
    triggerSmallYieldNearErr[i].resize( nPtBins );
    triggerSmallWidthNear[i].resize( nPtBins );
    triggerSmallWidthNearErr[i].resize( nPtBins );
    recoilLargeYieldNear[i].resize( nPtBins );
    recoilLargeYieldNearErr[i].resize( nPtBins );
    recoilLargeWidthNear[i].resize( nPtBins );
    recoilLargeWidthNearErr[i].resize( nPtBins );
    recoilSmallYieldNear[i].resize( nPtBins );
    recoilSmallYieldNearErr[i].resize( nPtBins );
    recoilSmallWidthNear[i].resize( nPtBins );
    recoilSmallWidthNearErr[i].resize( nPtBins );
    
    triggerLargeYieldAway[i].resize( nPtBins );
    triggerLargeYieldAwayErr[i].resize( nPtBins );
    triggerLargeWidthAway[i].resize( nPtBins );
    triggerLargeWidthAwayErr[i].resize( nPtBins );
    triggerSmallYieldAway[i].resize( nPtBins );
    triggerSmallYieldAwayErr[i].resize( nPtBins );
    triggerSmallWidthAway[i].resize( nPtBins );
    triggerSmallWidthAwayErr[i].resize( nPtBins );
    recoilLargeYieldAway[i].resize( nPtBins );
    recoilLargeYieldAwayErr[i].resize( nPtBins );
    recoilLargeWidthAway[i].resize( nPtBins );
    recoilLargeWidthAwayErr[i].resize( nPtBins );
    recoilSmallYieldAway[i].resize( nPtBins );
    recoilSmallYieldAwayErr[i].resize( nPtBins );
    recoilSmallWidthAway[i].resize( nPtBins );
    recoilSmallWidthAwayErr[i].resize( nPtBins );
    
    triggerSubYieldNear[i].resize( nPtBins );
    triggerSubYieldNearErr[i].resize( nPtBins );
    triggerSubWidthNear[i].resize( nPtBins );
    triggerSubWidthNearErr[i].resize( nPtBins );
    
    triggerSubYieldAway[i].resize( nPtBins );
    triggerSubYieldAwayErr[i].resize( nPtBins );
    triggerSubWidthAway[i].resize( nPtBins );
    triggerSubWidthAwayErr[i].resize( nPtBins );

    recoilSubYieldNear[i].resize( nPtBins );
    recoilSubYieldNearErr[i].resize( nPtBins );
    recoilSubWidthNear[i].resize( nPtBins );
    recoilSubWidthNearErr[i].resize( nPtBins );
    
    recoilSubYieldAway[i].resize( nPtBins );
    recoilSubYieldAwayErr[i].resize( nPtBins );
    recoilSubWidthAway[i].resize( nPtBins );
    recoilSubWidthAwayErr[i].resize( nPtBins );
    
    if ( ajSplit[i] ) {
      // for normalization
      double ajHighCount = nEvents[i]->Integral( 1, 1, 1, corrAnalysis::binsCentrality, 1, corrAnalysis::binsVz );
      double ajLowCount = nEvents[i]->Integral( 2, 2, 1, corrAnalysis::binsCentrality, 1, corrAnalysis::binsVz );
      
      for ( int j = 0; j < nPtBins; ++j ) {
        
        triggerSmallAjdPhi[i][j] = (TH1D*) ((TH1D*) recombinedPreSmall[i][j]->ProjectionY())->Clone();
        triggerLargeAjdPhi[i][j] = (TH1D*) ((TH1D*) recombinedPreLarge[i][j]->ProjectionY())->Clone();
        recoilSmallAjdPhi[i][j] = (TH1D*) ((TH1D*) recombinedSubPreSmall[i][j]->ProjectionY())->Clone();
        recoilLargeAjdPhi[i][j] = (TH1D*) ((TH1D*) recombinedSubPreLarge[i][j]->ProjectionY())->Clone();
        
        // Fit and subtract
        std::string triggerLargeName = "tmp_trigger_large_"; triggerLargeName += patch::to_string(i); triggerLargeName += patch::to_string(j);
        std::string triggerSmallName = "tmp_trigger_small_"; triggerSmallName += patch::to_string(i); triggerSmallName += patch::to_string(j);
        std::string recoilSmallName = "tmp_recoil_small_"; recoilSmallName += patch::to_string(i); recoilSmallName += patch::to_string(j);
        std::string recoilLargeName = "tmp_recoil_large_"; recoilLargeName += patch::to_string(i); recoilLargeName += patch::to_string(j);
        std::string phiForm = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)+[4]*exp(-0.5*((x-[5])/[6])**2)";
        double phiMin = -corrAnalysis::pi/2.0;
        double phiMax = 3.0*corrAnalysis::pi/2.0;

        TF1* triggerLargeTmpFit = new TF1(triggerLargeName.c_str(), phiForm.c_str(), phiMin, phiMax);
        triggerLargeTmpFit->FixParameter( 2, 0 );
        triggerLargeTmpFit->FixParameter( 5, corrAnalysis::pi );
        triggerLargeTmpFit->SetParameter( 3, 0.2 );
        triggerLargeTmpFit->SetParameter( 6, 0.2 );
        TF1* triggerSmallTmpFit = new TF1(triggerSmallName.c_str(), phiForm.c_str(), phiMin, phiMax);
        triggerSmallTmpFit->FixParameter( 2, 0 );
        triggerSmallTmpFit->FixParameter( 5, corrAnalysis::pi );
        triggerSmallTmpFit->SetParameter( 3, 0.2 );
        triggerSmallTmpFit->SetParameter( 6, 0.2 );
        TF1* recoilLargeTmpFit = new TF1(recoilLargeName.c_str(), phiForm.c_str(), phiMin, phiMax);
        recoilLargeTmpFit->FixParameter( 2, 0 );
        recoilLargeTmpFit->FixParameter( 5, corrAnalysis::pi );
        recoilLargeTmpFit->SetParameter( 3, 0.2 );
        recoilLargeTmpFit->SetParameter( 6, 0.2 );
        TF1* recoilSmallTmpFit = new TF1(recoilSmallName.c_str(), phiForm.c_str(), phiMin, phiMax);
        recoilSmallTmpFit->FixParameter( 2, 0 );
        recoilSmallTmpFit->FixParameter( 5, corrAnalysis::pi );
        recoilSmallTmpFit->SetParameter( 3, 0.2 );
        recoilSmallTmpFit->SetParameter( 6, 0.2 );
        
        
        // fit
        triggerSmallAjdPhi[i][j]->Fit( triggerSmallName.c_str(), "RM" );
        triggerLargeAjdPhi[i][j]->Fit( triggerLargeName.c_str(), "RM" );
        recoilSmallAjdPhi[i][j]->Fit( recoilSmallName.c_str(), "RM" );
        recoilLargeAjdPhi[i][j]->Fit( recoilLargeName.c_str(), "RM" );
        
        // subtract
        TF1* subFunction = new TF1("subFunc", "[0]", phiMin, phiMax );
        subFunction->SetParameter( 0, triggerLargeTmpFit->GetParameter(0) );
        triggerLargeAjdPhi[i][j]->Add( subFunction, -1 );
        subFunction->SetParameter( 0, triggerSmallTmpFit->GetParameter(0) );
        triggerSmallAjdPhi[i][j]->Add( subFunction, -1 );
        subFunction->SetParameter( 0, recoilSmallTmpFit->GetParameter(0) );
        recoilSmallAjdPhi[i][j]->Add( subFunction, -1 );
        subFunction->SetParameter( 0, recoilLargeTmpFit->GetParameter(0) );
        recoilLargeAjdPhi[i][j]->Add( subFunction, -1 );
        
        // now do final fits
        triggerLargeName = "trigger_large_"; triggerLargeName += patch::to_string(i); triggerLargeName += patch::to_string(j);
        triggerSmallName = "trigger_small_"; triggerSmallName += patch::to_string(i); triggerSmallName += patch::to_string(j);
        recoilSmallName = "recoil_small_"; recoilSmallName += patch::to_string(i); recoilSmallName += patch::to_string(j);
        recoilLargeName = "recoil_large_"; recoilLargeName += patch::to_string(i); recoilLargeName += patch::to_string(j);
        
        triggerLargeAjdPhiFit[i][j] = new TF1(triggerLargeName.c_str(), phiForm.c_str(), phiMin, phiMax);
        triggerLargeAjdPhiFit[i][j]->FixParameter( 2, 0 );
        triggerLargeAjdPhiFit[i][j]->FixParameter( 5, corrAnalysis::pi );
        triggerLargeAjdPhiFit[i][j]->SetParameter( 3, 0.2 );
        triggerLargeAjdPhiFit[i][j]->SetParameter( 6, 0.2 );
        triggerSmallAjdPhiFit[i][j] = new TF1(triggerSmallName.c_str(), phiForm.c_str(), phiMin, phiMax);
        triggerSmallAjdPhiFit[i][j]->FixParameter( 2, 0 );
        triggerSmallAjdPhiFit[i][j]->FixParameter( 5, corrAnalysis::pi );
        triggerSmallAjdPhiFit[i][j]->SetParameter( 3, 0.2 );
        triggerSmallAjdPhiFit[i][j]->SetParameter( 6, 0.2 );
        recoilLargeAjdPhiFit[i][j] = new TF1(recoilLargeName.c_str(), phiForm.c_str(), phiMin, phiMax);
        recoilLargeAjdPhiFit[i][j]->FixParameter( 2, 0 );
        recoilLargeAjdPhiFit[i][j]->FixParameter( 5, corrAnalysis::pi );
        recoilLargeAjdPhiFit[i][j]->SetParameter( 3, 0.2 );
        recoilLargeAjdPhiFit[i][j]->SetParameter( 6, 0.2 );
        recoilSmallAjdPhiFit[i][j] = new TF1(recoilSmallName.c_str(), phiForm.c_str(), phiMin, phiMax);
        recoilSmallAjdPhiFit[i][j]->FixParameter( 2, 0 );
        recoilSmallAjdPhiFit[i][j]->FixParameter( 5, corrAnalysis::pi );
        recoilSmallAjdPhiFit[i][j]->SetParameter( 3, 0.2 );
        recoilSmallAjdPhiFit[i][j]->SetParameter( 6, 0.2 );
        
        // now scale the histograms
        triggerSmallAjdPhi[i][j]->Scale( 1.0/triggerSmallAjdPhi[i][j]->GetXaxis()->GetBinWidth(1) );
        triggerSmallAjdPhi[i][j]->Scale( 1.0 / (double) ajLowCount );
        triggerLargeAjdPhi[i][j]->Scale( 1.0/triggerLargeAjdPhi[i][j]->GetXaxis()->GetBinWidth(1) );
        triggerLargeAjdPhi[i][j]->Scale( 1.0 / (double) ajHighCount );
        recoilSmallAjdPhi[i][j]->Scale( 1.0/recoilSmallAjdPhi[i][j]->GetXaxis()->GetBinWidth(1) );
        recoilSmallAjdPhi[i][j]->Scale( 1.0 / (double) ajLowCount );
        recoilLargeAjdPhi[i][j]->Scale( 1.0/recoilLargeAjdPhi[i][j]->GetXaxis()->GetBinWidth(1) );
        recoilLargeAjdPhi[i][j]->Scale( 1.0 / (double) ajHighCount );
        
        // do the final fits
        triggerSmallAjdPhi[i][j]->Fit( triggerSmallName.c_str(), "RM" );
        triggerLargeAjdPhi[i][j]->Fit( triggerLargeName.c_str(), "RM" );
        recoilSmallAjdPhi[i][j]->Fit( recoilSmallName.c_str(), "RM" );
        recoilLargeAjdPhi[i][j]->Fit( recoilLargeName.c_str(), "RM" );
        
        // now copy and subtract
        std::string subTriggerName = "ajsub_dphi_trigger_file_"; subTriggerName += patch::to_string(i); subTriggerName += patch::to_string(j);
        std::string subRecoilName = "ajsub_dphi_recoil_file_"; subRecoilName += patch::to_string(i); subRecoilName += patch::to_string(j);
        
        
        triggerSubtracted[i][j] = (TH1D*) triggerSmallAjdPhi[i][j]->Clone();
        triggerSubtracted[i][j]->SetName( subTriggerName.c_str() );
        recoilSubtracted[i][j] = (TH1D*) recoilSmallAjdPhi[i][j]->Clone();
        recoilSubtracted[i][j]->SetName( subRecoilName.c_str() );
        
        triggerSubtracted[i][j]->Add( triggerLargeAjdPhi[i][j], -1 );
        recoilSubtracted[i][j]->Add( recoilLargeAjdPhi[i][j], -1 );
        
        subTriggerName = "trigger_sub_"; subTriggerName += patch::to_string(i); subTriggerName += patch::to_string(j);
        subRecoilName = "recoil_sub_"; subRecoilName += patch::to_string(i); subRecoilName += patch::to_string(j);
        
        triggerSubFit[i][j] = new TF1(subTriggerName.c_str(), phiForm.c_str(), phiMin, phiMax);
        triggerSubFit[i][j]->FixParameter( 2, 0 );
        triggerSubFit[i][j]->FixParameter( 5, corrAnalysis::pi );
        triggerSubFit[i][j]->SetParameter( 3, 0.2 );
        triggerSubFit[i][j]->SetParameter( 6, 0.2 );
        
        recoilSubFit[i][j] = new TF1(subRecoilName.c_str(), phiForm.c_str(), phiMin, phiMax);
        recoilSubFit[i][j]->FixParameter( 2, 0 );
        recoilSubFit[i][j]->FixParameter( 5, corrAnalysis::pi );
        recoilSubFit[i][j]->SetParameter( 3, 0.2 );
        recoilSubFit[i][j]->SetParameter( 6, 0.2 );
        
        triggerSubtracted[i][j]->Fit( subTriggerName.c_str(), "RM" );
        recoilSubtracted[i][j]->Fit( subRecoilName.c_str(), "RM" );
        
        triggerLargeYieldNear[i][j] = triggerLargeAjdPhiFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(triggerLargeAjdPhiFit[i][j]->GetParameter(3))/ptBinWidth[j];
        triggerLargeYieldNearErr[i][j] = triggerLargeAjdPhiFit[i][j]->GetParError(1);
        triggerLargeWidthNear[i][j] = triggerLargeAjdPhiFit[i][j]->GetParameter(3);
        triggerLargeWidthNearErr[i][j] = triggerLargeAjdPhiFit[i][j]->GetParError(3);
        
        triggerSmallYieldNear[i][j] = triggerSmallAjdPhiFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(triggerSmallAjdPhiFit[i][j]->GetParameter(3))/ptBinWidth[j];
        triggerSmallYieldNearErr[i][j] = triggerSmallAjdPhiFit[i][j]->GetParError(1);
        triggerSmallWidthNear[i][j] = triggerSmallAjdPhiFit[i][j]->GetParameter(3);
        triggerSmallWidthNearErr[i][j] = triggerSmallAjdPhiFit[i][j]->GetParError(3);
        
        recoilLargeYieldNear[i][j] = recoilLargeAjdPhiFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(recoilLargeAjdPhiFit[i][j]->GetParameter(3))/ptBinWidth[j];
        recoilLargeYieldNearErr[i][j] = recoilLargeAjdPhiFit[i][j]->GetParError(1);
        recoilLargeWidthNear[i][j] = recoilLargeAjdPhiFit[i][j]->GetParameter(3);
        recoilLargeWidthNearErr[i][j] = recoilLargeAjdPhiFit[i][j]->GetParError(3);
        
        recoilSmallYieldNear[i][j] = recoilSmallAjdPhiFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(recoilSmallAjdPhiFit[i][j]->GetParameter(3))/ptBinWidth[j];
        recoilSmallYieldNearErr[i][j] = recoilSmallAjdPhiFit[i][j]->GetParError(1);
        recoilSmallWidthNear[i][j] = recoilSmallAjdPhiFit[i][j]->GetParameter(3);
        recoilSmallWidthNearErr[i][j] = recoilSmallAjdPhiFit[i][j]->GetParError(3);
        
        triggerLargeYieldAway[i][j] = triggerLargeAjdPhiFit[i][j]->GetParameter(4)*sqrt(2*corrAnalysis::pi)*fabs(triggerLargeAjdPhiFit[i][j]->GetParameter(6))/ptBinWidth[j];
        triggerLargeYieldAwayErr[i][j] = triggerLargeAjdPhiFit[i][j]->GetParError(4);
        triggerLargeWidthAway[i][j] = triggerLargeAjdPhiFit[i][j]->GetParameter(6);
        triggerLargeWidthAwayErr[i][j] = triggerLargeAjdPhiFit[i][j]->GetParError(6);
        
        triggerSmallYieldAway[i][j] = triggerSmallAjdPhiFit[i][j]->GetParameter(4)*sqrt(2*corrAnalysis::pi)*fabs(triggerSmallAjdPhiFit[i][j]->GetParameter(6))/ptBinWidth[j];
        triggerSmallYieldAwayErr[i][j] = triggerSmallAjdPhiFit[i][j]->GetParError(4);
        triggerSmallWidthAway[i][j] = triggerSmallAjdPhiFit[i][j]->GetParameter(6);
        triggerSmallWidthAwayErr[i][j] = triggerSmallAjdPhiFit[i][j]->GetParError(6);
        
        recoilLargeYieldAway[i][j] = recoilLargeAjdPhiFit[i][j]->GetParameter(4)*sqrt(2*corrAnalysis::pi)*fabs(recoilLargeAjdPhiFit[i][j]->GetParameter(6))/ptBinWidth[j];
        recoilLargeYieldAwayErr[i][j] = recoilLargeAjdPhiFit[i][j]->GetParError(4);
        recoilLargeWidthAway[i][j] = recoilLargeAjdPhiFit[i][j]->GetParameter(6);
        recoilLargeWidthAwayErr[i][j] = recoilLargeAjdPhiFit[i][j]->GetParError(6);
        
        recoilSmallYieldAway[i][j] = recoilSmallAjdPhiFit[i][j]->GetParameter(4)*sqrt(2*corrAnalysis::pi)*fabs(recoilSmallAjdPhiFit[i][j]->GetParameter(6))/ptBinWidth[j];
        recoilSmallYieldAwayErr[i][j] = recoilSmallAjdPhiFit[i][j]->GetParError(4);
        recoilSmallWidthAway[i][j] = recoilSmallAjdPhiFit[i][j]->GetParameter(6);
        recoilSmallWidthAwayErr[i][j] = recoilSmallAjdPhiFit[i][j]->GetParError(6);
        
        triggerSubYieldNear[i][j] = triggerSubFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(triggerSubFit[i][j]->GetParameter(3))/ptBinWidth[j];
        triggerSubYieldNearErr[i][j] = triggerSubFit[i][j]->GetParError(1);
        triggerSubWidthNear[i][j] = triggerSubFit[i][j]->GetParameter(3);
        triggerSubWidthNearErr[i][j] = triggerSubFit[i][j]->GetParError(3);
        
        triggerSubYieldAway[i][j] = triggerSubFit[i][j]->GetParameter(4)*sqrt(2*corrAnalysis::pi)*fabs(triggerSubFit[i][j]->GetParameter(6))/ptBinWidth[j];
        triggerSubYieldAwayErr[i][j] = triggerSubFit[i][j]->GetParError(4);
        triggerSubWidthAway[i][j] = triggerSubFit[i][j]->GetParameter(6);
        triggerSubWidthAwayErr[i][j] = triggerSubFit[i][j]->GetParError(6);

        recoilSubYieldNear[i][j] = recoilSubFit[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(recoilSubFit[i][j]->GetParameter(3))/ptBinWidth[j];
        recoilSubYieldNearErr[i][j] = recoilSubFit[i][j]->GetParError(1);
        recoilSubWidthNear[i][j] = recoilSubFit[i][j]->GetParameter(3);
        recoilSubWidthNearErr[i][j] = recoilSubFit[i][j]->GetParError(3);
        
        recoilSubYieldAway[i][j] = recoilSubFit[i][j]->GetParameter(4)*sqrt(2*corrAnalysis::pi)*fabs(recoilSubFit[i][j]->GetParameter(6))/ptBinWidth[j];
        recoilSubYieldAwayErr[i][j] = recoilSubFit[i][j]->GetParError(4);
        recoilSubWidthAway[i][j] = recoilSubFit[i][j]->GetParameter(6);
        recoilSubWidthAwayErr[i][j] = recoilSubFit[i][j]->GetParError(6);
        
      }
    }
  }
 
  // now do output
  
  
  
  // ***************************

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
  std::vector<std::vector<TH1D*> > dPhiLeadLarge;
  std::vector<std::vector<TH1D*> > dPhiLeadNearLarge;
  std::vector<std::vector<TH1D*> > dPhiLeadFarLarge;
  std::vector<std::vector<TH1D*> > dEtaLeadLarge;
  std::vector<std::vector<TH1D*> > dPhiLeadSmall;
  std::vector<std::vector<TH1D*> > dPhiLeadNearSmall;
  std::vector<std::vector<TH1D*> > dPhiLeadFarSmall;
  std::vector<std::vector<TH1D*> > dEtaLeadSmall;
  
  std::vector<std::vector<TH1D*> > dPhiSubLarge;
  std::vector<std::vector<TH1D*> > dPhiSubNearLarge;
  std::vector<std::vector<TH1D*> > dPhiSubFarLarge;
  std::vector<std::vector<TH1D*> > dEtaSubLarge;
  std::vector<std::vector<TH1D*> > dPhiSubSmall;
  std::vector<std::vector<TH1D*> > dPhiSubNearSmall;
  std::vector<std::vector<TH1D*> > dPhiSubFarSmall;
  std::vector<std::vector<TH1D*> > dEtaSubSmall;
  
  dPhiLeadLarge.resize( nFiles );
  dPhiLeadNearLarge.resize( nFiles );
  dPhiLeadFarLarge.resize( nFiles );
  dEtaLeadLarge.resize( nFiles );
  dPhiLeadSmall.resize( nFiles );
  dPhiLeadNearSmall.resize( nFiles );
  dPhiLeadFarSmall.resize( nFiles );
  dEtaLeadSmall.resize( nFiles );
  
  dPhiSubLarge.resize( nFiles );
  dPhiSubNearLarge.resize( nFiles );
  dPhiSubFarLarge.resize( nFiles );
  dEtaSubLarge.resize( nFiles );
  dPhiSubSmall.resize( nFiles );
  dPhiSubNearSmall.resize( nFiles );
  dPhiSubFarSmall.resize( nFiles );
  dEtaSubSmall.resize( nFiles );
  
  for ( int i = 0; i < nFiles; ++i ) {

    dPhiLeadLarge[i].resize( nPtBins );
    dPhiLeadNearLarge[i].resize( nPtBins );
    dPhiLeadFarLarge[i].resize( nPtBins );
    dEtaLeadLarge[i].resize( nPtBins );
    dPhiLeadSmall[i].resize( nPtBins );
    dPhiLeadNearSmall[i].resize( nPtBins );
    dPhiLeadFarSmall[i].resize( nPtBins );
    dEtaLeadSmall[i].resize( nPtBins );
    
    dPhiSubLarge[i].resize( nPtBins );
    dPhiSubNearLarge[i].resize( nPtBins );
    dPhiSubFarLarge[i].resize( nPtBins );
    dEtaSubLarge[i].resize( nPtBins );
    dPhiSubSmall[i].resize( nPtBins );
    dPhiSubNearSmall[i].resize( nPtBins );
    dPhiSubFarSmall[i].resize( nPtBins );
    dEtaSubSmall[i].resize( nPtBins );
    
    for ( int j = 0; j < nPtBins; ++j ) {
      // first restrict the eta range
      recombinedCorrLarge[i][j]->GetXaxis()->SetRangeUser( etaMin, etaMax  );
      recombinedSubLarge[i][j]->GetXaxis()->SetRangeUser( etaMin, etaMax );
      
      if ( ajSplit[i] ) {
        recombinedCorrSmall[i][j]->GetXaxis()->SetRangeUser( etaMin, etaMax  );
        recombinedSubSmall[i][j]->GetXaxis()->SetRangeUser( etaMin, etaMax );
      }
      
      // save the 2D histograms
      std::string leadOutName = "tmp/lead2d_" + analysisNames[i] +"_pt_"+ patch::to_string(j) + ".pdf";
      std::string subOutName = "tmp/sub2d_" + analysisNames[i] +"_pt_" + patch::to_string(j) + ".pdf";
      std::string leadOutNameSmall = "";
      std::string subOutNameSmall = "";
      
      if ( ajSplit[i] ) {
        leadOutName = "tmp/lead2d_largeAj_" + analysisNames[i] +"_pt_"+ patch::to_string(j) + ".pdf";
        subOutName = "tmp/sub2d_largeAj_" + analysisNames[i] +"_pt_" + patch::to_string(j) + ".pdf";
        leadOutNameSmall = "tmp/lead2d_smallAj_" + analysisNames[i] +"_pt_"+ patch::to_string(j) + ".pdf";
        subOutNameSmall = "tmp/sub2d_smallAj_" + analysisNames[i] +"_pt_" + patch::to_string(j) + ".pdf";
      }
      
      TCanvas c1;
      recombinedCorrLarge[i][j]->Draw("surf1");
      c1.SaveAs( leadOutName.c_str() );
      recombinedSubLarge[i][j]->Draw("surf1");
      c1.SaveAs( subOutName.c_str() );
      
      if ( ajSplit[i] ) {
        recombinedCorrSmall[i][j]->Draw("surf1");
        c1.SaveAs( leadOutNameSmall.c_str() );
        recombinedSubSmall[i][j]->Draw("surf1");
        c1.SaveAs( subOutNameSmall.c_str() );
      }
      
      dPhiLeadLarge[i][j] = (TH1D*) ((TH1D*) recombinedCorrLarge[i][j]->ProjectionY())->Clone();
      dPhiSubLarge[i][j] = (TH1D*) ((TH1D*) recombinedSubLarge[i][j]->ProjectionY())->Clone();
      
      if ( ajSplit[i] ) {
        dPhiLeadSmall[i][j] = (TH1D*) ((TH1D*) recombinedCorrSmall[i][j]->ProjectionY())->Clone();
        dPhiSubSmall[i][j] = (TH1D*) ((TH1D*) recombinedSubSmall[i][j]->ProjectionY())->Clone();
      }
      
      recombinedCorrLarge[i][j]->GetYaxis()->SetRangeUser( phiMinClose, phiMaxClose );
      recombinedSubLarge[i][j]->GetYaxis()->SetRangeUser( phiMinClose, phiMaxClose );
      if ( ajSplit[i] ) {
        recombinedCorrSmall[i][j]->GetYaxis()->SetRangeUser( phiMinClose, phiMaxClose );
        recombinedSubSmall[i][j]->GetYaxis()->SetRangeUser( phiMinClose, phiMaxClose );
      }

      dEtaLeadLarge[i][j] = (TH1D*) ((TH1D*) recombinedCorrLarge[i][j]->ProjectionX())->Clone();
      dEtaSubLarge[i][j] = (TH1D*) ((TH1D*) recombinedSubLarge[i][j]->ProjectionX())->Clone();
      if ( ajSplit[i] ) {
        dEtaLeadSmall[i][j] = (TH1D*) ((TH1D*) recombinedCorrSmall[i][j]->ProjectionX())->Clone();
        dEtaSubSmall[i][j] = (TH1D*) ((TH1D*) recombinedSubSmall[i][j]->ProjectionX())->Clone();
      }
      
      recombinedCorrLarge[i][j]->GetYaxis()->SetRangeUser( phiMin, phiMaxFar );
      recombinedSubLarge[i][j]->GetYaxis()->SetRangeUser( phiMin, phiMaxFar );
      if ( ajSplit[i] ) {
        recombinedCorrSmall[i][j]->GetYaxis()->SetRangeUser( phiMin, phiMaxFar );
        recombinedSubSmall[i][j]->GetYaxis()->SetRangeUser( phiMin, phiMaxFar );
      }

      // now get dphi in "near" and "far" eta ranges
      recombinedCorrLarge[i][j]->GetXaxis()->SetRange( etaNearMinBin, etaNearMaxBin  );
      recombinedSubLarge[i][j]->GetXaxis()->SetRange( etaNearMinBin, etaNearMaxBin );
      if ( ajSplit[i] ) {
        recombinedCorrSmall[i][j]->GetXaxis()->SetRange( etaNearMinBin, etaNearMaxBin  );
        recombinedSubSmall[i][j]->GetXaxis()->SetRange( etaNearMinBin, etaNearMaxBin );
      }
      
      dPhiLeadNearLarge[i][j] = (TH1D*) ((TH1D*) recombinedCorrLarge[i][j]->ProjectionY())->Clone();
      dPhiSubNearLarge[i][j] = (TH1D*) ((TH1D*) recombinedSubLarge[i][j]->ProjectionY())->Clone();
      if ( ajSplit[i] ) {
        dPhiLeadNearSmall[i][j] = (TH1D*) ((TH1D*) recombinedCorrSmall[i][j]->ProjectionY())->Clone();
        dPhiSubNearSmall[i][j] = (TH1D*) ((TH1D*) recombinedSubSmall[i][j]->ProjectionY())->Clone();
      }

      recombinedCorrLarge[i][j]->GetXaxis()->SetRange( etaFarMinBin, etaFarMaxBin  );
      recombinedSubLarge[i][j]->GetXaxis()->SetRange( etaFarMinBin, etaFarMaxBin );
      if ( ajSplit[i] ) {
        recombinedCorrSmall[i][j]->GetXaxis()->SetRange( etaFarMinBin, etaFarMaxBin  );
        recombinedSubSmall[i][j]->GetXaxis()->SetRange( etaFarMinBin, etaFarMaxBin );
      }
      
      dPhiLeadFarLarge[i][j] = (TH1D*) ((TH1D*) recombinedCorrLarge[i][j]->ProjectionY())->Clone();
      dPhiSubFarLarge[i][j] = (TH1D*) ((TH1D*) recombinedSubLarge[i][j]->ProjectionY())->Clone();
      if ( ajSplit[i] ) {
        dPhiLeadFarSmall[i][j] = (TH1D*) ((TH1D*) recombinedCorrSmall[i][j]->ProjectionY())->Clone();
        dPhiSubFarSmall[i][j] = (TH1D*) ((TH1D*) recombinedSubSmall[i][j]->ProjectionY())->Clone();
        
      }

      recombinedCorrLarge[i][j]->GetXaxis()->SetRange( etaMinBin, etaMaxBin  );
      recombinedSubLarge[i][j]->GetXaxis()->SetRange( etaMinBin, etaMaxBin );
      if ( ajSplit[i] ) {
        recombinedCorrSmall[i][j]->GetXaxis()->SetRange( etaMinBin, etaMaxBin  );
        recombinedSubSmall[i][j]->GetXaxis()->SetRange( etaMinBin, etaMaxBin );
      }
      
      dPhiLeadFarLarge[i][j]->Add( recombinedCorrLarge[i][j]->ProjectionY() );
      dPhiSubFarLarge[i][j]->Add( recombinedSubLarge[i][j]->ProjectionY() );
      if ( ajSplit[i] ) {
        dPhiLeadFarSmall[i][j]->Add( recombinedCorrSmall[i][j]->ProjectionY() );
        dPhiSubFarSmall[i][j]->Add( recombinedSubSmall[i][j]->ProjectionY() );
      }

    }
  }
  

  // now  first overlay and output,
  // then subtract near from far eta regions
  for ( int i = 0; i < nFiles; ++i )
    for ( int j = 0; j < nPtBins; ++j ) {
      std::string leadPhiName = "tmp/lead_phi_near_far_"+analysisNames[i]+"_pt_"+patch::to_string(j)+".pdf";
      std::string subPhiName = "tmp/sub_phi_near_far_"+analysisNames[i]+"_pt_"+patch::to_string(j)+".pdf";
      
      std::string leadPhiNameSmall = "";
      std::string subPhiNameSmall = "";
      
      if ( ajSplit[i] ) {
        leadPhiName = "tmp/largeAj_lead_phi_near_far_"+analysisNames[i]+"_pt_"+patch::to_string(j)+".pdf";
        subPhiName = "tmp/largeAj_sub_phi_near_far_"+analysisNames[i]+"_pt_"+patch::to_string(j)+".pdf";
        leadPhiNameSmall = "tmp/smallAj_lead_phi_near_far_"+analysisNames[i]+"_pt_"+patch::to_string(j)+".pdf";
        subPhiNameSmall = "tmp/smallAj_sub_phi_near_far_"+analysisNames[i]+"_pt_"+patch::to_string(j)+".pdf";
      }
      
      
      TCanvas c1;
      dPhiLeadNearLarge[i][j]->SetLineColor(kBlack);
      dPhiLeadNearLarge[i][j]->SetMarkerStyle(29);
      dPhiLeadNearLarge[i][j]->SetMarkerSize(3);
      dPhiLeadNearLarge[i][j]->SetMarkerColor(kBlack);
      dPhiLeadNearLarge[i][j]->Draw();
      dPhiLeadFarLarge[i][j]->SetLineColor(kRed);
      dPhiLeadFarLarge[i][j]->SetMarkerStyle(29);
      dPhiLeadFarLarge[i][j]->SetMarkerSize(3);
      dPhiLeadFarLarge[i][j]->SetMarkerColor(kRed);
      dPhiLeadFarLarge[i][j]->Draw("SAME");
      c1.SaveAs( leadPhiName.c_str() );
      dPhiSubNearLarge[i][j]->SetLineColor(kBlack);
      dPhiSubNearLarge[i][j]->SetMarkerStyle(29);
      dPhiSubNearLarge[i][j]->SetMarkerSize(3);
      dPhiSubNearLarge[i][j]->SetMarkerColor(kBlack);
      dPhiSubNearLarge[i][j]->Draw();
      dPhiSubFarLarge[i][j]->SetLineColor(kRed);
      dPhiSubFarLarge[i][j]->SetMarkerStyle(29);
      dPhiSubFarLarge[i][j]->SetMarkerSize(3);
      dPhiSubFarLarge[i][j]->SetMarkerColor(kRed);
      dPhiSubFarLarge[i][j]->Draw("SAME");
      c1.SaveAs( subPhiName.c_str() );
      
      if ( ajSplit[i] ) {
        dPhiLeadNearSmall[i][j]->SetLineColor(kBlack);
        dPhiLeadNearSmall[i][j]->SetMarkerStyle(29);
        dPhiLeadNearSmall[i][j]->SetMarkerSize(3);
        dPhiLeadNearSmall[i][j]->SetMarkerColor(kBlack);
        dPhiLeadNearSmall[i][j]->Draw();
        dPhiLeadFarSmall[i][j]->SetLineColor(kRed);
        dPhiLeadFarSmall[i][j]->SetMarkerStyle(29);
        dPhiLeadFarSmall[i][j]->SetMarkerSize(3);
        dPhiLeadFarSmall[i][j]->SetMarkerColor(kRed);
        dPhiLeadFarSmall[i][j]->Draw("SAME");
        c1.SaveAs( leadPhiNameSmall.c_str() );
        dPhiSubNearSmall[i][j]->SetLineColor(kBlack);
        dPhiSubNearSmall[i][j]->SetMarkerStyle(29);
        dPhiSubNearSmall[i][j]->SetMarkerSize(3);
        dPhiSubNearSmall[i][j]->SetMarkerColor(kBlack);
        dPhiSubNearSmall[i][j]->Draw();
        dPhiSubFarSmall[i][j]->SetLineColor(kRed);
        dPhiSubFarSmall[i][j]->SetMarkerStyle(29);
        dPhiSubFarSmall[i][j]->SetMarkerSize(3);
        dPhiSubFarSmall[i][j]->SetMarkerColor(kRed);
        dPhiSubFarSmall[i][j]->Draw("SAME");
        c1.SaveAs( subPhiNameSmall.c_str() );

      }

      
      // Now do the subtraction
      dPhiLeadNearLarge[i][j]->Add( dPhiLeadFarLarge[i][j], -1 );
      dPhiSubNearLarge[i][j]->Add( dPhiSubFarLarge[i][j], -1 );
      
      if ( ajSplit[i] ) {
        dPhiLeadNearSmall[i][j]->Add( dPhiLeadFarSmall[i][j], -1 );
        dPhiSubNearSmall[i][j]->Add( dPhiSubFarSmall[i][j], -1 );

      }
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
      
      std::string dPhiLeadNameSmall, dPhiSubNameSmall, dPhiLeadNameDifSmall, dPhiSubNameDifSmall, dEtaLeadNameSmall, dEtaSubNameSmall;

      dPhiLeadNameSmall = dPhiLeadName +"_small";
      dPhiSubNameSmall = dPhiSubName + "_small";
      dPhiLeadNameDifSmall = dPhiLeadNameDif + "_small";
      dPhiSubNameDifSmall = dPhiSubNameDif + "_small";
      dEtaLeadNameSmall = dEtaLeadName + "_small";
      dEtaSubNameSmall = dEtaSubName + "_small";
      
      TF1* leadPhiInitFit = new TF1( dPhiLeadName.c_str(), phiForm.c_str(), phiMin, phiMax );
      leadPhiInitFit->FixParameter( 2, 0 );
      leadPhiInitFit->FixParameter( 5, corrAnalysis::pi );
      leadPhiInitFit->SetParameter( 3, 0.2 );
      leadPhiInitFit->SetParameter( 6, 0.2 );
      
      TF1* leadPhiSmallInitFit = new TF1( dPhiLeadNameSmall.c_str(), phiForm.c_str(), phiMin, phiMax );
      leadPhiSmallInitFit->FixParameter( 2, 0 );
      leadPhiSmallInitFit->FixParameter( 5, corrAnalysis::pi );
      leadPhiSmallInitFit->SetParameter( 3, 0.2 );
      leadPhiSmallInitFit->SetParameter( 6, 0.2 );
      
      TF1* leadPhiDifInitFit = new TF1( dPhiLeadNameDif.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      leadPhiDifInitFit->FixParameter( 2, 0 );
      leadPhiDifInitFit->SetParameter( 3, 0.2 );
      
      TF1* leadPhiSmallDifInitFit = new TF1( dPhiLeadNameDifSmall.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      leadPhiSmallDifInitFit->FixParameter( 2, 0 );
      leadPhiSmallDifInitFit->SetParameter( 3, 0.2 );
      
      TF1* subPhiInitFit = new TF1( dPhiSubName.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiInitFit->FixParameter( 2, 0 );
      subPhiInitFit->FixParameter( 5, corrAnalysis::pi );
      subPhiInitFit->SetParameter( 3, 0.2 );
      subPhiInitFit->SetParameter( 6, 0.2 );
      
      TF1* subPhiSmallInitFit = new TF1( dPhiSubNameSmall.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiSmallInitFit->FixParameter( 2, 0 );
      subPhiSmallInitFit->FixParameter( 5, corrAnalysis::pi );
      subPhiSmallInitFit->SetParameter( 3, 0.2 );
      subPhiSmallInitFit->SetParameter( 6, 0.2 );
      
      TF1* subPhiDifInitFit = new TF1( dPhiSubNameDif.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      subPhiDifInitFit->FixParameter( 2, 0 );
      subPhiDifInitFit->SetParameter( 3, 0.2 );
      
      TF1* subPhiDifSmallInitFit = new TF1( dPhiSubNameDifSmall.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      subPhiDifSmallInitFit->FixParameter( 2, 0 );
      subPhiDifSmallInitFit->SetParameter( 3, 0.2 );
      
      TF1* leadEtaInitFit = new TF1( dEtaLeadName.c_str(), etaForm.c_str(), etaMin, etaMax );
      leadEtaInitFit->FixParameter( 2, 0 );
      leadEtaInitFit->SetParameter( 3, 0.2 );
      
      TF1* leadEtaSmallInitFit = new TF1( dEtaLeadNameSmall.c_str(), etaForm.c_str(), etaMin, etaMax );
      leadEtaSmallInitFit->FixParameter( 2, 0 );
      leadEtaSmallInitFit->SetParameter( 3, 0.2 );
      
      TF1* subEtaInitFit = new TF1( dEtaSubName.c_str(), etaForm.c_str(), etaMin, etaMax );
      subEtaInitFit->FixParameter( 2, 0 );
      subEtaInitFit->SetParameter( 3, 0.2 );
      
      TF1* subEtaSmallInitFit = new TF1( dEtaSubNameSmall.c_str(), etaForm.c_str(), etaMin, etaMax );
      subEtaSmallInitFit->FixParameter( 2, 0 );
      subEtaSmallInitFit->SetParameter( 3, 0.2 );
      
      dPhiLeadLarge[i][j]->Fit( dPhiLeadName.c_str(), "RM" );
      dPhiSubLarge[i][j]->Fit( dPhiSubName.c_str(), "RM" );
      dPhiLeadNearLarge[i][j]->Fit( dPhiLeadNameDif.c_str(), "RM" );
      dPhiSubNearLarge[i][j]->Fit( dPhiSubNameDif.c_str(), "RM" );
      dEtaLeadLarge[i][j]->Fit( dEtaLeadName.c_str(), "RM" );
      dEtaSubLarge[i][j]->Fit( dEtaSubName.c_str(), "RM" );
      
      if ( ajSplit[i] ) {
        dPhiLeadSmall[i][j]->Fit( dPhiLeadNameSmall.c_str(), "RM" );
        dPhiSubSmall[i][j]->Fit( dPhiSubNameSmall.c_str(), "RM" );
        dPhiLeadNearSmall[i][j]->Fit( dPhiLeadNameDifSmall.c_str(), "RM" );
        dPhiSubNearSmall[i][j]->Fit( dPhiSubNameDifSmall.c_str(), "RM" );
        dEtaLeadSmall[i][j]->Fit( dEtaLeadNameSmall.c_str(), "RM" );
        dEtaSubSmall[i][j]->Fit( dEtaSubNameSmall.c_str(), "RM" );
      }
      
      // Now to subtract the constants
      TF1* subConst = new TF1( "subConst", "[0]", phiMin, phiMax);
      TF1* subConstEta = new TF1("subConstEta", "[0]", etaMin, etaMax);
      subConst->SetParameter( 0, leadPhiInitFit->GetParameter(0) );
      dPhiLeadLarge[i][j]->Add( subConst, -1 );
      subConst->SetParameter( 0, leadPhiDifInitFit->GetParameter(0));
      dPhiLeadNearLarge[i][j]->Add( subConst, -1 );
      subConstEta->SetParameter( 0, leadEtaInitFit->GetParameter(0));
      dEtaLeadLarge[i][j]->Add( subConstEta, -1 );
      
      subConst->SetParameter( 0, subPhiInitFit->GetParameter(0) );
      dPhiSubLarge[i][j]->Add( subConst, -1 );
      subConst->SetParameter( 0, subPhiDifInitFit->GetParameter(0));
      dPhiSubNearLarge[i][j]->Add( subConst, -1 );
      subConstEta->SetParameter( 0, subEtaInitFit->GetParameter(0));
      dEtaSubLarge[i][j]->Add( subConstEta, -1 );
      
      if ( ajSplit[i] ) {
        subConst->SetParameter( 0, leadPhiSmallInitFit->GetParameter(0) );
        dPhiLeadSmall[i][j]->Add( subConst, -1 );
        subConst->SetParameter( 0, leadPhiSmallDifInitFit->GetParameter(0));
        dPhiLeadNearSmall[i][j]->Add( subConst, -1 );
        subConstEta->SetParameter( 0, leadEtaSmallInitFit->GetParameter(0));
        dEtaLeadSmall[i][j]->Add( subConstEta, -1 );
        
        subConst->SetParameter( 0, subPhiSmallInitFit->GetParameter(0) );
        dPhiSubSmall[i][j]->Add( subConst, -1 );
        subConst->SetParameter( 0, subPhiDifSmallInitFit->GetParameter(0));
        dPhiSubNearSmall[i][j]->Add( subConst, -1 );
        subConstEta->SetParameter( 0, subEtaSmallInitFit->GetParameter(0));
        dEtaSubSmall[i][j]->Add( subConstEta, -1 );
      }
      
      // now scale the histograms
      if ( !ajSplit[i] ) {
        dPhiLeadLarge[i][j]->Scale( 1.0 / dPhiLeadLarge[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiLeadLarge[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
        dPhiLeadNearLarge[i][j]->Scale( 1.0 / dPhiLeadNearLarge[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiLeadNearLarge[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
        dEtaLeadLarge[i][j]->Scale( 1.0 / dEtaLeadLarge[i][j]->GetXaxis()->GetBinWidth(1) );
        dEtaLeadLarge[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
        dPhiSubLarge[i][j]->Scale( 1.0 / dPhiSubLarge[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiSubLarge[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
        dPhiSubNearLarge[i][j]->Scale( 1.0 / dPhiSubNearLarge[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiSubNearLarge[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
        dEtaSubLarge[i][j]->Scale( 1.0 / dEtaSubLarge[i][j]->GetXaxis()->GetBinWidth(1) );
        dEtaSubLarge[i][j]->Scale( 1.0 / (double) nEvents[i]->GetEntries() );
      }
      else {
        
        // we have to find the integral depending on small or large aj
        double ajHighCount = nEvents[i]->Integral( 1, 1, 1, corrAnalysis::binsCentrality, 1, corrAnalysis::binsVz );
        double ajLowCount = nEvents[i]->Integral( 2, 2, 1, corrAnalysis::binsCentrality, 1, corrAnalysis::binsVz );
        
        dPhiLeadLarge[i][j]->Scale( 1.0 / dPhiLeadSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiLeadLarge[i][j]->Scale( 1.0 / ajHighCount );
        dPhiLeadNearLarge[i][j]->Scale( 1.0 / dPhiLeadNearSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiLeadNearLarge[i][j]->Scale( 1.0 / ajHighCount );
        dEtaLeadLarge[i][j]->Scale( 1.0 / dEtaLeadSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dEtaLeadLarge[i][j]->Scale( 1.0 / ajHighCount );
        dPhiSubLarge[i][j]->Scale( 1.0 / dPhiSubSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiSubLarge[i][j]->Scale( 1.0 / ajHighCount );
        dPhiSubNearLarge[i][j]->Scale( 1.0 / dPhiSubNearSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiSubNearLarge[i][j]->Scale( 1.0 / ajHighCount );
        dEtaSubLarge[i][j]->Scale( 1.0 / dEtaSubSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dEtaSubLarge[i][j]->Scale( 1.0 / ajHighCount );
        
        dPhiLeadSmall[i][j]->Scale( 1.0 / dPhiLeadSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiLeadSmall[i][j]->Scale( 1.0 / ajLowCount );
        dPhiLeadNearSmall[i][j]->Scale( 1.0 / dPhiLeadNearSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiLeadNearSmall[i][j]->Scale( 1.0 / ajLowCount );
        dEtaLeadSmall[i][j]->Scale( 1.0 / dEtaLeadSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dEtaLeadSmall[i][j]->Scale( 1.0 / ajLowCount );
        dPhiSubSmall[i][j]->Scale( 1.0 / dPhiSubSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiSubSmall[i][j]->Scale( 1.0 / ajLowCount );
        dPhiSubNearSmall[i][j]->Scale( 1.0 / dPhiSubNearSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dPhiSubNearSmall[i][j]->Scale( 1.0 / ajLowCount );
        dEtaSubSmall[i][j]->Scale( 1.0 / dEtaSubSmall[i][j]->GetXaxis()->GetBinWidth(1) );
        dEtaSubSmall[i][j]->Scale( 1.0 / ajLowCount );
      }
    }
  }

  // final fitting
  std::vector<std::vector<TF1*> > leadPhiFitLarge;
  leadPhiFitLarge.resize( nFiles );
  std::vector<std::vector<TF1*> > leadPhiFitSmall;
  leadPhiFitSmall.resize( nFiles );
  std::vector<std::vector<TF1*> > leadPhiDifFitLarge;
  leadPhiDifFitLarge.resize( nFiles );
  std::vector<std::vector<TF1*> > leadPhiDifFitSmall;
  leadPhiDifFitSmall.resize( nFiles );
  std::vector<std::vector<TF1*> > leadEtaFitLarge;
  leadEtaFitLarge.resize( nFiles );
  std::vector<std::vector<TF1*> > leadEtaFitSmall;
  leadEtaFitSmall.resize( nFiles );
  std::vector<std::vector<TF1*> > subPhiFitLarge;
  subPhiFitLarge.resize( nFiles );
  std::vector<std::vector<TF1*> > subPhiFitSmall;
  subPhiFitSmall.resize( nFiles );
  std::vector<std::vector<TF1*> > subPhiDifFitLarge;
  subPhiDifFitLarge.resize( nFiles );
  std::vector<std::vector<TF1*> > subPhiDifFitSmall;
  subPhiDifFitSmall.resize( nFiles );
  std::vector<std::vector<TF1*> > subEtaFitLarge;
  subEtaFitLarge.resize( nFiles );
  std::vector<std::vector<TF1*> > subEtaFitSmall;
  subEtaFitSmall.resize( nFiles );
  
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiFitLarge[i].resize( nPtBins );
    leadPhiDifFitLarge[i].resize( nPtBins );
    leadEtaFitLarge[i].resize( nPtBins );
    subPhiFitLarge[i].resize( nPtBins );
    subPhiDifFitLarge[i].resize( nPtBins );
    subEtaFitLarge[i].resize( nPtBins );

    leadPhiFitSmall[i].resize( nPtBins );
    leadPhiDifFitSmall[i].resize( nPtBins );
    leadEtaFitSmall[i].resize( nPtBins );
    subPhiFitSmall[i].resize( nPtBins );
    subPhiDifFitSmall[i].resize( nPtBins );
    subEtaFitSmall[i].resize( nPtBins );

    
    for ( int j = 0; j < nPtBins; ++j ) {
      std::string dPhiLeadName = "fit_lead_phi_" + patch::to_string(i) + patch::to_string(j);
       std::string dPhiLeadNameSmall = "small_fit_lead_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubName = "fit_sub_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubNameSmall = "small_fit_sub_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiLeadNameDif = "fit_lead_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiLeadNameDifSmall = "small_fit_lead_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubNameDif = "fit_sub_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubNameDifSmall = "small_fit_sub_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaLeadName = "fit_lead_eta_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaLeadNameSmall = "small_fit_lead_eta_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaSubName = "fit_sub_eta_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaSubNameSmall = "small_fit_sub_eta_" + patch::to_string(i) + patch::to_string(j);
      
      leadPhiFitLarge[i][j] = new TF1( dPhiLeadName.c_str(), phiForm.c_str(), phiMin, phiMax );
      leadPhiFitLarge[i][j]->FixParameter( 2, 0 );
      leadPhiFitLarge[i][j]->FixParameter( 5, corrAnalysis::pi );
      leadPhiFitLarge[i][j]->SetParameter( 3, 0.2 );
      leadPhiFitLarge[i][j]->SetParameter( 6, 0.2 );
      leadPhiFitLarge[i][j]->SetLineColor( i + 1 );
      
      leadPhiFitSmall[i][j] = new TF1( dPhiLeadNameSmall.c_str(), phiForm.c_str(), phiMin, phiMax );
      leadPhiFitSmall[i][j]->FixParameter( 2, 0 );
      leadPhiFitSmall[i][j]->FixParameter( 5, corrAnalysis::pi );
      leadPhiFitSmall[i][j]->SetParameter( 3, 0.2 );
      leadPhiFitSmall[i][j]->SetParameter( 6, 0.2 );
      leadPhiFitSmall[i][j]->SetLineColor( i + 1 + nFiles );
      
      leadPhiDifFitLarge[i][j] = new TF1( dPhiLeadNameDif.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      leadPhiDifFitLarge[i][j]->FixParameter( 2, 0 );
      leadPhiDifFitLarge[i][j]->SetParameter( 3, 0.2 );
      leadPhiDifFitLarge[i][j]->SetLineColor( i + 1 );
      
      leadPhiDifFitSmall[i][j] = new TF1( dPhiLeadNameDifSmall.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      leadPhiDifFitSmall[i][j]->FixParameter( 2, 0 );
      leadPhiDifFitSmall[i][j]->SetParameter( 3, 0.2 );
      leadPhiDifFitSmall[i][j]->SetLineColor( i + 1 + nFiles );
      
      subPhiFitLarge[i][j] = new TF1( dPhiSubName.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiFitLarge[i][j]->FixParameter( 2, 0 );
      subPhiFitLarge[i][j]->FixParameter( 5, corrAnalysis::pi );
      subPhiFitLarge[i][j]->SetParameter( 3, 0.2 );
      subPhiFitLarge[i][j]->SetParameter( 6, 0.2 );
      subPhiFitLarge[i][j]->SetLineColor( i + 1 );
      
      subPhiFitSmall[i][j] = new TF1( dPhiSubNameSmall.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiFitSmall[i][j]->FixParameter( 2, 0 );
      subPhiFitSmall[i][j]->FixParameter( 5, corrAnalysis::pi );
      subPhiFitSmall[i][j]->SetParameter( 3, 0.2 );
      subPhiFitSmall[i][j]->SetParameter( 6, 0.2 );
      subPhiFitSmall[i][j]->SetLineColor( i + 1 + nFiles );
      
      subPhiDifFitLarge[i][j] = new TF1( dPhiSubNameDif.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      subPhiDifFitLarge[i][j]->FixParameter( 2, 0 );
      subPhiDifFitLarge[i][j]->SetParameter( 3, 0.2 );
      subPhiDifFitLarge[i][j]->SetLineColor( i + 1 );

      subPhiDifFitSmall[i][j] = new TF1( dPhiSubNameDifSmall.c_str(), phiDifForm.c_str(), phiMin, phiDifMax );
      subPhiDifFitSmall[i][j]->FixParameter( 2, 0 );
      subPhiDifFitSmall[i][j]->SetParameter( 3, 0.2 );
      subPhiDifFitSmall[i][j]->SetLineColor( i + 1 + nFiles );

      
      leadEtaFitLarge[i][j] = new TF1( dEtaLeadName.c_str(), etaForm.c_str(), etaMin, etaMax );
      leadEtaFitLarge[i][j]->FixParameter( 2, 0 );
      leadEtaFitLarge[i][j]->SetParameter( 3, 0.2 );
      leadEtaFitLarge[i][j]->SetLineColor( i + 1 );
      
      leadEtaFitSmall[i][j] = new TF1( dEtaLeadNameSmall.c_str(), etaForm.c_str(), etaMin, etaMax );
      leadEtaFitSmall[i][j]->FixParameter( 2, 0 );
      leadEtaFitSmall[i][j]->SetParameter( 3, 0.2 );
      leadEtaFitSmall[i][j]->SetLineColor( i + 1 + nFiles );
      
      subEtaFitLarge[i][j] = new TF1( dEtaSubName.c_str(), etaForm.c_str(), etaMin, etaMax );
      subEtaFitLarge[i][j]->FixParameter( 2, 0 );
      subEtaFitLarge[i][j]->SetParameter( 3, 0.2 );
      subEtaFitLarge[i][j]->SetLineColor( i + 1 );
      
      subEtaFitSmall[i][j] = new TF1( dEtaSubNameSmall.c_str(), etaForm.c_str(), etaMin, etaMax );
      subEtaFitSmall[i][j]->FixParameter( 2, 0 );
      subEtaFitSmall[i][j]->SetParameter( 3, 0.2 );
      subEtaFitSmall[i][j]->SetLineColor( i + 1 + nFiles );
      
      // Now set same colors and fit
      dPhiLeadLarge[i][j]->SetLineColor( i + 1 );
      dPhiLeadLarge[i][j]->Fit( dPhiLeadName.c_str(), "RM" );
      dPhiSubLarge[i][j]->SetLineColor( i + 1  );
      dPhiSubLarge[i][j]->Fit( dPhiSubName.c_str(), "RM" );
      dPhiLeadNearLarge[i][j]->SetLineColor( i + 1 );
      dPhiLeadNearLarge[i][j]->Fit( dPhiLeadNameDif.c_str(), "RM" );
      dPhiSubNearLarge[i][j]->SetLineColor( i + 1 );
      dPhiSubNearLarge[i][j]->Fit( dPhiSubNameDif.c_str(), "RM" );
      dEtaLeadLarge[i][j]->SetLineColor( i + 1 );
      dEtaLeadLarge[i][j]->Fit( dEtaLeadName.c_str(), "RM" );
      dEtaSubLarge[i][j]->SetLineColor( i + 1 );
      dEtaSubLarge[i][j]->Fit( dEtaSubName.c_str(), "RM" );
      
      if ( ajSplit[i] ) {
        dPhiLeadSmall[i][j]->SetLineColor( i + 1 + nFiles );
        dPhiLeadSmall[i][j]->Fit( dPhiLeadNameSmall.c_str(), "RM" );
        dPhiSubSmall[i][j]->SetLineColor( i + 1 + nFiles);
        dPhiSubSmall[i][j]->Fit( dPhiSubNameSmall.c_str(), "RM" );
        dPhiLeadNearSmall[i][j]->SetLineColor( i + 1 + nFiles);
        dPhiLeadNearSmall[i][j]->Fit( dPhiLeadNameDifSmall.c_str(), "RM" );
        dPhiSubNearSmall[i][j]->SetLineColor( i + 1 + nFiles);
        dPhiSubNearSmall[i][j]->Fit( dPhiSubNameDifSmall.c_str(), "RM" );
        dEtaLeadSmall[i][j]->SetLineColor( i + 1 + nFiles);
        dEtaLeadSmall[i][j]->Fit( dEtaLeadNameSmall.c_str(), "RM" );
        dEtaSubSmall[i][j]->SetLineColor( i + 1 + nFiles);
        dEtaSubSmall[i][j]->Fit( dEtaSubNameSmall.c_str(), "RM" );
      }
      
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
      dPhiLeadLarge[j][i]->SetLineColor(j+1);
      dPhiLeadLarge[j][i]->SetMarkerStyle(29);
      dPhiLeadLarge[j][i]->SetMarkerSize(3);
      dPhiLeadLarge[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        std::string outTitle = "Trigger Jet #Delta#phi " + ptBinString[i];
        dPhiLeadLarge[j][i]->SetTitle( outTitle.c_str() );
        dPhiLeadLarge[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiLeadLarge[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiLeadLarge[j][i]->Draw();
      }
      else {
        dPhiLeadLarge[j][i]->Draw("same");
      }
      if ( ajSplit[j] ) {
        dPhiLeadSmall[j][i]->SetLineColor(j+1+nFiles);
        dPhiLeadSmall[j][i]->SetMarkerStyle(29);
        dPhiLeadSmall[j][i]->SetMarkerSize(3);
        dPhiLeadSmall[j][i]->SetMarkerColor(j+1+nFiles);
        dPhiLeadSmall[j][i]->Draw("same");
        leadPhiFitSmall[j][i]->Draw("same");
      }
    }
    c1.SaveAs( leadPhiOut.c_str() );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string leadPhiDifOut = leadPhiDifOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      dPhiLeadNearLarge[j][i]->SetLineColor(j+1);
      dPhiLeadNearLarge[j][i]->SetMarkerStyle(29);
      dPhiLeadNearLarge[j][i]->SetMarkerSize(3);
      dPhiLeadNearLarge[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        std::string outTitle = "Trigger Jet #Delta#eta subtracted #Delta#phi " + ptBinString[i];
        dPhiLeadNearLarge[j][i]->SetTitle( outTitle.c_str() );
        dPhiLeadNearLarge[j][i]->GetXaxis()->SetRangeUser(-corrAnalysis::pi/2.0, corrAnalysis::pi/2.0);
        dPhiLeadNearLarge[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiLeadNearLarge[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiLeadNearLarge[j][i]->Draw();
      }
      else {
        dPhiLeadNearLarge[j][i]->Draw("same");
      }
      if ( ajSplit[j] ) {
        dPhiLeadNearSmall[j][i]->SetLineColor(j+1+nFiles);
        dPhiLeadNearSmall[j][i]->SetMarkerStyle(29);
        dPhiLeadNearSmall[j][i]->SetMarkerSize(3);
        dPhiLeadNearSmall[j][i]->SetMarkerColor(j+1+nFiles);
        dPhiLeadNearSmall[j][i]->Draw("same");
        leadPhiDifFitSmall[j][i]->Draw("same");
      }
    }
    c1.SaveAs( leadPhiDifOut.c_str() );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string leadEtaOut = leadEtaOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      dEtaLeadLarge[j][i]->SetLineColor(j+1);
      dEtaLeadLarge[j][i]->SetMarkerStyle(29);
      dEtaLeadLarge[j][i]->SetMarkerSize(3);
      dEtaLeadLarge[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        std::string outTitle = "Trigger Jet #Delta#eta " + ptBinString[i];
        dEtaLeadLarge[j][i]->SetTitle( outTitle.c_str() );
        dEtaLeadLarge[j][i]->GetXaxis()->SetTitle("#Delta#eta");
        dEtaLeadLarge[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#eta");
        dEtaLeadLarge[j][i]->Draw();
      }
      else {
        dEtaLeadLarge[j][i]->Draw("same");
      }
      if ( ajSplit[j] ) {
        dEtaLeadSmall[j][i]->SetLineColor(j+1+nFiles);
        dEtaLeadSmall[j][i]->SetMarkerStyle(29);
        dEtaLeadSmall[j][i]->SetMarkerSize(3);
        dEtaLeadSmall[j][i]->SetMarkerColor(j+1+nFiles);
        dEtaLeadSmall[j][i]->Draw("same");
        leadEtaFitSmall[j][i]->Draw("same");
      }
    }
    c1.SaveAs( leadEtaOut.c_str() );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string subPhiOut = subPhiOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      dPhiSubLarge[j][i]->SetLineColor(j+1);
      dPhiSubLarge[j][i]->SetMarkerStyle(29);
      dPhiSubLarge[j][i]->SetMarkerSize(3);
      dPhiSubLarge[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        std::string outTitle = "Recoil Jet #Delta#phi " + ptBinString[i];
        dPhiSubLarge[j][i]->SetTitle( outTitle.c_str() );
        dPhiSubLarge[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiSubLarge[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiSubLarge[j][i]->Draw();
      }
      else {
        dPhiSubLarge[j][i]->Draw("same");
      }
      if ( ajSplit[j] ) {
        dPhiSubSmall[j][i]->SetLineColor(j+1+nFiles);
        dPhiSubSmall[j][i]->SetMarkerStyle(29);
        dPhiSubSmall[j][i]->SetMarkerSize(3);
        dPhiSubSmall[j][i]->SetMarkerColor(j+1+nFiles);
        dPhiSubSmall[j][i]->Draw("same");
        subPhiFitSmall[j][i]->Draw("same");
      }
    }
    c1.SaveAs( subPhiOut.c_str() );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string subPhiDifOut = subPhiDifOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      dPhiSubNearLarge[j][i]->SetLineColor(j+1);
      dPhiSubNearLarge[j][i]->SetMarkerStyle(29);
      dPhiSubNearLarge[j][i]->SetMarkerSize(3);
      dPhiSubNearLarge[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        std::string outTitle = "Recoil Jet #Delta#eta subtracted #Delta#phi " + ptBinString[i];
        dPhiSubNearLarge[j][i]->SetTitle( outTitle.c_str() );
        dPhiSubNearLarge[j][i]->GetXaxis()->SetRangeUser(-corrAnalysis::pi/2.0, corrAnalysis::pi/2.0);
        dPhiSubNearLarge[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiSubNearLarge[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiSubNearLarge[j][i]->Draw();
      }
      else {
        dPhiSubNearLarge[j][i]->Draw("same");
      }
      if ( ajSplit[j] ) {
        dPhiSubNearSmall[j][i]->SetLineColor(j+1+nFiles);
        dPhiSubNearSmall[j][i]->SetMarkerStyle(29);
        dPhiSubNearSmall[j][i]->SetMarkerSize(3);
        dPhiSubNearSmall[j][i]->SetMarkerColor(j+1+nFiles);
        dPhiSubNearSmall[j][i]->Draw("same");
        subPhiDifFitSmall[j][i]->Draw("same");
      }
    }
    c1.SaveAs( subPhiDifOut.c_str() );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    TCanvas c1;
    
    std::string subEtaOut = subEtaOutBase + patch::to_string(i) + outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      dEtaSubLarge[j][i]->SetLineColor(j+1);
      dEtaSubLarge[j][i]->SetMarkerStyle(29);
      dEtaSubLarge[j][i]->SetMarkerSize(3);
      dEtaSubLarge[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        std::string outTitle = "Recoil Jet #Delta#eta " + ptBinString[i];
        dEtaSubLarge[j][i]->SetTitle( outTitle.c_str() );
        dEtaSubLarge[j][i]->GetXaxis()->SetTitle("#Delta#eta");
        dEtaSubLarge[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#eta");
        dEtaSubLarge[j][i]->Draw();
      }
      else {
        dEtaSubLarge[j][i]->Draw("same");
      }
      if ( ajSplit[j] ) {
        dEtaSubSmall[j][i]->SetLineColor(j+1+nFiles);
        dEtaSubSmall[j][i]->SetMarkerStyle(29);
        dEtaSubSmall[j][i]->SetMarkerSize(3);
        dEtaSubSmall[j][i]->SetMarkerColor(j+1+nFiles);
        dEtaSubSmall[j][i]->Draw("same");
        subEtaFitSmall[j][i]->Draw("same");
      }
    }
    c1.SaveAs( subEtaOut.c_str() );
  }
  
  // now to get yields
  std::vector<std::vector<double> > leadPhiYieldLarge( nFiles );
  std::vector<std::vector<double> > leadPhiWidthLarge( nFiles );
  std::vector<std::vector<double> > leadPhiErrorLarge( nFiles );
  std::vector<std::vector<double> > leadPhiWidthErrorLarge( nFiles );
  std::vector<std::vector<double> > leadPhiDifYieldLarge( nFiles );
  std::vector<std::vector<double> > leadPhiDifWidthLarge( nFiles );
  std::vector<std::vector<double> > leadPhiDifErrorLarge( nFiles );
  std::vector<std::vector<double> > leadPhiDifWidthErrorLarge( nFiles );
  std::vector<std::vector<double> > leadEtaYieldLarge( nFiles );
  std::vector<std::vector<double> > leadEtaWidthLarge( nFiles );
  std::vector<std::vector<double> > leadEtaErrorLarge( nFiles );
  std::vector<std::vector<double> > leadEtaWidthErrorLarge( nFiles );
  std::vector<std::vector<double> > subPhiYieldLarge( nFiles );
  std::vector<std::vector<double> > subPhiWidthLarge( nFiles );
  std::vector<std::vector<double> > subPhiErrorLarge( nFiles );
  std::vector<std::vector<double> > subPhiWidthErrorLarge( nFiles );
  std::vector<std::vector<double> > subPhiDifYieldLarge( nFiles );
  std::vector<std::vector<double> > subPhiDifWidthLarge( nFiles );
  std::vector<std::vector<double> > subPhiDifErrorLarge( nFiles );
  std::vector<std::vector<double> > subPhiDifWidthErrorLarge( nFiles );
  std::vector<std::vector<double> > subEtaYieldLarge( nFiles );
  std::vector<std::vector<double> > subEtaWidthLarge( nFiles );
  std::vector<std::vector<double> > subEtaErrorLarge( nFiles );
  std::vector<std::vector<double> > subEtaWidthErrorLarge( nFiles );
  
  // and for possible split values...
  std::vector<std::vector<double> > leadPhiYieldSmall( nFiles );
  std::vector<std::vector<double> > leadPhiWidthSmall( nFiles );
  std::vector<std::vector<double> > leadPhiErrorSmall( nFiles );
  std::vector<std::vector<double> > leadPhiWidthErrorSmall( nFiles );
  std::vector<std::vector<double> > leadPhiDifYieldSmall( nFiles );
  std::vector<std::vector<double> > leadPhiDifWidthSmall( nFiles );
  std::vector<std::vector<double> > leadPhiDifErrorSmall( nFiles );
  std::vector<std::vector<double> > leadPhiDifWidthErrorSmall( nFiles );
  std::vector<std::vector<double> > leadEtaYieldSmall( nFiles );
  std::vector<std::vector<double> > leadEtaWidthSmall( nFiles );
  std::vector<std::vector<double> > leadEtaErrorSmall( nFiles );
  std::vector<std::vector<double> > leadEtaWidthErrorSmall( nFiles );
  std::vector<std::vector<double> > subPhiYieldSmall( nFiles );
  std::vector<std::vector<double> > subPhiWidthSmall( nFiles );
  std::vector<std::vector<double> > subPhiErrorSmall( nFiles );
  std::vector<std::vector<double> > subPhiWidthErrorSmall( nFiles );
  std::vector<std::vector<double> > subPhiDifYieldSmall( nFiles );
  std::vector<std::vector<double> > subPhiDifWidthSmall( nFiles );
  std::vector<std::vector<double> > subPhiDifErrorSmall( nFiles );
  std::vector<std::vector<double> > subPhiDifWidthErrorSmall( nFiles );
  std::vector<std::vector<double> > subEtaYieldSmall( nFiles );
  std::vector<std::vector<double> > subEtaWidthSmall( nFiles );
  std::vector<std::vector<double> > subEtaErrorSmall( nFiles );
  std::vector<std::vector<double> > subEtaWidthErrorSmall( nFiles );
  
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiYieldLarge[i].resize( nPtBins );
    leadPhiWidthLarge[i].resize( nPtBins );
    leadPhiErrorLarge[i].resize( nPtBins );
    leadPhiWidthErrorLarge[i].resize( nPtBins );
    leadPhiDifYieldLarge[i].resize( nPtBins );
    leadPhiDifWidthLarge[i].resize( nPtBins );
    leadPhiDifErrorLarge[i].resize( nPtBins );
    leadPhiDifWidthErrorLarge[i].resize( nPtBins );
    leadEtaYieldLarge[i].resize( nPtBins );
    leadEtaWidthLarge[i].resize( nPtBins );
    leadEtaErrorLarge[i].resize( nPtBins );
    leadEtaWidthErrorLarge[i].resize( nPtBins );
    subPhiYieldLarge[i].resize( nPtBins );
    subPhiWidthLarge[i].resize( nPtBins );
    subPhiErrorLarge[i].resize( nPtBins );
    subPhiWidthErrorLarge[i].resize( nPtBins );
    subPhiDifYieldLarge[i].resize( nPtBins );
    subPhiDifWidthLarge[i].resize( nPtBins );
    subPhiDifErrorLarge[i].resize( nPtBins );
    subPhiDifWidthErrorLarge[i].resize( nPtBins );
    subEtaYieldLarge[i].resize( nPtBins );
    subEtaWidthLarge[i].resize( nPtBins );
    subEtaErrorLarge[i].resize( nPtBins );
    subEtaWidthErrorLarge[i].resize( nPtBins );
    
    leadPhiYieldSmall[i].resize( nPtBins );
    leadPhiWidthSmall[i].resize( nPtBins );
    leadPhiErrorSmall[i].resize( nPtBins );
    leadPhiWidthErrorSmall[i].resize( nPtBins );
    leadPhiDifYieldSmall[i].resize( nPtBins );
    leadPhiDifWidthSmall[i].resize( nPtBins );
    leadPhiDifErrorSmall[i].resize( nPtBins );
    leadPhiDifWidthErrorSmall[i].resize( nPtBins );
    leadEtaYieldSmall[i].resize( nPtBins );
    leadEtaWidthSmall[i].resize( nPtBins );
    leadEtaErrorSmall[i].resize( nPtBins );
    leadEtaWidthErrorSmall[i].resize( nPtBins );
    subPhiYieldSmall[i].resize( nPtBins );
    subPhiWidthSmall[i].resize( nPtBins );
    subPhiErrorSmall[i].resize( nPtBins );
    subPhiWidthErrorSmall[i].resize( nPtBins );
    subPhiDifYieldSmall[i].resize( nPtBins );
    subPhiDifWidthSmall[i].resize( nPtBins );
    subPhiDifErrorSmall[i].resize( nPtBins );
    subPhiDifWidthErrorSmall[i].resize( nPtBins );
    subEtaYieldSmall[i].resize( nPtBins );
    subEtaWidthSmall[i].resize( nPtBins );
    subEtaErrorSmall[i].resize( nPtBins );
    subEtaWidthErrorSmall[i].resize( nPtBins );
    for ( int j = 0; j < nPtBins; ++j ) {
      
      // for default
      leadPhiYieldLarge[i][j] = leadPhiFitLarge[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(leadPhiFitLarge[i][j]->GetParameter(3))/ptBinWidth[j];
      leadPhiErrorLarge[i][j] = leadPhiFitLarge[i][j]->GetParError(1);
      leadPhiWidthLarge[i][j] = fabs(leadPhiFitLarge[i][j]->GetParameter(3));
      leadPhiWidthErrorLarge[i][j] = leadPhiFitLarge[i][j]->GetParError(3);
      leadPhiDifYieldLarge[i][j] = leadPhiDifFitLarge[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(leadPhiDifFitLarge[i][j]->GetParameter(3))/ptBinWidth[j];
      leadPhiDifErrorLarge[i][j] = leadPhiDifFitLarge[i][j]->GetParError(1);
      leadPhiDifWidthLarge[i][j] = fabs(leadPhiDifFitLarge[i][j]->GetParameter(3));
      leadPhiDifWidthErrorLarge[i][j] = leadPhiDifFitLarge[i][j]->GetParError(3);
      leadEtaYieldLarge[i][j] = leadEtaFitLarge[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(leadEtaFitLarge[i][j]->GetParameter(3))/ptBinWidth[j];
      leadEtaErrorLarge[i][j] = leadEtaFitLarge[i][j]->GetParError(1);
      leadEtaWidthLarge[i][j] = fabs(leadEtaFitLarge[i][j]->GetParameter(3));
      leadEtaWidthErrorLarge[i][j] = leadEtaFitLarge[i][j]->GetParError(3);
      subPhiYieldLarge[i][j] = subPhiFitLarge[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(subPhiFitLarge[i][j]->GetParameter(3))/ptBinWidth[j];
      subPhiErrorLarge[i][j] = subPhiFitLarge[i][j]->GetParError(1);
      subPhiWidthLarge[i][j] = fabs(subPhiFitLarge[i][j]->GetParameter(3));
      subPhiWidthErrorLarge[i][j] = subPhiFitLarge[i][j]->GetParError(3);
      subPhiDifYieldLarge[i][j] = subPhiDifFitLarge[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(subPhiDifFitLarge[i][j]->GetParameter(3))/ptBinWidth[j];
      subPhiDifErrorLarge[i][j] = subPhiDifFitLarge[i][j]->GetParError(1);
      subPhiDifWidthLarge[i][j] = fabs(subPhiDifFitLarge[i][j]->GetParameter(3));
      subPhiDifWidthErrorLarge[i][j] = subPhiDifFitLarge[i][j]->GetParError(3);
      subEtaYieldLarge[i][j] = subEtaFitLarge[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(subEtaFitLarge[i][j]->GetParameter(3))/ptBinWidth[j];
      subEtaErrorLarge[i][j] = subEtaFitLarge[i][j]->GetParError(1);
      subEtaWidthLarge[i][j] = fabs(subEtaFitLarge[i][j]->GetParameter(3));
      subEtaWidthErrorLarge[i][j] = subEtaFitLarge[i][j]->GetParError(3);
      
      // if we split on aj we have this as well....
      if ( ajSplit[i] ) {
        leadPhiYieldSmall[i][j] = leadPhiFitSmall[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(leadPhiFitSmall[i][j]->GetParameter(3))/ptBinWidth[j];
        leadPhiErrorSmall[i][j] = leadPhiFitSmall[i][j]->GetParError(1);
        leadPhiWidthSmall[i][j] = fabs(leadPhiFitSmall[i][j]->GetParameter(3));
        leadPhiWidthErrorSmall[i][j] = leadPhiFitSmall[i][j]->GetParError(3);
        leadPhiDifYieldSmall[i][j] = leadPhiDifFitSmall[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(leadPhiDifFitSmall[i][j]->GetParameter(3))/ptBinWidth[j];
        leadPhiDifErrorSmall[i][j] = leadPhiDifFitSmall[i][j]->GetParError(1);
        leadPhiDifWidthSmall[i][j] = fabs(leadPhiDifFitSmall[i][j]->GetParameter(3));
        leadPhiDifWidthErrorSmall[i][j] = leadPhiDifFitSmall[i][j]->GetParError(3);
        leadEtaYieldSmall[i][j] = leadEtaFitSmall[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(leadEtaFitSmall[i][j]->GetParameter(3))/ptBinWidth[j];
        leadEtaErrorSmall[i][j] = leadEtaFitSmall[i][j]->GetParError(1);
        leadEtaWidthSmall[i][j] = fabs(leadEtaFitSmall[i][j]->GetParameter(3));
        leadEtaWidthErrorSmall[i][j] = leadEtaFitSmall[i][j]->GetParError(3);
        subPhiYieldSmall[i][j] = subPhiFitSmall[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(subPhiFitSmall[i][j]->GetParameter(3))/ptBinWidth[j];
        subPhiErrorSmall[i][j] = subPhiFitSmall[i][j]->GetParError(1);
        subPhiWidthSmall[i][j] = fabs(subPhiFitSmall[i][j]->GetParameter(3));
        subPhiWidthErrorSmall[i][j] = subPhiFitSmall[i][j]->GetParError(3);
        subPhiDifYieldSmall[i][j] = subPhiDifFitSmall[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(subPhiDifFitSmall[i][j]->GetParameter(3))/ptBinWidth[j];
        subPhiDifErrorSmall[i][j] = subPhiDifFitSmall[i][j]->GetParError(1);
        subPhiDifWidthSmall[i][j] = fabs(subPhiDifFitSmall[i][j]->GetParameter(3));
        subPhiDifWidthErrorSmall[i][j] = subPhiDifFitSmall[i][j]->GetParError(3);
        subEtaYieldSmall[i][j] = subEtaFitSmall[i][j]->GetParameter(1)*sqrt(2*corrAnalysis::pi)*fabs(subEtaFitSmall[i][j]->GetParameter(3))/ptBinWidth[j];
        subEtaErrorSmall[i][j] = subEtaFitSmall[i][j]->GetParError(1);
        subEtaWidthSmall[i][j] = fabs(subEtaFitSmall[i][j]->GetParameter(3));
        subEtaWidthErrorSmall[i][j] = subEtaFitSmall[i][j]->GetParError(3);
      }
    }
  }
  
  // finding the weighted center of the pt bins for future graphs
  std::vector<std::vector<double> > ptBinCenter(nFiles);
  for ( int i = 0; i < nFiles; ++i ) {
    ptBinCenter[i].resize(nPtBins);
    
    for ( int j = 0; j < nPtBins; ++j ) {
      double weightedTotal = 0;
      double entries = 0;
      
      for ( int k = ptBinLo[j]; k <= ptBinHi[j]; ++k ) {
        entries += recombinedPtLead[i]->GetBinContent(k);
        weightedTotal += recombinedPtLead[i]->GetBinContent(k) * recombinedPtLead[i]->GetBinCenter(k);
      }
      ptBinCenter[i][j] = weightedTotal / entries;
    }
  }
  
  
  std::vector<TGraphErrors*> leadPhiGraphLarge( nFiles );
  std::vector<TGraphErrors*> leadPhiWidthGraphLarge( nFiles );
  std::vector<TGraphErrors*> leadPhiDifGraphLarge( nFiles );
  std::vector<TGraphErrors*> leadPhiDifWidthGraphLarge( nFiles );
  std::vector<TGraphErrors*> leadEtaGraphLarge( nFiles );
  std::vector<TGraphErrors*> leadEtaWidthGraphLarge( nFiles );
  std::vector<TGraphErrors*> subPhiGraphLarge( nFiles );
  std::vector<TGraphErrors*> subPhiWidthGraphLarge( nFiles );
  std::vector<TGraphErrors*> subPhiDifGraphLarge( nFiles );
  std::vector<TGraphErrors*> subPhiDifWidthGraphLarge( nFiles );
  std::vector<TGraphErrors*> subEtaGraphLarge( nFiles );
  std::vector<TGraphErrors*> subEtaWidthGraphLarge( nFiles );
  
  std::vector<TGraphErrors*> leadPhiGraphSmall( nFiles );
  std::vector<TGraphErrors*> leadPhiWidthGraphSmall( nFiles );
  std::vector<TGraphErrors*> leadPhiDifGraphSmall( nFiles );
  std::vector<TGraphErrors*> leadPhiDifWidthGraphSmall( nFiles );
  std::vector<TGraphErrors*> leadEtaGraphSmall( nFiles );
  std::vector<TGraphErrors*> leadEtaWidthGraphSmall( nFiles );
  std::vector<TGraphErrors*> subPhiGraphSmall( nFiles );
  std::vector<TGraphErrors*> subPhiWidthGraphSmall( nFiles );
  std::vector<TGraphErrors*> subPhiDifGraphSmall( nFiles );
  std::vector<TGraphErrors*> subPhiDifWidthGraphSmall( nFiles );
  std::vector<TGraphErrors*> subEtaGraphSmall( nFiles );
  std::vector<TGraphErrors*> subEtaWidthGraphSmall( nFiles );
  
  for ( int i = 0; i < nFiles; ++i ) {
    double leadPhiTmpLarge[nPtBins-startPtBin];
    double leadPhiErrLarge[nPtBins-startPtBin];
    double leadPhiWidthTmpLarge[nPtBins-startPtBin];
    double leadPhiWidthErrTmpLarge[nPtBins-startPtBin];
    double leadPhiDifTmpLarge[nPtBins-startPtBin];
    double leadPhiDifErrLarge[nPtBins-startPtBin];
    double leadPhiDifWidthTmpLarge[nPtBins-startPtBin];
    double leadPhiDifWidthErrTmpLarge[nPtBins-startPtBin];
    double leadEtaTmpLarge[nPtBins-startPtBin];
    double leadEtaErrLarge[nPtBins-startPtBin];
    double leadEtaWidthTmpLarge[nPtBins-startPtBin];
    double leadEtaWidthErrTmpLarge[nPtBins-startPtBin];
    double subPhiTmpLarge[nPtBins-startPtBin];
    double subPhiErrLarge[nPtBins-startPtBin];
    double subPhiWidthTmpLarge[nPtBins-startPtBin];
    double subPhiWidthErrTmpLarge[nPtBins-startPtBin];
    double subPhiDifTmpLarge[nPtBins-startPtBin];
    double subPhiDifErrLarge[nPtBins-startPtBin];
    double subPhiDifWidthTmpLarge[nPtBins-startPtBin];
    double subPhiDifWidthErrTmpLarge[nPtBins-startPtBin];
    double subEtaTmpLarge[nPtBins-startPtBin];
    double subEtaErrLarge[nPtBins-startPtBin];
    double subEtaWidthTmpLarge[nPtBins-startPtBin];
    double subEtaWidthErrTmpLarge[nPtBins-startPtBin];
    
    double leadPhiTmpSmall[nPtBins-startPtBin];
    double leadPhiErrSmall[nPtBins-startPtBin];
    double leadPhiWidthTmpSmall[nPtBins-startPtBin];
    double leadPhiWidthErrTmpSmall[nPtBins-startPtBin];
    double leadPhiDifTmpSmall[nPtBins-startPtBin];
    double leadPhiDifErrSmall[nPtBins-startPtBin];
    double leadPhiDifWidthTmpSmall[nPtBins-startPtBin];
    double leadPhiDifWidthErrTmpSmall[nPtBins-startPtBin];
    double leadEtaTmpSmall[nPtBins-startPtBin];
    double leadEtaErrSmall[nPtBins-startPtBin];
    double leadEtaWidthTmpSmall[nPtBins-startPtBin];
    double leadEtaWidthErrTmpSmall[nPtBins-startPtBin];
    double subPhiTmpSmall[nPtBins-startPtBin];
    double subPhiErrSmall[nPtBins-startPtBin];
    double subPhiWidthTmpSmall[nPtBins-startPtBin];
    double subPhiWidthErrTmpSmall[nPtBins-startPtBin];
    double subPhiDifTmpSmall[nPtBins-startPtBin];
    double subPhiDifErrSmall[nPtBins-startPtBin];
    double subPhiDifWidthTmpSmall[nPtBins-startPtBin];
    double subPhiDifWidthErrTmpSmall[nPtBins-startPtBin];
    double subEtaTmpSmall[nPtBins-startPtBin];
    double subEtaErrSmall[nPtBins-startPtBin];
    double subEtaWidthTmpSmall[nPtBins-startPtBin];
    double subEtaWidthErrTmpSmall[nPtBins-startPtBin];
    
    double errX[nPtBins - startPtBin];
    double tmpPtBin[nPtBins - startPtBin];
    for ( int j = 0; j < nPtBins-startPtBin; ++j ) {
      errX[j] = 0;
      tmpPtBin[j] = ptBinCenter[i][j+startPtBin];
      
      leadPhiTmpLarge[j] = leadPhiYieldLarge[i][j+startPtBin];
      leadPhiErrLarge[j] = leadPhiErrorLarge[i][j+startPtBin];
      leadPhiWidthTmpLarge[j] = leadPhiWidthLarge[i][j+startPtBin];
      leadPhiWidthErrTmpLarge[j] = leadPhiWidthErrorLarge[i][j+startPtBin];
      leadPhiDifTmpLarge[j] = leadPhiDifYieldLarge[i][j+startPtBin];
      leadPhiDifErrLarge[j] = leadPhiDifErrorLarge[i][j+startPtBin];
      leadPhiDifWidthTmpLarge[j] = leadPhiDifWidthLarge[i][j+startPtBin];
      leadPhiDifWidthErrTmpLarge[j] = leadPhiDifWidthErrorLarge[i][j+startPtBin];
      leadEtaTmpLarge[j] = leadEtaYieldLarge[i][j+startPtBin];
      leadEtaErrLarge[j] = leadEtaErrorLarge[i][j+startPtBin];
      leadEtaWidthTmpLarge[j] = leadEtaWidthLarge[i][j+startPtBin];
      leadEtaWidthErrTmpLarge[j] = leadEtaWidthErrorLarge[i][j+startPtBin];
      subPhiTmpLarge[j] = subPhiYieldLarge[i][j+startPtBin];
      subPhiErrLarge[j] = subPhiErrorLarge[i][j+startPtBin];
      subPhiWidthTmpLarge[j] = subPhiWidthLarge[i][j+startPtBin];
      subPhiWidthErrTmpLarge[j] = subPhiWidthErrorLarge[i][j+startPtBin];
      subPhiDifTmpLarge[j] = subPhiDifYieldLarge[i][j+startPtBin];
      subPhiDifErrLarge[j] = subPhiDifErrorLarge[i][j+startPtBin];
      subPhiDifWidthTmpLarge[j] = subPhiDifWidthLarge[i][j+startPtBin];
      subPhiDifWidthErrTmpLarge[j] = subPhiDifWidthErrorLarge[i][j+startPtBin];
      subEtaTmpLarge[j] = subEtaYieldLarge[i][j+startPtBin];
      subEtaErrLarge[j] = subEtaErrorLarge[i][j+startPtBin];
      subEtaWidthTmpLarge[j] = subEtaWidthLarge[i][j+startPtBin];
      subEtaWidthErrTmpLarge[j] = subEtaWidthErrorLarge[i][j+startPtBin];
      
      leadPhiTmpSmall[j] = leadPhiYieldSmall[i][j+startPtBin];
      leadPhiErrSmall[j] = leadPhiErrorSmall[i][j+startPtBin];
      leadPhiWidthTmpSmall[j] = leadPhiWidthSmall[i][j+startPtBin];
      leadPhiWidthErrTmpSmall[j] = leadPhiWidthErrorSmall[i][j+startPtBin];
      leadPhiDifTmpSmall[j] = leadPhiDifYieldSmall[i][j+startPtBin];
      leadPhiDifErrSmall[j] = leadPhiDifErrorSmall[i][j+startPtBin];
      leadPhiDifWidthTmpSmall[j] = leadPhiDifWidthSmall[i][j+startPtBin];
      leadPhiDifWidthErrTmpSmall[j] = leadPhiDifWidthErrorSmall[i][j+startPtBin];
      leadEtaTmpSmall[j] = leadEtaYieldSmall[i][j+startPtBin];
      leadEtaErrSmall[j] = leadEtaErrorSmall[i][j+startPtBin];
      leadEtaWidthTmpSmall[j] = leadEtaWidthSmall[i][j+startPtBin];
      leadEtaWidthErrTmpSmall[j] = leadEtaWidthErrorSmall[i][j+startPtBin];
      subPhiTmpSmall[j] = subPhiYieldSmall[i][j+startPtBin];
      subPhiErrSmall[j] = subPhiErrorSmall[i][j+startPtBin];
      subPhiWidthTmpSmall[j] = subPhiWidthSmall[i][j+startPtBin];
      subPhiWidthErrTmpSmall[j] = subPhiWidthErrorSmall[i][j+startPtBin];
      subPhiDifTmpSmall[j] = subPhiDifYieldSmall[i][j+startPtBin];
      subPhiDifErrSmall[j] = subPhiDifErrorSmall[i][j+startPtBin];
      subPhiDifWidthTmpSmall[j] = subPhiDifWidthSmall[i][j+startPtBin];
      subPhiDifWidthErrTmpSmall[j] = subPhiDifWidthErrorSmall[i][j+startPtBin];
      subEtaTmpSmall[j] = subEtaYieldSmall[i][j+startPtBin];
      subEtaErrSmall[j] = subEtaErrorSmall[i][j+startPtBin];
      subEtaWidthTmpSmall[j] = subEtaWidthSmall[i][j+startPtBin];
      subEtaWidthErrTmpSmall[j] = subEtaWidthErrorSmall[i][j+startPtBin];
    }
    
    leadPhiGraphLarge[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadPhiTmpLarge, errX, leadPhiErrLarge );
    leadPhiWidthGraphLarge[i] = new TGraphErrors( nPtBins - startPtBin, tmpPtBin, leadPhiWidthTmpLarge, errX, leadPhiWidthErrTmpLarge );
    leadPhiDifGraphLarge[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadPhiDifTmpLarge, errX, leadPhiDifErrLarge );
    leadPhiDifWidthGraphLarge[i] = new TGraphErrors( nPtBins - startPtBin, tmpPtBin, leadPhiDifWidthTmpLarge, errX, leadPhiDifWidthErrTmpLarge );
    leadEtaGraphLarge[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadEtaTmpLarge, errX, leadEtaErrLarge );
    leadEtaWidthGraphLarge[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadEtaWidthTmpLarge, errX, leadEtaWidthErrTmpLarge );
    subPhiGraphLarge[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subPhiTmpLarge, errX, subPhiErrLarge );
    subPhiWidthGraphLarge[i] = new TGraphErrors(nPtBins-startPtBin, tmpPtBin, subPhiWidthTmpLarge, errX, subPhiWidthErrTmpLarge );
    subPhiDifGraphLarge[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subPhiDifTmpLarge, errX, subPhiDifErrLarge );
    subPhiDifWidthGraphLarge[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subPhiDifWidthTmpLarge, errX, subPhiDifWidthErrTmpLarge );
    subEtaGraphLarge[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subEtaTmpLarge, errX, subEtaErrLarge );
    subEtaWidthGraphLarge[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subEtaWidthTmpLarge, errX, subEtaWidthErrTmpLarge );

    
    leadPhiGraphSmall[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadPhiTmpSmall, errX, leadPhiErrSmall );
    leadPhiWidthGraphSmall[i] = new TGraphErrors( nPtBins - startPtBin, tmpPtBin, leadPhiWidthTmpSmall, errX, leadPhiWidthErrTmpSmall );
    leadPhiDifGraphSmall[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadPhiDifTmpSmall, errX, leadPhiDifErrSmall );
    leadPhiDifWidthGraphSmall[i] = new TGraphErrors( nPtBins - startPtBin, tmpPtBin, leadPhiDifWidthTmpSmall, errX, leadPhiDifWidthErrTmpSmall );
    leadEtaGraphSmall[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadEtaTmpSmall, errX, leadEtaErrSmall );
    leadEtaWidthGraphSmall[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, leadEtaWidthTmpSmall, errX, leadEtaWidthErrTmpSmall );
    subPhiGraphSmall[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subPhiTmpSmall, errX, subPhiErrSmall );
    subPhiWidthGraphSmall[i] = new TGraphErrors(nPtBins-startPtBin, tmpPtBin, subPhiWidthTmpSmall, errX, subPhiWidthErrTmpSmall );
    subPhiDifGraphSmall[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subPhiDifTmpSmall, errX, subPhiDifErrSmall );
    subPhiDifWidthGraphSmall[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subPhiDifWidthTmpSmall, errX, subPhiDifWidthErrTmpSmall );
    subEtaGraphSmall[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subEtaTmpSmall, errX, subEtaErrSmall );
    subEtaWidthGraphSmall[i] = new TGraphErrors(nPtBins - startPtBin, tmpPtBin, subEtaWidthTmpSmall, errX, subEtaWidthErrTmpSmall );

  }
  
  // print out the yields
  
  TCanvas* c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiGraphLarge[i]->SetLineColor(i+1);
    leadPhiGraphLarge[i]->SetMarkerStyle(29);
    leadPhiGraphLarge[i]->SetMarkerSize(3);
    leadPhiGraphLarge[i]->SetMarkerColor(i+1);
    leadPhiGraphLarge[i]->SetTitle("Trigger Jet - #Delta#phi Fit Yield");
    leadPhiGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    leadPhiGraphLarge[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    leadPhiGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0)
      leadPhiGraphLarge[i]->Draw();
    else
      leadPhiGraphLarge[i]->Draw("P");
    
    if ( ajSplit[i] ) {
      leadPhiGraphSmall[i]->SetLineColor(i+1 + nFiles);
      leadPhiGraphSmall[i]->SetMarkerStyle(29);
      leadPhiGraphSmall[i]->SetMarkerSize(3);
      leadPhiGraphSmall[i]->SetMarkerColor(i+1 + nFiles);
      leadPhiGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/leadphiyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiDifGraphLarge[i]->SetLineColor(i+1);
    leadPhiDifGraphLarge[i]->SetMarkerStyle(29);
    leadPhiDifGraphLarge[i]->SetMarkerSize(3);
    leadPhiDifGraphLarge[i]->SetMarkerColor(i+1);
    leadPhiDifGraphLarge[i]->SetTitle("Trigger Jet - #eta Subtracted #Delta#phi Fit Yield");
    leadPhiDifGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    leadPhiDifGraphLarge[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    leadPhiDifGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0)
      leadPhiDifGraphLarge[i]->Draw();
    else
      leadPhiDifGraphLarge[i]->Draw("P");
    if ( ajSplit[i] ) {
      leadPhiDifGraphSmall[i]->SetLineColor(i+1+nFiles);
      leadPhiDifGraphSmall[i]->SetMarkerStyle(29);
      leadPhiDifGraphSmall[i]->SetMarkerSize(3);
      leadPhiDifGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      leadPhiDifGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/leadphidifyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadEtaGraphLarge[i]->SetLineColor(i+1);
    leadEtaGraphLarge[i]->SetMarkerStyle(29);
    leadEtaGraphLarge[i]->SetMarkerSize(3);
    leadEtaGraphLarge[i]->SetMarkerColor(i+1);
    leadEtaGraphLarge[i]->SetTitle("Trigger Jet - #Delta#eta Fit Yield");
    leadEtaGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    leadEtaGraphLarge[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    leadEtaGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0)
      leadEtaGraphLarge[i]->Draw();
    else
      leadEtaGraphLarge[i]->Draw("P");
    if ( ajSplit[i] ) {
      leadEtaGraphSmall[i]->SetLineColor(i+1+nFiles);
      leadEtaGraphSmall[i]->SetMarkerStyle(29);
      leadEtaGraphSmall[i]->SetMarkerSize(3);
      leadEtaGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      leadEtaGraphSmall[i]->Draw("P");

    }
  }
  c1->SaveAs("tmp/leadetayield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiGraphLarge[i]->SetLineColor(i+1);
    subPhiGraphLarge[i]->SetMarkerStyle(29);
    subPhiGraphLarge[i]->SetMarkerSize(3);
    subPhiGraphLarge[i]->SetMarkerColor(i+1);
    subPhiGraphLarge[i]->SetTitle("Recoil Jet - #Delta#phi Fit Yield");
    subPhiGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    subPhiGraphLarge[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    subPhiGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0)
      subPhiGraphLarge[i]->Draw();
    else
      subPhiGraphLarge[i]->Draw("P");
    if ( ajSplit[i] ) {
      subPhiGraphSmall[i]->SetLineColor(i+1+nFiles);
      subPhiGraphSmall[i]->SetMarkerStyle(29);
      subPhiGraphSmall[i]->SetMarkerSize(3);
      subPhiGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      subPhiGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/subphiyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiDifGraphLarge[i]->SetLineColor(i+1);
    subPhiDifGraphLarge[i]->SetMarkerStyle(29);
    subPhiDifGraphLarge[i]->SetMarkerSize(3);
    subPhiDifGraphLarge[i]->SetMarkerColor(i+1);
    subPhiDifGraphLarge[i]->SetTitle("Recoil Jet - #eta Subtracted #Delta#phi Fit Yield");
    subPhiDifGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    subPhiDifGraphLarge[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    subPhiDifGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0)
      subPhiDifGraphLarge[i]->Draw();
    else
      subPhiDifGraphLarge[i]->Draw("P");
    if ( ajSplit[i] ) {
      subPhiDifGraphSmall[i]->SetLineColor(i+1+nFiles);
      subPhiDifGraphSmall[i]->SetMarkerStyle(29);
      subPhiDifGraphSmall[i]->SetMarkerSize(3);
      subPhiDifGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      subPhiDifGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/subphidifyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subEtaGraphLarge[i]->SetLineColor(i+1);
    subEtaGraphLarge[i]->SetMarkerStyle(29);
    subEtaGraphLarge[i]->SetMarkerSize(3);
    subEtaGraphLarge[i]->SetMarkerColor(i+1);
    subEtaGraphLarge[i]->SetTitle("Recoil Jet - #Delta#eta Fit Yield");
    subEtaGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    subEtaGraphLarge[i]->GetYaxis()->SetTitle("1/N_{dijet}dN/dp_{T}");
    subEtaGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 5 );
    if ( i == 0) {
      subEtaGraphLarge[i]->Draw();
    }
    else {
      subEtaGraphLarge[i]->Draw("P");
    }
    if ( ajSplit[i] ) {
      subEtaGraphSmall[i]->SetLineColor(i+1+nFiles);
      subEtaGraphSmall[i]->SetMarkerStyle(29);
      subEtaGraphSmall[i]->SetMarkerSize(3);
      subEtaGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      subEtaGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/subetayield.pdf");
  
  // now print out widths
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiWidthGraphLarge[i]->SetLineColor(i+1);
    leadPhiWidthGraphLarge[i]->SetMarkerStyle(29);
    leadPhiWidthGraphLarge[i]->SetMarkerSize(3);
    leadPhiWidthGraphLarge[i]->SetMarkerColor(i+1);
    leadPhiWidthGraphLarge[i]->SetTitle("Trigger Jet #Delta#phi - Gaussian Widths");
    leadPhiWidthGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    leadPhiWidthGraphLarge[i]->GetYaxis()->SetTitle("Width");
    leadPhiWidthGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      leadPhiWidthGraphLarge[i]->Draw();
    else
      leadPhiWidthGraphLarge[i]->Draw("P");
    if ( ajSplit[i] ) {
      leadPhiWidthGraphSmall[i]->SetLineColor(i+1+nFiles);
      leadPhiWidthGraphSmall[i]->SetMarkerStyle(29);
      leadPhiWidthGraphSmall[i]->SetMarkerSize(3);
      leadPhiWidthGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      leadPhiWidthGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/leadphiwidth.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiDifWidthGraphLarge[i]->SetLineColor(i+1);
    leadPhiDifWidthGraphLarge[i]->SetMarkerStyle(29);
    leadPhiDifWidthGraphLarge[i]->SetMarkerSize(3);
    leadPhiDifWidthGraphLarge[i]->SetMarkerColor(i+1);
    leadPhiDifWidthGraphLarge[i]->SetTitle("Trigger Jet #Delta#eta Subtracted #Delta#phi - Gaussian Widths");
    leadPhiDifWidthGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    leadPhiDifWidthGraphLarge[i]->GetYaxis()->SetTitle("Width");
    leadPhiDifWidthGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      leadPhiDifWidthGraphLarge[i]->Draw();
    else
      leadPhiDifWidthGraphLarge[i]->Draw("P");
    if ( ajSplit[i] ) {
      leadPhiDifWidthGraphSmall[i]->SetLineColor(i+1+nFiles);
      leadPhiDifWidthGraphSmall[i]->SetMarkerStyle(29);
      leadPhiDifWidthGraphSmall[i]->SetMarkerSize(3);
      leadPhiDifWidthGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      leadPhiDifWidthGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/leadphidifwidth.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadEtaWidthGraphLarge[i]->SetLineColor(i+1);
    leadEtaWidthGraphLarge[i]->SetMarkerStyle(29);
    leadEtaWidthGraphLarge[i]->SetMarkerSize(3);
    leadEtaWidthGraphLarge[i]->SetMarkerColor(i+1);
    leadEtaWidthGraphLarge[i]->SetTitle("Trigger Jet #Delta#eta - Gaussian Widths");
    leadEtaWidthGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    leadEtaWidthGraphLarge[i]->GetYaxis()->SetTitle("Width");
    leadEtaWidthGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      leadEtaWidthGraphLarge[i]->Draw();
    else
      leadEtaWidthGraphLarge[i]->Draw("P");
    if ( ajSplit[i] ) {
      leadEtaWidthGraphSmall[i]->SetLineColor(i+1+nFiles);
      leadEtaWidthGraphSmall[i]->SetMarkerStyle(29);
      leadEtaWidthGraphSmall[i]->SetMarkerSize(3);
      leadEtaWidthGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      leadEtaWidthGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/leadetawidth.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiWidthGraphLarge[i]->SetLineColor(i+1);
    subPhiWidthGraphLarge[i]->SetMarkerStyle(29);
    subPhiWidthGraphLarge[i]->SetMarkerSize(3);
    subPhiWidthGraphLarge[i]->SetMarkerColor(i+1);
    subPhiWidthGraphLarge[i]->SetTitle("Recoil Jet #Delta#phi - Gaussian Widths");
    subPhiWidthGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    subPhiWidthGraphLarge[i]->GetYaxis()->SetTitle("Width");
    subPhiWidthGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      subPhiWidthGraphLarge[i]->Draw();
    else
      subPhiWidthGraphLarge[i]->Draw("P");
    if ( ajSplit[i] ) {
      subPhiWidthGraphSmall[i]->SetLineColor(i+1+nFiles);
      subPhiWidthGraphSmall[i]->SetMarkerStyle(29);
      subPhiWidthGraphSmall[i]->SetMarkerSize(3);
      subPhiWidthGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      subPhiWidthGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/subphiwidth.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiDifWidthGraphLarge[i]->SetLineColor(i+1);
    subPhiDifWidthGraphLarge[i]->SetMarkerStyle(29);
    subPhiDifWidthGraphLarge[i]->SetMarkerSize(3);
    subPhiDifWidthGraphLarge[i]->SetMarkerColor(i+1);
    subPhiDifWidthGraphLarge[i]->SetTitle("Recoil Jet #Delta#eta Subtracted #Delta#phi - Gaussian Widths");
    subPhiDifWidthGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    subPhiDifWidthGraphLarge[i]->GetYaxis()->SetTitle("Width");
    subPhiDifWidthGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      subPhiDifWidthGraphLarge[i]->Draw();
    else
      subPhiDifWidthGraphLarge[i]->Draw("P");
    if ( ajSplit[i] ) {
      subPhiDifWidthGraphSmall[i]->SetLineColor(i+1+nFiles);
      subPhiDifWidthGraphSmall[i]->SetMarkerStyle(29);
      subPhiDifWidthGraphSmall[i]->SetMarkerSize(3);
      subPhiDifWidthGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      subPhiDifWidthGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/subphidifwidth.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subEtaWidthGraphLarge[i]->SetLineColor(i+1);
    subEtaWidthGraphLarge[i]->SetMarkerStyle(29);
    subEtaWidthGraphLarge[i]->SetMarkerSize(3);
    subEtaWidthGraphLarge[i]->SetMarkerColor(i+1);
    subEtaWidthGraphLarge[i]->SetTitle("Recoil Jet #Delta#eta - Gaussian Widths");
    subEtaWidthGraphLarge[i]->GetXaxis()->SetTitle("p_{T}");
    subEtaWidthGraphLarge[i]->GetYaxis()->SetTitle("Width");
    subEtaWidthGraphLarge[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      subEtaWidthGraphLarge[i]->Draw();
    else
      subEtaWidthGraphLarge[i]->Draw("P");
    if ( ajSplit[i] ) {
      subEtaWidthGraphSmall[i]->SetLineColor(i+1+nFiles);
      subEtaWidthGraphSmall[i]->SetMarkerStyle(29);
      subEtaWidthGraphSmall[i]->SetMarkerSize(3);
      subEtaWidthGraphSmall[i]->SetMarkerColor(i+1+nFiles);
      subEtaWidthGraphSmall[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/subetawidth.pdf");
  
  
  // now put leading and subleading yields on their own plots
  for ( int i = 0; i < nFiles; ++i ) {
    c1 = new TCanvas();
    leadPhiGraphLarge[i]->SetLineColor(1);
    leadPhiGraphLarge[i]->SetMarkerStyle(29);
    leadPhiGraphLarge[i]->SetMarkerSize(2);
    leadPhiGraphLarge[i]->SetMarkerColor(1);
    leadPhiDifGraphLarge[i]->SetLineColor(2);
    leadPhiDifGraphLarge[i]->SetMarkerStyle(20);
    leadPhiDifGraphLarge[i]->SetMarkerSize(2);
    leadPhiDifGraphLarge[i]->SetMarkerColor(2);
    leadEtaGraphLarge[i]->SetLineColor(3);
    leadEtaGraphLarge[i]->SetMarkerStyle(21);
    leadEtaGraphLarge[i]->SetMarkerSize(2);
    leadEtaGraphLarge[i]->SetMarkerColor(3);
    
    leadPhiGraphSmall[i]->SetLineColor(4);
    leadPhiGraphSmall[i]->SetMarkerStyle(29);
    leadPhiGraphSmall[i]->SetMarkerSize(2);
    leadPhiGraphSmall[i]->SetMarkerColor(4);
    leadPhiDifGraphSmall[i]->SetLineColor(5);
    leadPhiDifGraphSmall[i]->SetMarkerStyle(20);
    leadPhiDifGraphSmall[i]->SetMarkerSize(2);
    leadPhiDifGraphSmall[i]->SetMarkerColor(5);
    leadEtaGraphSmall[i]->SetLineColor(6);
    leadEtaGraphSmall[i]->SetMarkerStyle(21);
    leadEtaGraphSmall[i]->SetMarkerSize(2);
    leadEtaGraphSmall[i]->SetMarkerColor(7);
    
    leadPhiGraphLarge[i]->Draw();
    leadPhiDifGraphLarge[i]->Draw("P");
    leadEtaGraphLarge[i]->Draw("P");
    leadPhiGraphSmall[i]->Draw("P");
    leadPhiDifGraphSmall[i]->Draw("P");
    leadEtaGraphSmall[i]->Draw("P");
    
    TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->SetHeader("Trigger Jet Yields"); // option "C" allows to center the header
    leg->AddEntry(leadPhiGraphLarge[i],"#Delta#phi |Aj| > 0.2","lep");
    leg->AddEntry(leadPhiDifGraphLarge[i],"#Delta#phi #Delta#eta subtracted |Aj| > 0.2","lep");
    leg->AddEntry(leadEtaGraphLarge[i],"#Delta#eta |Aj| > 0.2","lep");
    leg->AddEntry(leadPhiGraphSmall[i],"#Delta#phi |Aj| < 0.2","lep");
    leg->AddEntry(leadPhiDifGraphSmall[i],"#Delta#phi #Delta#eta subtracted |Aj| < 0.2","lep");
    leg->AddEntry(leadEtaGraphSmall[i],"#Delta#eta |Aj| < 0.2","lep");
    leg->Draw();
    
    std::string graphOutName = "tmp/graph_out_"+analysisNames[i]+"_lead.pdf";
    
    c1->SaveAs( graphOutName.c_str() );
    
    c1 = new TCanvas();
    subPhiGraphLarge[i]->SetLineColor(1);
    subPhiGraphLarge[i]->SetMarkerStyle(29);
    subPhiGraphLarge[i]->SetMarkerSize(2);
    subPhiGraphLarge[i]->SetMarkerColor(1);
    subPhiDifGraphLarge[i]->SetLineColor(2);
    subPhiDifGraphLarge[i]->SetMarkerStyle(20);
    subPhiDifGraphLarge[i]->SetMarkerSize(2);
    subPhiDifGraphLarge[i]->SetMarkerColor(2);
    subEtaGraphLarge[i]->SetLineColor(3);
    subEtaGraphLarge[i]->SetMarkerStyle(21);
    subEtaGraphLarge[i]->SetMarkerSize(2);
    subEtaGraphLarge[i]->SetMarkerColor(3);

    subPhiGraphSmall[i]->SetLineColor(4);
    subPhiGraphSmall[i]->SetMarkerStyle(29);
    subPhiGraphSmall[i]->SetMarkerSize(2);
    subPhiGraphSmall[i]->SetMarkerColor(4);
    subPhiDifGraphSmall[i]->SetLineColor(5);
    subPhiDifGraphSmall[i]->SetMarkerStyle(20);
    subPhiDifGraphSmall[i]->SetMarkerSize(2);
    subPhiDifGraphSmall[i]->SetMarkerColor(5);
    subEtaGraphSmall[i]->SetLineColor(6);
    subEtaGraphSmall[i]->SetMarkerStyle(21);
    subEtaGraphSmall[i]->SetMarkerSize(2);
    subEtaGraphSmall[i]->SetMarkerColor(6);

    subPhiGraphLarge[i]->Draw();
    subPhiDifGraphLarge[i]->Draw("P");
    subEtaGraphLarge[i]->Draw("P");
    subPhiGraphSmall[i]->Draw("P");
    subPhiDifGraphSmall[i]->Draw("P");
    subEtaGraphSmall[i]->Draw("P");
    
    leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->SetHeader("Recoil Jet Yields"); // option "C" allows to center the header
    leg->AddEntry(subPhiGraphLarge[i],"#Delta#phi |Aj| > 0.2","lep");
    leg->AddEntry(subPhiDifGraphLarge[i],"#Delta#phi #Delta#eta subtracted |Aj| > 0.2","lep");
    leg->AddEntry(subEtaGraphLarge[i],"#Delta#eta |Aj| > 0.2","lep");
    leg->AddEntry(subPhiGraphSmall[i],"#Delta#phi |Aj| < 0.2","lep");
    leg->AddEntry(subPhiDifGraphSmall[i],"#Delta#phi #Delta#eta subtracted |Aj| < 0.2","lep");
    leg->AddEntry(subEtaGraphSmall[i],"#Delta#eta |Aj| < 0.2","lep");
    leg->Draw();
    
    graphOutName = "tmp/graph_out_"+analysisNames[i]+"_sub.pdf";
    
    c1->SaveAs( graphOutName.c_str() );
    
    // joern asked for dphi without subtraction overlayed
    graphOutName = "tmp/dphi_trig_recoil_yield_"+analysisNames[i]+".pdf";
    c1 = new TCanvas();
    leadPhiGraphLarge[i]->SetTitle("");
    leadPhiGraphLarge[i]->SetLineColor(1);
    leadPhiGraphLarge[i]->SetMarkerColor(1);
    subPhiGraphLarge[i]->SetLineColor(7);
    subPhiGraphLarge[i]->SetMarkerColor(7);
    subPhiGraphLarge[i]->SetMarkerStyle(20);
    
    leadPhiGraphSmall[i]->SetTitle("");
    leadPhiGraphSmall[i]->SetLineColor(2);
    leadPhiGraphSmall[i]->SetMarkerColor(2);
    subPhiGraphSmall[i]->SetLineColor(8);
    subPhiGraphSmall[i]->SetMarkerColor(8);
    subPhiGraphSmall[i]->SetMarkerStyle(21);
    
    leadPhiGraphLarge[i]->Draw();
    subPhiGraphLarge[i]->Draw("P");
    leadPhiGraphSmall[i]->Draw("P");
    subPhiGraphSmall[i]->Draw("P");

    
    leg = new TLegend(0.65,0.7,0.9,0.9);
    leg->SetHeader("#Delta#phi Jet Yields");
    leg->AddEntry(leadPhiGraphLarge[i],"Trigger |Aj| > 0.2","lep");
    leg->AddEntry(subPhiGraphLarge[i],"Recoil |Aj| > 0.2","lep");
    leg->AddEntry(leadPhiGraphSmall[i],"Trigger |Aj| < 0.2","lep");
    leg->AddEntry(subPhiGraphSmall[i],"Recoil |Aj| < 0.2","lep");
    leg->Draw();
    
    c1->SaveAs( graphOutName.c_str() );

  }

  return 0;
}
