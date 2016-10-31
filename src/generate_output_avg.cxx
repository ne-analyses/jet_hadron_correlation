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
  std::vector<std::vector<std::vector<TH3D*> > > mixCentVz;
  std::vector<std::vector<std::vector<TH3D*> > > subCentVz;
  std::vector<std::vector<std::vector<TH3D*> > > mixSubCentVz;
  corrCentVz.resize( nFiles );
  mixCentVz.resize( nFiles );
  subCentVz.resize( nFiles );
  mixSubCentVz.resize( nFiles );
  for ( int i = 0; i < nFiles; ++i ) {
    corrCentVz[i].resize( corrAnalysis::binsCentrality );
    mixCentVz[i].resize( corrAnalysis::binsCentrality );
    subCentVz[i].resize( corrAnalysis::binsCentrality );
    mixSubCentVz[i].resize( corrAnalysis::binsCentrality );
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      corrCentVz[i][j].resize( corrAnalysis::binsVz );
      mixCentVz[i][j].resize( corrAnalysis::binsVz );
      subCentVz[i][j].resize( corrAnalysis::binsVz );
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
        corrCentVz[i][j][k] = (TH3D*) corrFiles[i]->Get( corrDifInitName.c_str() );
        corrCentVz[i][j][k]->SetName( corrDifBaseName.c_str() );
        mixCentVz[i][j][k] = (TH3D*) mixFiles[i]->Get( mixDifInitName.c_str() );
        mixCentVz[i][j][k]->SetName( mixDifBaseName.c_str() );
        subCentVz[i][j][k] = (TH3D*) corrFiles[i]->Get( subDifInitName.c_str() );
        subCentVz[i][j][k]->SetName( subDifBaseName.c_str() );
        mixSubCentVz[i][j][k] = (TH3D*) mixFiles[i]->Get( mixSubDifInitName.c_str() );
        mixSubCentVz[i][j][k]->SetName( mixSubDifBaseName.c_str() );
      }
  }
  
  // setup for 2d projections along pt axis
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > corrCentVzPt;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > mixCentVzPt;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > subCentVzPt;
  std::vector<std::vector<std::vector<std::vector<TH2D*> > > > mixSubCentVzPt;
  corrCentVzPt.resize( nFiles );
  mixCentVzPt.resize( nFiles );
  subCentVzPt.resize( nFiles );
  mixSubCentVzPt.resize( nFiles );
  for ( int i = 0; i < nFiles; ++i ) {
    corrCentVzPt[i].resize( corrAnalysis::binsCentrality );
    mixCentVzPt[i].resize( corrAnalysis::binsCentrality );
    subCentVzPt[i].resize( corrAnalysis::binsCentrality );
    mixSubCentVzPt[i].resize( corrAnalysis::binsCentrality );
    
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      corrCentVzPt[i][j].resize( corrAnalysis::binsVz );
      mixCentVzPt[i][j].resize( corrAnalysis::binsVz );
      subCentVzPt[i][j].resize( corrAnalysis::binsVz );
      mixSubCentVzPt[i][j].resize( corrAnalysis::binsVz );
      for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
        corrCentVzPt[i][j][k].resize( nPtBins );
        mixCentVzPt[i][j][k].resize( nPtBins );
        subCentVzPt[i][j][k].resize( nPtBins );
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
          mixCentVz[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          subCentVz[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          mixSubCentVz[i][j][k]->GetZaxis()->SetRange( ptBinLo[l], ptBinHi[l] );
          
          corrCentVzPt[i][j][k][l] = (TH2D*) ((TH2D*) corrCentVz[i][j][k]->Project3D( "YX" ))->Clone();
          mixCentVzPt[i][j][k][l] = (TH2D*) ((TH2D*) mixCentVz[i][j][k]->Project3D( "YX" ))->Clone();
          subCentVzPt[i][j][k][l] = (TH2D*) ((TH2D*) subCentVz[i][j][k]->Project3D("YX"))->Clone();
          mixSubCentVzPt[i][j][k][l] = (TH2D*) ((TH2D*) mixSubCentVz[i][j][k]->Project3D("YX"))->Clone();
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
    TH1D* vzBins = nEvents[i]->ProjectionY();
    
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
  std::vector<std::vector<TH2D*> > recombinedCorr;
  std::vector<std::vector<TH2D*> > recombinedPre;
  std::vector<std::vector<TH2D*> > recombinedSub;
  std::vector<std::vector<TH2D*> > recombinedSubPre;
  recombinedCorr.resize( nFiles );
  recombinedPre.resize( nFiles );
  recombinedSub.resize( nFiles );
  recombinedSubPre.resize( nFiles );
  
  for (int i = 0; i < nFiles; ++ i ) {
    
    recombinedCorr[i].resize( nPtBins );
    recombinedPre[i].resize( nPtBins );
    recombinedSub[i].resize( nPtBins );
    recombinedSubPre[i].resize( nPtBins );
    
    for ( int l = 0; l < nPtBins; ++l ) {
      
      std::string corrName = analysisNames[i] + " " + ptBinString[l];
      std::string preName = "pre_" + analysisNames[i] + " " + ptBinString[l];
      std::string subName = analysisNames[i] + "_sub " + ptBinString[l];
      std::string subPreName = "pre_" + analysisNames[i] + "_sub " + ptBinString[l];
      
      recombinedCorr[i][l] = new TH2D( corrName.c_str(), corrName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      recombinedPre[i][l] = new TH2D( preName.c_str(), preName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      recombinedSub[i][l] = new TH2D( subName.c_str(), subName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      recombinedSubPre[i][l] = new TH2D( subPreName.c_str(), subPreName.c_str(), corrAnalysis::binsEta, corrAnalysis::dEtaLowEdge, corrAnalysis::dEtaHighEdge, corrAnalysis::binsPhi, corrAnalysis::phiLowEdge, corrAnalysis::phiHighEdge );
      
      if ( l <= 2 ) {
      
        for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
          for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
            if ( weightedMix[i][l]->GetEntries() != 0 && corrCentVzPt[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedPre[i][l]->Add( corrCentVzPt[i][j][k][l] );
              
              corrCentVzPt[i][j][k][l]->Divide( weightedMix[i][l] );
            
              recombinedCorr[i][l]->Add( corrCentVzPt[i][j][k][l] );
            }
            if ( subCentVzPt[i][j][k][l]->GetEntries() != 0 && weightedSub[i][l]->GetEntries() != 0 ) {
              
              recombinedSubPre[i][l]->Add( subCentVzPt[i][j][k][l] );
              
              subCentVzPt[i][j][k][l]->Divide( weightedSub[i][l] );
              
              recombinedSub[i][l]->Add( subCentVzPt[i][j][k][l] );
              
            }
          }
        }
      }
      
      else {
        for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
          for ( int k = 0; k < corrAnalysis::binsVz; ++k ) {
            if ( weightedMix[i][2]->GetEntries() != 0 && corrCentVzPt[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedPre[i][l]->Add( corrCentVzPt[i][j][k][l] );
              
              corrCentVzPt[i][j][k][l]->Divide( weightedMix[i][2] );
              
              recombinedCorr[i][l]->Add( corrCentVzPt[i][j][k][l] );
            }
            if ( weightedSub[i][2]->GetEntries() != 0 && subCentVzPt[i][j][k][l]->GetEntries() != 0 ) {
              
              recombinedSubPre[i][l]->Add( subCentVzPt[i][j][k][l] );
              
              subCentVzPt[i][j][k][l]->Divide( weightedSub[i][2] );
              
              recombinedSub[i][l]->Add( subCentVzPt[i][j][k][l] );
            }
          }
        }
      }
    }
  }
  
  // get the reduced eta and phi ranges for projections
  double etaMax = 1.3;
  double etaMin = -1.3;
  double etaNearMin = etaMin/2.0;
  double etaNearMax = etaMax/2.0;
  double phiMin = -corrAnalysis::pi/2.0;
  double phiMinClose = -0.8;
  double phiMaxClose = 0.8;
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
      
      // save the 2D histograms
      std::string leadOutName = "tmp/lead2d_" + analysisNames[i] +"_pt_"+ patch::to_string(j) + ".pdf";
      std::string subOutName = "tmp/sub2d_" + analysisNames[i] +"_pt_" + patch::to_string(j) + ".pdf";
      
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
      
    }
  }
  
  // now  first overlay and output,
  // then subtract near from far eta regions
  for ( int i = 0; i < nFiles; ++i )
    for ( int j = 0; j < nPtBins; ++j ) {
      std::string leadPhiName = "tmp/lead_phi_near_far_"+analysisNames[i]+"_pt_"+patch::to_string(j)+".pdf";
      std::string subPhiName = "tmp/sub_phi_near_far_"+analysisNames[i]+"_pt_"+patch::to_string(j)+".pdf";
      
      TCanvas c1;
      dPhiLeadNear[i][j]->Draw();
      dPhiLeadFar[i][j]->Draw("SAME");
      c1.SaveAs( leadPhiName.c_str() );
      dPhiSubNear[i][j]->Draw();
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
  
  // do a first, temporary fit to remove background
  for ( int i = 0; i < nFiles; ++i ) {
    for ( int j = 0; j < nPtBins; ++j ) {
      std::string dPhiLeadName = "tmp_fit_lead_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubName = "tmp_fit_sub_phi_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiLeadNameDif = "tmp_fit_lead_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dPhiSubNameDif = "tmp_fit_sub_phi_dif_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaLeadName = "tmp_fit_lead_eta_" + patch::to_string(i) + patch::to_string(j);
      std::string dEtaSubName = "tmp_fit_sub_eta_" + patch::to_string(i) + patch::to_string(j);
      
      TF1* leadPhiInitFit = new TF1( dPhiLeadName.c_str(), phiForm.c_str(), phiMin, phiMax );
      leadPhiInitFit->FixParameter( 2, 0 );
      leadPhiInitFit->SetParameter( 5, corrAnalysis::pi );
      leadPhiInitFit->SetParameter( 3, 0.2 );
      leadPhiInitFit->SetParameter( 6, 0.2 );
      
      TF1* leadPhiDifInitFit = new TF1( dPhiLeadNameDif.c_str(), phiForm.c_str(), phiMin, phiMax );
      leadPhiDifInitFit->FixParameter( 2, 0 );
      leadPhiDifInitFit->SetParameter( 5, corrAnalysis::pi );
      leadPhiDifInitFit->SetParameter( 3, 0.2 );
      leadPhiDifInitFit->SetParameter( 6, 0.2 );
      
      TF1* subPhiInitFit = new TF1( dPhiSubName.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiInitFit->FixParameter( 2, 0 );
      subPhiInitFit->SetParameter( 5, corrAnalysis::pi );
      subPhiInitFit->SetParameter( 3, 0.2 );
      subPhiInitFit->SetParameter( 6, 0.2 );
      
      TF1* subPhiDifInitFit = new TF1( dPhiSubNameDif.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiDifInitFit->FixParameter( 2, 0 );
      subPhiDifInitFit->SetParameter( 5, corrAnalysis::pi );
      subPhiDifInitFit->SetParameter( 3, 0.2 );
      subPhiDifInitFit->SetParameter( 6, 0.2 );
      
      TF1* leadEtaInitFit = new TF1( dEtaLeadName.c_str(), etaForm.c_str(), etaMin, etaMax );
      leadEtaInitFit->FixParameter( 2, 0 );
      leadEtaInitFit->SetParameter( 3, 0.2 );
      
      TF1* subEtaInitFit = new TF1( dEtaSubName.c_str(), etaForm.c_str(), etaMin, etaMax );
      subEtaInitFit->FixParameter( 2, 0 );
      subEtaInitFit->SetParameter( 3, 0.2 );
      
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
      leadPhiFit[i][j]->SetParameter( 5, corrAnalysis::pi );
      leadPhiFit[i][j]->SetParameter( 3, 0.2 );
      leadPhiFit[i][j]->SetParameter( 6, 0.2 );
      leadPhiFit[i][j]->SetLineColor( i + 1 );
      
      leadPhiDifFit[i][j] = new TF1( dPhiLeadNameDif.c_str(), phiForm.c_str(), phiMin, phiMax );
      leadPhiDifFit[i][j]->FixParameter( 2, 0 );
      leadPhiDifFit[i][j]->SetParameter( 5, corrAnalysis::pi );
      leadPhiDifFit[i][j]->SetParameter( 3, 0.2 );
      leadPhiDifFit[i][j]->SetParameter( 6, 0.2 );
      leadPhiDifFit[i][j]->SetLineColor( i + 1 );
      
      subPhiFit[i][j] = new TF1( dPhiSubName.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiFit[i][j]->FixParameter( 2, 0 );
      subPhiFit[i][j]->SetParameter( 5, corrAnalysis::pi );
      subPhiFit[i][j]->SetParameter( 3, 0.2 );
      subPhiFit[i][j]->SetParameter( 6, 0.2 );
      subPhiFit[i][j]->SetLineColor( i + 1 );
      
      subPhiDifFit[i][j] = new TF1( dPhiSubNameDif.c_str(), phiForm.c_str(), phiMin, phiMax );
      subPhiDifFit[i][j]->FixParameter( 2, 0 );
      subPhiDifFit[i][j]->SetParameter( 5, corrAnalysis::pi );
      subPhiDifFit[i][j]->SetParameter( 3, 0.2 );
      subPhiDifFit[i][j]->SetParameter( 6, 0.2 );
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
      dPhiLead[i][j]->Fit( dPhiLeadName.c_str(), "R" );
      dPhiSub[i][j]->SetLineColor( i + 1 );
      dPhiSub[i][j]->Fit( dPhiSubName.c_str(), "R" );
      dPhiLeadNear[i][j]->SetLineColor( i + 1 );
      dPhiLeadNear[i][j]->Fit( dPhiLeadNameDif.c_str(), "R" );
      dPhiSubNear[i][j]->SetLineColor( i + 1 );
      dPhiSubNear[i][j]->Fit( dPhiSubNameDif.c_str(), "R" );
      dEtaLead[i][j]->SetLineColor( i + 1 );
      dEtaLead[i][j]->Fit( dEtaLeadName.c_str(), "R" );
      dEtaSub[i][j]->SetLineColor( i + 1 );
      dEtaSub[i][j]->Fit( dEtaSubName.c_str(), "R" );
      
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
  std::vector<std::vector<double> > leadPhiError( nFiles );
  std::vector<std::vector<double> > leadPhiDifYield( nFiles );
  std::vector<std::vector<double> > leadPhiDifError( nFiles );
  std::vector<std::vector<double> > leadEtaYield( nFiles );
  std::vector<std::vector<double> > leadEtaError( nFiles );
  std::vector<std::vector<double> > subPhiYield( nFiles );
  std::vector<std::vector<double> > subPhiError( nFiles );
  std::vector<std::vector<double> > subPhiDifYield( nFiles );
  std::vector<std::vector<double> > subPhiDifError( nFiles );
  std::vector<std::vector<double> > subEtaYield( nFiles );
  std::vector<std::vector<double> > subEtaError( nFiles );
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiYield[i].resize( nPtBins );
    leadPhiError[i].resize( nPtBins );
    leadPhiDifYield[i].resize( nPtBins );
    leadPhiDifError[i].resize( nPtBins );
    leadEtaYield[i].resize( nPtBins );
    leadEtaError[i].resize( nPtBins );
    subPhiYield[i].resize( nPtBins );
    subPhiError[i].resize( nPtBins );
    subPhiDifYield[i].resize( nPtBins );
    subPhiDifError[i].resize( nPtBins );
    subEtaYield[i].resize( nPtBins );
    subEtaError[i].resize( nPtBins );
    for ( int j = 0; j < nPtBins; ++j ) {
      leadPhiYield[i][j] = leadPhiFit[i][j]->GetParameter(1);
      leadPhiError[i][j] = leadPhiFit[i][j]->GetParError(1);
      leadPhiDifYield[i][j] = leadPhiDifFit[i][j]->GetParameter(1);
      leadPhiDifError[i][j] = leadPhiDifFit[i][j]->GetParError(1);
      leadEtaYield[i][j] = leadEtaFit[i][j]->GetParameter(1);
      leadEtaError[i][j] = leadEtaFit[i][j]->GetParError(1);
      subPhiYield[i][j] = subPhiFit[i][j]->GetParameter(1);
      subPhiError[i][j] = subPhiFit[i][j]->GetParError(1);
      subPhiDifYield[i][j] = subPhiDifFit[i][j]->GetParameter(1);
      subPhiDifError[i][j] = subPhiDifFit[i][j]->GetParError(1);
      subEtaYield[i][j] = subEtaFit[i][j]->GetParameter(1);
      subEtaError[i][j] = subEtaFit[i][j]->GetParError(1);
    }
  }
  
  double ptBins[5] = { 0.75, 1.5, 2.5, 3.5, 5 };
  
  std::vector<TGraphErrors*> leadPhiGraph( nFiles );
  std::vector<TGraphErrors*> leadPhiDifGraph( nFiles );
  std::vector<TGraphErrors*> leadEtaGraph( nFiles );
  std::vector<TGraphErrors*> subPhiGraph( nFiles );
  std::vector<TGraphErrors*> subPhiDifGraph( nFiles );
  std::vector<TGraphErrors*> subEtaGraph( nFiles );
  
  for ( int i = 0; i < nFiles; ++i ) {
    double leadPhiTmp[nPtBins];
    double leadPhiErr[nPtBins];
    double leadPhiDifTmp[nPtBins];
    double leadPhiDifErr[nPtBins];
    double leadEtaTmp[nPtBins];
    double leadEtaErr[nPtBins];
    double subPhiTmp[nPtBins];
    double subPhiErr[nPtBins];
    double subPhiDifTmp[nPtBins];
    double subPhiDifErr[nPtBins];
    double subEtaTmp[nPtBins];
    double subEtaErr[nPtBins];
    
    double errX[nPtBins];
    for ( int j = 0; j < nPtBins; ++j ) {
      errX[j] = 0;
      leadPhiTmp[j] = leadPhiYield[i][j];
      leadPhiErr[j] = leadPhiError[i][j];
      leadPhiDifTmp[j] = leadPhiDifYield[i][j];
      leadPhiDifErr[j] = leadPhiDifError[i][j];
      leadEtaTmp[j] = leadEtaYield[i][j];
      leadEtaErr[j] = leadEtaError[i][j];
      subPhiTmp[j] = subPhiYield[i][j];
      subPhiErr[j] = subPhiError[i][j];
      subPhiDifTmp[j] = subPhiDifYield[i][j];
      subPhiDifErr[j] = subPhiDifError[i][j];
      subEtaTmp[j] = subEtaYield[i][j];
      subEtaErr[j] = subEtaError[i][j];
    }
    
    leadPhiGraph[i] = new TGraphErrors(nPtBins, ptBins, leadPhiTmp, errX, leadPhiErr );
    
    leadPhiDifGraph[i] = new TGraphErrors(nPtBins, ptBins, leadPhiDifTmp, errX, leadPhiDifErr );
    
    leadEtaGraph[i] = new TGraphErrors(nPtBins, ptBins, leadEtaTmp, errX, leadEtaErr );
    
    subPhiGraph[i] = new TGraphErrors(nPtBins, ptBins, subPhiTmp, errX, subPhiErr );
    
    subPhiDifGraph[i] = new TGraphErrors(nPtBins, ptBins, subPhiDifTmp, errX, subPhiDifErr );
    
    subEtaGraph[i] = new TGraphErrors(nPtBins, ptBins, subEtaTmp, errX, subEtaErr );

    
  }
  
  TCanvas* c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiGraph[i]->SetLineColor(i+1);
    if ( i == 0)
      leadPhiGraph[i]->Draw();
    else
      leadPhiGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/leadphiyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiDifGraph[i]->SetLineColor(i+1);
    if ( i == 0)
      leadPhiDifGraph[i]->Draw();
    else
      leadPhiGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/leadphidifyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    leadEtaGraph[i]->SetLineColor(i+1);
    if ( i == 0)
      leadPhiGraph[i]->Draw();
    else
      leadPhiGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/leadetayield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiGraph[i]->SetLineColor(i+1);
    if ( i == 0)
      subPhiGraph[i]->Draw();
    else
      subPhiGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/subphiyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiDifGraph[i]->SetLineColor(i+1);
    if ( i == 0)
      subPhiDifGraph[i]->Draw();
    else
      subPhiDifGraph[i]->Draw("SAME");
  }
  c1->SaveAs("tmp/subphidifyield.pdf");
  c1 = new TCanvas;
  for ( int i = 0; i < nFiles; ++i ) {
    subEtaGraph[i]->SetLineColor(i+1);
    if ( i == 0) {
      subEtaGraph[i]->GetYaxis()->SetRangeUser( 0, 13 );
      subEtaGraph[i]->Draw();
    }
    else {
      subEtaGraph[i]->Draw("SAME");
    }
  }
  c1->SaveAs("tmp/subetayield.pdf");

  return 0;
}

