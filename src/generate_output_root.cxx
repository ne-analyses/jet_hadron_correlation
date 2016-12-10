// macro to produce output of dijet hadron
// correlations with event mixing, efficiency
// correction and pt bin dependence.
// Nick Elsey 10/04/2016

// All reader and histogram settings
// Are located in corrParameters.hh
#include "src/corrParameters.hh"

// The majority of the jetfinding
// And correlation code is located in
// corrFunctions.hh
//#include "src/corrFunctions.hh"

// ROOT is used for histograms and
// As a base for the TStarJetPico library
// ROOT Headers
//#include "TH1.h"
//#include "TH2.h"
//#include "TH3.h"
//#include "TF1.h"
//#include "TF2.h"
//#include "TProfile.h"
//#include "TProfile2D.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TLegend.h"
//#include "TObjArray.h"
//#include "TString.h"
//#include "TFile.h"
//#include "TLorentzVector.h"
//#include "TClonesArray.h"
//#include "TChain.h"
//#include "TBranch.h"
//#include "TMath.h"
//#include "TRandom.h"
//#include "TRandom3.h"
//#include "TCanvas.h"
//#include "TStopwatch.h"
//#include "TSystem.h"
//#include "TStyle.h"

// Make use of std::vector,
// std::string, IO and algorithm
// STL Headers
//#include <stdio.h>
//#include <stdlib.h>
//#include <sstream>
//#include <iostream>
//#include <fstream>
//#include <cstring>
//#include <algorithm>
//#include <cstring>
//#include <vector>
//#include <string>
//#include <limits.h>
//#include <unistd.h>

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

int generate_output_root() {
  
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
  defaultCorrNames.resize(2);
  defaultCorrNames[0] = "dijet";
  defaultCorrNames[1] = "ppdijet";
  //defaultCorrNames[2] = "15 < Jet < 20";
  //defaultCorrNames[3] = "Jet > 20";
  
  
  // First check to make sure we're located properly
  //std::string currentDirectory = jetHadron::getPWD( );
  
  // If we arent in the analysis directory, exit
//  if ( !(jetHadron::HasEnding ( currentDirectory, "jet_hadron_corr" ) || jetHadron::HasEnding ( currentDirectory, "jet_hadron_correlation" )) ) {
//    std::cerr << "Error: Need to be in jet_hadron_corr directory" << std::endl;
//    return -1;
//  }
  
  
  // files and naming
  TFile** corrFiles;
  TFile** mixFiles;
  std::vector<std::string> analysisNames;
  
  switch ( 1 ) {
    case 1: { // Default case
      //__OUT( "Using Default Settings" )
      corrFiles = new TFile*[2];
      mixFiles = new TFile*[2];
      
      // default files
      corrFiles[0] = new TFile( "out/tmp/dijet_corr24.root", "READ" );
      corrFiles[1] = new TFile( "out/tmp/ppdijet_corr24.root", "READ" );
      //corrFiles[2] = new TFile( "out/tmp/jet15_corr.root", "READ" );
      //corrFiles[3] = new TFile( "out/tmp/jet20_corr.root", "READ" );
      mixFiles[0] = new TFile( "out/tmp/dijet_mix24.root", "READ" );
      mixFiles[1] = new TFile( "out/tmp/ppdijet_mix24.root", "READ" );
      //mixFiles[2] = new TFile( "out/tmp/jet15_mix.root", "READ" );
      //mixFiles[3] = new TFile( "out/tmp/jet20_mix.root", "READ" );
      analysisNames = defaultCorrNames;
      
      break;
    }
    default: {
      if ( (argc-1)%3 != 0 ) {
        //__ERR("Need correlation file, mixing file, and analysis name for each entry")
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
  
  const int nFiles = analysisNames.size();
  
  // Load in the histograms and get the pt spectra
  TH2D* nEvents[ nFiles ];
  TH1D* hVz[ nFiles ];
  TH3D* corrHist[ nFiles ];
  TH3D* mixHist[ nFiles ];
  TH3D* corrCentVz[nFiles][jetHadron::binsCentrality][jetHadron::binsVz];
  TH3D* subCentVz[nFiles][jetHadron::binsCentrality][jetHadron::binsVz];
  TH3D* mixCentVz[nFiles][jetHadron::binsCentrality][jetHadron::binsVz];
  TH3D* mixSubCentVz[nFiles][jetHadron::binsCentrality][jetHadron::binsVz];
  TH1D* recombinedPtLead[nFiles];
  TH1D* recombinedPtSub[nFiles];
  for ( int i = 0; i < nFiles; ++i ) {
    
    std::string ptLeadName = analysisNames[i] + "_pt_lead";
    std::string ptSubName = analysisNames[i] + "_pt_sub";
    
    recombinedPtLead[i] = new TH1D( ptLeadName.c_str(), "p_{T} Spectrum Trigger Jet", jetHadron::binsPt, jetHadron::ptLowEdge, jetHadron::ptHighEdge );
    recombinedPtSub[i] = new TH1D( ptSubName.c_str(), "p_{T} Spectrum Recoil Jet", jetHadron::binsPt, jetHadron::ptLowEdge, jetHadron::ptHighEdge );
    
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
    for ( int j = 0; j < jetHadron::binsCentrality; ++j )
      for ( int k = 0; k < jetHadron::binsVz; ++k ) {
        // make the initial name
        TString corrDifInitName = "lead_cent_"; corrDifInitName += j;
        TString subDifInitName = "sub_cent_"; subDifInitName += j;
        TString mixDifInitName = "mix_lead_cent_"; mixDifInitName += j;
        TString mixSubDifInitName = "mix_sub_cent_";
        mixSubDifInitName += j;
        
        corrDifInitName += "_vz_"; corrDifInitName += k;
        subDifInitName += "_vz_"; subDifInitName += k;
        mixDifInitName += "_vz_"; mixDifInitName += k;
        mixSubDifInitName += "_vz_"; mixSubDifInitName += k;
        
        
        // make the new histogram name
        TString corrDifBaseName = "corr_file_"; corrDifBaseName += i;
        TString subDifBaseName = "sub_file_"; subDifBaseName += i;
        TString mixDifBaseName = "mix_file_"; mixDifBaseName += i;
        TString mixSubDifBaseName = "mix_file_"; mixSubDifBaseName += i;
        
        corrDifBaseName += "_cent_"; corrDifBaseName += j;
        corrDifBaseName += "_vz_"; corrDifBaseName += k;
        
        subDifBaseName += "_cent_"; subDifBaseName += j;
        subDifBaseName += "_vz_"; subDifBaseName += k;
        
        mixDifBaseName += "_cent_"; mixDifBaseName += j;
        mixDifBaseName += "_vz_"; mixDifBaseName += k;
        
        mixSubDifBaseName += "_cent_"; mixSubDifBaseName += j;
        mixSubDifBaseName += "_vz_"; mixSubDifBaseName += k;
        
        // get the histograms
        corrCentVz[i][j][k] = (TH3D*) corrFiles[i]->Get( corrDifInitName );
        corrCentVz[i][j][k]->SetName( corrDifBaseName );
        subCentVz[i][j][k] = (TH3D*) corrFiles[i]->Get( subDifInitName );
        subCentVz[i][j][k]->SetName( subDifBaseName );
        mixCentVz[i][j][k] = (TH3D*) mixFiles[i]->Get( mixDifInitName );
        mixCentVz[i][j][k]->SetName( mixDifBaseName );
        mixSubCentVz[i][j][k] = (TH3D*) mixFiles[i]->Get( mixSubDifInitName );
        mixSubCentVz[i][j][k]->SetName( mixSubDifBaseName );
        
        // now we can get the pt spectrum as well
        recombinedPtLead[i]->Add( (TH1D*) corrCentVz[i][j][k]->Project3D("Z") );
        recombinedPtSub[i]->Add( (TH1D*) subCentVz[i][j][k]->Project3D("Z") );
      }
  }
  
  std::cout<<"loaded all histograms"<<std::endl;
  
  
  // setup for 2d projections along pt axis
  TH2D* corrCentVzPt[nFiles][jetHadron::binsCentrality][jetHadron::binsVz][nPtBins];
  TH2D* subCentVzPt[nFiles][jetHadron::binsCentrality][jetHadron::binsVz][nPtBins];
  TH2D* mixCentVzPt[nFiles][jetHadron::binsCentrality][jetHadron::binsVz][nPtBins];
  TH2D* mixSubCentVzPt[nFiles][jetHadron::binsCentrality][jetHadron::binsVz][nPtBins];
  
  // now get the pt projections
  for ( int i = 0; i < nFiles; ++i ) {
    for ( int j = 0; j < jetHadron::binsCentrality; ++j ) {
      for ( int k = 0; k < jetHadron::binsVz; ++ k ) {
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
    for ( int j = 0; j < jetHadron::binsCentrality; ++j ) {
      for ( int k = 0; k < jetHadron::binsVz; ++k ) {
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
  TH2D* recombinedCorr[nFiles][nPtBins];
  TH2D* recombinedSub[nFiles][nPtBins];
  TH2D* recombinedPre[nFiles][nPtBins];
  TH2D* recombinedSubPre[nFiles][nPtBins];
  
  for (int i = 0; i < nFiles; ++ i ) {
    

    
    for ( int l = 0; l < nPtBins; ++l ) {
      
      TString corrName = analysisNames[i] + " " + ptBinString[l];
      TString subName = analysisNames[i] + "_sub " + ptBinString[l];
      TString preName = "pre_" + analysisNames[i] + " " + ptBinString[l];
      TString subPreName = "pre_" + analysisNames[i] + "_sub " + ptBinString[l];
      
      recombinedCorr[i][l] = new TH2D( corrName, corrName, jetHadron::binsEta, jetHadron::dEtaLowEdge, jetHadron::dEtaHighEdge, jetHadron::binsPhi, jetHadron::phiLowEdge, jetHadron::phiHighEdge );
      
      recombinedSub[i][l] = new TH2D( subName, subName, jetHadron::binsEta, jetHadron::dEtaLowEdge, jetHadron::dEtaHighEdge, jetHadron::binsPhi, jetHadron::phiLowEdge, jetHadron::phiHighEdge );
      
      recombinedPre[i][l] = new TH2D( preName, preName, jetHadron::binsEta, jetHadron::dEtaLowEdge, jetHadron::dEtaHighEdge, jetHadron::binsPhi, jetHadron::phiLowEdge, jetHadron::phiHighEdge );
      
      recombinedSubPre[i][l] = new TH2D( subPreName, subPreName, jetHadron::binsEta, jetHadron::dEtaLowEdge, jetHadron::dEtaHighEdge, jetHadron::binsPhi, jetHadron::phiLowEdge, jetHadron::phiHighEdge );
      
      if ( l <= 2 ) {
      
        for ( int j = 0; j < jetHadron::binsCentrality; ++j ) {
          for ( int k = 0; k < jetHadron::binsVz; ++k ) {
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
        for ( int j = 0; j < jetHadron::binsCentrality; ++j ) {
          for ( int k = 0; k < jetHadron::binsVz; ++k ) {
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
  double phiMin = -jetHadron::pi/2.0;
  double phiMinClose = -0.6;
  double phiMaxClose = 0.6;
  double phiMaxFar = 3.0*jetHadron::pi/2.0;
  double phiMax = 3.0*jetHadron::pi/2.0;
  double phiDifMax = jetHadron::pi/2.0;
  
  
  // going to get the 1D projections
  TH1D* dPhiLead[nFiles][nPtBins];
  TH1D* dPhiLeadNear[nFiles][nPtBins];
  TH1D* dPhiLeadFar[nFiles][nPtBins];
  TH1D* dEtaLead[nFiles][nPtBins];
  TH1D* dPhiSub[nFiles][nPtBins];
  TH1D* dPhiSubNear[nFiles][nPtBins];
  TH1D* dPhiSubFar[nFiles][nPtBins];
  TH1D* dEtaSub[nFiles][nPtBins];
  
  TCanvas* c1;
  for ( int i = 0; i < nFiles; ++i ) {
    
    for ( int j = 0; j < nPtBins; ++j ) {
      // first restrict the eta range
      recombinedCorr[i][j]->GetXaxis()->SetRangeUser( etaMin, etaMax  );
      recombinedSub[i][j]->GetXaxis()->SetRangeUser( etaMin, etaMax );
      
      // save the 2D histograms
      TString leadOutName = "tmp/lead2d_" + analysisNames[i] +"_pt_" + j + ".pdf";
      TString subOutName = "tmp/sub2d_" + analysisNames[i] +"_pt_"+ j + ".pdf";
      
      c1 = new TCanvas();
      recombinedCorr[i][j]->DrawCopy("surf1");
      c1->SaveAs( leadOutName );
      c1 = new TCanvas();
      recombinedSub[i][j]->DrawCopy("surf1");
      c1->SaveAs( subOutName );
      
      
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
      TString leadPhiName = "tmp/lead_phi_near_far_"+analysisNames[i]+"_pt_"+ j +".pdf";
      TString subPhiName = "tmp/sub_phi_near_far_"+analysisNames[i]+"_pt_"+ j +".pdf";
      
      dPhiLeadNear[i][j]->SetLineColor(kBlack);
      dPhiLeadNear[i][j]->SetMarkerStyle(29);
      dPhiLeadNear[i][j]->SetMarkerSize(3);
      dPhiLeadNear[i][j]->SetMarkerColor(kBlack);
      dPhiLeadNear[i][j]->DrawCopy();
      dPhiLeadFar[i][j]->SetLineColor(kRed);
      dPhiLeadFar[i][j]->SetMarkerStyle(29);
      dPhiLeadFar[i][j]->SetMarkerSize(3);
      dPhiLeadFar[i][j]->SetMarkerColor(kRed);
      dPhiLeadFar[i][j]->DrawCopy("SAME");
      c1 = new TCanvas();
      c1->SaveAs( leadPhiName );
      dPhiSubNear[i][j]->SetLineColor(kBlack);
      dPhiSubNear[i][j]->SetMarkerStyle(29);
      dPhiSubNear[i][j]->SetMarkerSize(3);
      dPhiSubNear[i][j]->SetMarkerColor(kBlack);
      dPhiSubNear[i][j]->DrawCopy();
      dPhiSubFar[i][j]->SetLineColor(kRed);
      dPhiSubFar[i][j]->SetMarkerStyle(29);
      dPhiSubFar[i][j]->SetMarkerSize(3);
      dPhiSubFar[i][j]->SetMarkerColor(kRed);
      dPhiSubFar[i][j]->DrawCopy("SAME");
      c1 = new TCanvas();
      c1->SaveAs( subPhiName );
      
      // Now do the subtraction
      dPhiLeadNear[i][j]->Add( dPhiLeadFar[i][j], -1 );
      dPhiSubNear[i][j]->Add( dPhiSubFar[i][j], -1 );
    }

  
  // Now to do some fitting and subtract the background
  // define the fits
  // ---------------
  TString phiForm = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)+[4]*exp(-0.5*((x-[5])/[6])**2)";
  TString etaForm = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)";
  TString phiDifForm = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)";

  // do a first, temporary fit to remove background
  for ( int i = 0; i < nFiles; ++i ) {
    for ( int j = 0; j < nPtBins; ++j ) {
      TString dPhiLeadName = "tmp_fit_lead_phi_" + i + j;
      TString dPhiSubName = "tmp_fit_sub_phi_" + i + j;
      TString dPhiLeadNameDif = "tmp_fit_lead_phi_dif_" + i + j;
      TString dPhiSubNameDif = "tmp_fit_sub_phi_dif_" + i + j;
      TString dEtaLeadName = "tmp_fit_lead_eta_" + i + j;
      TString dEtaSubName = "tmp_fit_sub_eta_" + i + j;
      
      TF1* leadPhiInitFit = new TF1( dPhiLeadName, phiForm, phiMin, phiDifMax );
      leadPhiInitFit->FixParameter( 2, 0 );
      leadPhiInitFit->FixParameter( 5, jetHadron::pi );
      leadPhiInitFit->SetParameter( 3, 0.2 );
      leadPhiInitFit->SetParameter( 6, 0.2 );
      
      TF1* leadPhiDifInitFit = new TF1( dPhiLeadNameDif, phiDifForm, phiMin, phiMax );
      leadPhiDifInitFit->FixParameter( 2, 0 );
      leadPhiDifInitFit->SetParameter( 3, 0.2 );
      
      TF1* subPhiInitFit = new TF1( dPhiSubName, phiForm, phiMin, phiMax );
      subPhiInitFit->FixParameter( 2, 0 );
      subPhiInitFit->FixParameter( 5, jetHadron::pi );
      subPhiInitFit->SetParameter( 3, 0.2 );
      subPhiInitFit->SetParameter( 6, 0.2 );
      
      TF1* subPhiDifInitFit = new TF1( dPhiSubNameDif, phiDifForm, phiMin, phiDifMax );
      subPhiDifInitFit->FixParameter( 2, 0 );
      subPhiDifInitFit->SetParameter( 3, 0.2 );
      
      TF1* leadEtaInitFit = new TF1( dEtaLeadName, etaForm, etaMin, etaMax );
      leadEtaInitFit->FixParameter( 2, 0 );
      leadEtaInitFit->SetParameter( 3, 0.2 );
      
      TF1* subEtaInitFit = new TF1( dEtaSubName, etaForm, etaMin, etaMax );
      subEtaInitFit->FixParameter( 2, 0 );
      subEtaInitFit->SetParameter( 3, 0.2 );
      
      dPhiLead[i][j]->Fit( dPhiLeadName, "RMI" );
      dPhiSub[i][j]->Fit( dPhiSubName, "RMI" );
      dPhiLeadNear[i][j]->Fit( dPhiLeadNameDif, "RMI" );
      dPhiSubNear[i][j]->Fit( dPhiSubNameDif, "RMI" );
      dEtaLead[i][j]->Fit( dEtaLeadName, "RMI" );
      dEtaSub[i][j]->Fit( dEtaSubName, "RMI" );
      
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
  TF1* leadPhiFit[nFiles][nPtBins];
  TF1* leadPhiDifFit[nFiles][nPtBins];
  TF1* leadEtaFit[nFiles][nPtBins];
  TF1* subPhiFit[nFiles][nPtBins];
  TF1* subPhiDifFit[nFiles][nPtBins];
  TF1* subEtaFit[nFiles][nPtBins];
  
  for ( int i = 0; i < nFiles; ++i ) {
    
    for ( int j = 0; j < nPtBins; ++j ) {
      TString dPhiLeadName = "fit_lead_phi_" + i + j;
      TString dPhiSubName = "fit_sub_phi_" + i + j;
      TString dPhiLeadNameDif = "fit_lead_phi_dif_" + i + j;
      TString dPhiSubNameDif = "fit_sub_phi_dif_" + i + j;
      TString dEtaLeadName = "fit_lead_eta_" + i + j;
      TString dEtaSubName = "fit_sub_eta_" + i + j;
    
      leadPhiFit[i][j] = new TF1( dPhiLeadName, phiForm, phiMin, phiMax );
      leadPhiFit[i][j]->FixParameter( 2, 0 );
      leadPhiFit[i][j]->FixParameter( 5, jetHadron::pi );
      leadPhiFit[i][j]->SetParameter( 3, 0.2 );
      leadPhiFit[i][j]->SetParameter( 6, 0.2 );
      leadPhiFit[i][j]->SetLineColor( i + 1 );
      
      leadPhiDifFit[i][j] = new TF1( dPhiLeadNameDif, phiDifForm, phiMin, phiDifMax );
      leadPhiDifFit[i][j]->FixParameter( 2, 0 );
      leadPhiDifFit[i][j]->SetParameter( 3, 0.2 );
      leadPhiDifFit[i][j]->SetLineColor( i + 1 );
      
      subPhiFit[i][j] = new TF1( dPhiSubName, phiForm, phiMin, phiMax );
      subPhiFit[i][j]->FixParameter( 2, 0 );
      subPhiFit[i][j]->FixParameter( 5, jetHadron::pi );
      subPhiFit[i][j]->SetParameter( 3, 0.2 );
      subPhiFit[i][j]->SetParameter( 6, 0.2 );
      subPhiFit[i][j]->SetLineColor( i + 1 );
      
      subPhiDifFit[i][j] = new TF1( dPhiSubNameDif, phiDifForm, phiMin, phiDifMax );
      subPhiDifFit[i][j]->FixParameter( 2, 0 );
      subPhiDifFit[i][j]->SetParameter( 3, 0.2 );
      subPhiDifFit[i][j]->SetLineColor( i + 1 );
      
      leadEtaFit[i][j] = new TF1( dEtaLeadName, etaForm, etaMin, etaMax );
      leadEtaFit[i][j]->FixParameter( 2, 0 );
      leadEtaFit[i][j]->SetParameter( 3, 0.2 );
      leadEtaFit[i][j]->SetLineColor( i + 1 );
      
      subEtaFit[i][j] = new TF1( dEtaSubName, etaForm, etaMin, etaMax );
      subEtaFit[i][j]->FixParameter( 2, 0 );
      subEtaFit[i][j]->SetParameter( 3, 0.2 );
      subEtaFit[i][j]->SetLineColor( i + 1 );

      // Now set same colors and fit
      dPhiLead[i][j]->SetLineColor( i + 1 );
      dPhiLead[i][j]->Fit( dPhiLeadName, "RMI" );
      dPhiSub[i][j]->SetLineColor( i + 1 );
      dPhiSub[i][j]->Fit( dPhiSubName, "RMI" );
      dPhiLeadNear[i][j]->SetLineColor( i + 1 );
      dPhiLeadNear[i][j]->Fit( dPhiLeadNameDif, "RMI" );
      dPhiSubNear[i][j]->SetLineColor( i + 1 );
      dPhiSubNear[i][j]->Fit( dPhiSubNameDif, "RMI" );
      dEtaLead[i][j]->SetLineColor( i + 1 );
      dEtaLead[i][j]->Fit( dEtaLeadName, "RMI" );
      dEtaSub[i][j]->SetLineColor( i + 1 );
      dEtaSub[i][j]->Fit( dEtaSubName, "RMI" );
      
    }
  }
  
  // Now start making output
  TString outBase = "tmp/";
  TString leadPhiOutBase = outBase + "leadphi_pt";
  TString leadPhiDifOutBase = outBase + "leadphidif_pt";
  TString leadEtaOutBase = outBase + "leadeta_pt";
  TString subPhiOutBase = outBase + "subphi_pt";
  TString subPhiDifOutBase = outBase + "subphidif_pt";
  TString subEtaOutBase = outBase + "subeta_pt";
  TString outExt = ".pdf";
  
  for ( int i = 0; i < nPtBins; ++i ) {
    c1 = new TCanvas();
    TString leadPhiOut = leadPhiOutBase;
    leadPhiOut += i;
    leadPhiOut += outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      dPhiLead[j][i]->SetLineColor(j+1);
      dPhiLead[j][i]->SetMarkerStyle(29);
      dPhiLead[j][i]->SetMarkerSize(3);
      dPhiLead[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        TString outTitle = "Trigger Jet #Delta#phi " + ptBinString[i];
        dPhiLead[j][i]->SetTitle( outTitle );
        dPhiLead[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiLead[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiLead[j][i]->DrawCopy();
      }
      else {
        dPhiLead[j][i]->DrawCopy("same");
      }
    }
    c1->SaveAs( leadPhiOut );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    c1 = new TCanvas();
    TString leadPhiDifOut = leadPhiDifOutBase;
    leadPhiDifOut += i;
    leadPhiDifOut += outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      dPhiLeadNear[j][i]->SetLineColor(j+1);
      dPhiLeadNear[j][i]->SetMarkerStyle(29);
      dPhiLeadNear[j][i]->SetMarkerSize(3);
      dPhiLeadNear[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        TString outTitle = "Trigger Jet #Delta#eta subtracted #Delta#phi " + ptBinString[i];
        dPhiLeadNear[j][i]->SetTitle( outTitle );
        dPhiLeadNear[j][i]->GetXaxis()->SetRangeUser( -jetHadron::pi/2.0, jetHadron::pi/2.0 );
        dPhiLeadNear[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiLeadNear[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiLeadNear[j][i]->DrawCopy();
      }
      else {
        dPhiLeadNear[j][i]->DrawCopy("same");
      }
    }
    c1->SaveAs( leadPhiDifOut );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    c1 = new TCanvas();
    
    TString leadEtaOut = leadEtaOutBase;
    leadEtaOut += i;
    leadEtaOut += outExt;
    for ( int j = 0; j < nFiles; ++ j ) {
      dEtaLead[j][i]->SetLineColor(j+1);
      dEtaLead[j][i]->SetMarkerStyle(29);
      dEtaLead[j][i]->SetMarkerSize(3);
      dEtaLead[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        TString outTitle = "Trigger Jet #Delta#eta " + ptBinString[i];
        dEtaLead[j][i]->SetTitle( outTitle );
        dEtaLead[j][i]->GetXaxis()->SetTitle("#Delta#eta");
        dEtaLead[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#eta");
        dEtaLead[j][i]->DrawCopy();
      }
      else {
        dEtaLead[j][i]->DrawCopy("same");
      }
    }
    c1->SaveAs( leadEtaOut );
  }
  
  for ( int i = 0; i < nPtBins; ++i ) {
    c1 = new TCanvas();
    
    TString subPhiOut = subPhiOutBase;
    subPhiOut += i;
    subPhiOut += outExt;
    for ( int j = 0; j < nFiles; ++j ) {
      dPhiSub[j][i]->SetLineColor(j+1);
      dPhiSub[j][i]->SetMarkerStyle(29);
      dPhiSub[j][i]->SetMarkerSize(3);
      dPhiSub[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        TString outTitle = "Recoil Jet #Delta#phi " + ptBinString[i];
        dPhiSub[j][i]->SetTitle( outTitle );
        dPhiSub[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiSub[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiSub[j][i]->DrawCopy();
      }
      else {
        dPhiSub[j][i]->DrawCopy("same");
      }
    }
    c1->SaveAs( subPhiOut );
  }

  for ( int i = 0; i < nPtBins; ++i ) {
    c1 = new TCanvas();
    
    TString subPhiDifOut = subPhiDifOutBase;
    subPhiDifOut += i;
    subPhiDifOut += outExt;
    for ( int j = 0; j < nFiles; ++j ) {
      dPhiSubNear[j][i]->SetLineColor(j+1);
      dPhiSubNear[j][i]->SetMarkerStyle(29);
      dPhiSubNear[j][i]->SetMarkerSize(3);
      dPhiSubNear[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        TString outTitle = "Recoil Jet #Delta#eta subtracted #Delta#phi " + ptBinString[i];
        dPhiSubNear[j][i]->SetTitle( outTitle );
        dPhiSubNear[j][i]->GetXaxis()->SetRangeUser( -jetHadron::pi/2.0, jetHadron::pi/2.0 );
        dPhiSubNear[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        dPhiSubNear[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#phi");
        dPhiSubNear[j][i]->DrawCopy();
      }
      else {
        dPhiSubNear[j][i]->DrawCopy("same");
      }
    }
    c1->SaveAs( subPhiDifOut );
  }

  for ( int i = 0; i < nPtBins; ++i ) {
    c1 = new TCanvas();
    
    TString subEtaOut = subEtaOutBase;
    subEtaOut += i;
    subEtaOut += outExt;
    for ( int j = 0; j < nFiles; ++j ) {
      dEtaSub[j][i]->SetLineColor(j+1);
      dEtaSub[j][i]->SetMarkerStyle(29);
      dEtaSub[j][i]->SetMarkerSize(3);
      dEtaSub[j][i]->SetMarkerColor(j+1);
      if ( j == 0 ) {
        TString outTitle = "Recoil Jet #Delta#eta " + ptBinString[i];
        dEtaSub[j][i]->SetTitle( outTitle );
        dEtaSub[j][i]->GetXaxis()->SetTitle("#Delta#eta");
        dEtaSub[j][i]->GetYaxis()->SetTitle("1/N_{dijet}dN/d#eta");
        dEtaSub[j][i]->DrawCopy();
      }
      else {
        dEtaSub[j][i]->DrawCopy("same");
      }
    }
    c1->SaveAs( subEtaOut );
  }
  
  // now to get yields
  double leadPhiYield[nFiles][nPtBins];
  double leadPhiWidth[nFiles][nPtBins];
  double leadPhiError[nFiles][nPtBins];
  double leadPhiWidthError[nFiles][nPtBins];
  double leadPhiDifYield[nFiles][nPtBins];
  double leadPhiDifWidth[nFiles][nPtBins];
  double leadPhiDifError[nFiles][nPtBins];
  double leadPhiDifWidthError[nFiles][nPtBins];
  double leadEtaYield[nFiles][nPtBins];
  double leadEtaWidth[nFiles][nPtBins];
  double leadEtaError[nFiles][nPtBins];
  double leadEtaWidthError[nFiles][nPtBins];
  double subPhiYield[nFiles][nPtBins];
  double subPhiWidth[nFiles][nPtBins];
  double subPhiError[nFiles][nPtBins];
  double subPhiWidthError[nFiles][nPtBins];
  double subPhiDifYield[nFiles][nPtBins];
  double subPhiDifWidth[nFiles][nPtBins];
  double subPhiDifError[nFiles][nPtBins];
  double subPhiDifWidthError[nFiles][nPtBins];
  double subEtaYield[nFiles][nPtBins];
  double subEtaWidth[nFiles][nPtBins];
  double subEtaError[nFiles][nPtBins];
  double subEtaWidthError[nFiles][nPtBins];
  for ( int i = 0; i < nFiles; ++i ) {
    for ( int j = 0; j < nPtBins; ++j ) {
      leadPhiYield[i][j] = leadPhiFit[i][j]->GetParameter(1)*sqrt(2*jetHadron::pi)*fabs(leadPhiFit[i][j]->GetParameter(3))/ptBinWidth[j];
      leadPhiError[i][j] = leadPhiFit[i][j]->GetParError(1);
      leadPhiWidth[i][j] = fabs(leadPhiFit[i][j]->GetParameter(3));
      leadPhiWidthError[i][j] = leadPhiFit[i][j]->GetParError(3);
      leadPhiDifYield[i][j] = leadPhiDifFit[i][j]->GetParameter(1)*sqrt(2*jetHadron::pi)*fabs(leadPhiDifFit[i][j]->GetParameter(3))/ptBinWidth[j];
      leadPhiDifError[i][j] = leadPhiDifFit[i][j]->GetParError(1);
      leadPhiDifWidth[i][j] = fabs(leadPhiDifFit[i][j]->GetParameter(3));
      leadPhiDifWidthError[i][j] = leadPhiDifFit[i][j]->GetParError(3);
      leadEtaYield[i][j] = leadEtaFit[i][j]->GetParameter(1)*sqrt(2*jetHadron::pi)*fabs(leadEtaFit[i][j]->GetParameter(3))/ptBinWidth[j];
      leadEtaError[i][j] = leadEtaFit[i][j]->GetParError(1);
      leadEtaWidth[i][j] = fabs(leadEtaFit[i][j]->GetParameter(3));
      leadEtaWidthError[i][j] = leadEtaFit[i][j]->GetParError(3);
      subPhiYield[i][j] = subPhiFit[i][j]->GetParameter(1)*sqrt(2*jetHadron::pi)*fabs(subPhiFit[i][j]->GetParameter(3))/ptBinWidth[j];
      subPhiError[i][j] = subPhiFit[i][j]->GetParError(1);
      subPhiWidth[i][j] = fabs(subPhiFit[i][j]->GetParameter(3));
      subPhiWidthError[i][j] = subPhiFit[i][j]->GetParError(3);
      subPhiDifYield[i][j] = subPhiDifFit[i][j]->GetParameter(1)*sqrt(2*jetHadron::pi)*fabs(subPhiDifFit[i][j]->GetParameter(3))/ptBinWidth[j];
      subPhiDifError[i][j] = subPhiDifFit[i][j]->GetParError(1);
      subPhiDifWidth[i][j] = fabs(subPhiDifFit[i][j]->GetParameter(3));
      subPhiDifWidthError[i][j] = subPhiDifFit[i][j]->GetParError(3);
      subEtaYield[i][j] = subEtaFit[i][j]->GetParameter(1)*sqrt(2*jetHadron::pi)*fabs(subEtaFit[i][j]->GetParameter(3))/ptBinWidth[j];
      subEtaError[i][j] = subEtaFit[i][j]->GetParError(1);
      subEtaWidth[i][j] = fabs(subEtaFit[i][j]->GetParameter(3));
      subEtaWidthError[i][j] = subEtaFit[i][j]->GetParError(3);
    }
  }
  
  // finding the weighted center of the pt bins for future graphs
  double ptBinCenter[nFiles][nPtBins];
  for ( int i = 0; i < nFiles; ++i ) {
    
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
  
  TGraphErrors* leadPhiGraph[nFiles];
  TGraphErrors* leadPhiWidthGraph[nFiles];
  TGraphErrors* leadPhiDifGraph[nFiles];
  TGraphErrors* leadPhiDifWidthGraph[nFiles];
  TGraphErrors* leadEtaGraph[nFiles];
  TGraphErrors* leadEtaWidthGraph[nFiles];
  TGraphErrors* subPhiGraph[nFiles];
  TGraphErrors* subPhiWidthGraph[nFiles];
  TGraphErrors* subPhiDifGraph[nFiles];
  TGraphErrors* subPhiDifWidthGraph[nFiles];
  TGraphErrors* subEtaGraph[nFiles];
  TGraphErrors* subEtaWidthGraph[nFiles];
  
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
      tmpPtBin[j] = ptBinCenter[i][j+startPtBin];
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
  
  c1 = new TCanvas();
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
      leadPhiGraph[i]->Draw("P");
    else
      leadPhiGraph[i]->Draw("P");
  }
  c1->SaveAs("tmp/leadphiyield.pdf");
  c1 = new TCanvas();
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
      leadPhiDifGraph[i]->Draw("P");
    else
      leadPhiDifGraph[i]->Draw("P");
  }
  c1->SaveAs("tmp/leadphidifyield.pdf");
  c1 = new TCanvas();
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
      leadEtaGraph[i]->Draw("P");
    else
      leadEtaGraph[i]->Draw("P");
  }
  c1->SaveAs("tmp/leadetayield.pdf");
  c1 = new TCanvas();
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
      subPhiGraph[i]->Draw("P");
    else
      subPhiGraph[i]->Draw("P");
  }
  c1->SaveAs("tmp/subphiyield.pdf");
  c1 = new TCanvas();
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
      subPhiDifGraph[i]->Draw("P");
    else
      subPhiDifGraph[i]->Draw("P");
  }
  c1->SaveAs("tmp/subphidifyield.pdf");
  c1 = new TCanvas();
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
      subEtaGraph[i]->Draw("P");
    }
    else {
      subEtaGraph[i]->Draw("P");
    }
  }
  c1->SaveAs("tmp/subetayield.pdf");
  
  // now print out widths
  c1 = new TCanvas();
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiWidthGraph[i]->SetLineColor(i+1);
    leadPhiWidthGraph[i]->SetMarkerStyle(29);
    leadPhiWidthGraph[i]->SetMarkerSize(3);
    leadPhiWidthGraph[i]->SetMarkerColor(i+1);
    leadPhiWidthGraph[i]->SetTitle("Trigger Jet #Delta#phi - Gaussian Widths");
    leadPhiWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    leadPhiWidthGraph[i]->GetYaxis()->SetTitle("Width");
    leadPhiWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      leadPhiWidthGraph[i]->Draw("P");
    else
      leadPhiWidthGraph[i]->Draw("P");
  }
  c1->SaveAs("tmp/leadphiwidth.pdf");
  c1 = new TCanvas();
  for ( int i = 0; i < nFiles; ++i ) {
    leadPhiDifWidthGraph[i]->SetLineColor(i+1);
    leadPhiDifWidthGraph[i]->SetMarkerStyle(29);
    leadPhiDifWidthGraph[i]->SetMarkerSize(3);
    leadPhiDifWidthGraph[i]->SetMarkerColor(i+1);
    leadPhiDifWidthGraph[i]->SetTitle("Trigger Jet #Delta#eta Subtracted #Delta#phi - Gaussian Widths");
    leadPhiDifWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    leadPhiDifWidthGraph[i]->GetYaxis()->SetTitle("Width");
    leadPhiDifWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      leadPhiDifWidthGraph[i]->Draw("P");
    else
      leadPhiDifWidthGraph[i]->Draw("P");
  }
  c1->SaveAs("tmp/leadphidifwidth.pdf");
  c1 = new TCanvas();
  for ( int i = 0; i < nFiles; ++i ) {
    leadEtaWidthGraph[i]->SetLineColor(i+1);
    leadEtaWidthGraph[i]->SetMarkerStyle(29);
    leadEtaWidthGraph[i]->SetMarkerSize(3);
    leadEtaWidthGraph[i]->SetMarkerColor(i+1);
    leadEtaWidthGraph[i]->SetTitle("Trigger Jet #Delta#eta - Gaussian Widths");
    leadEtaWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    leadEtaWidthGraph[i]->GetYaxis()->SetTitle("Width");
    leadEtaWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      leadEtaWidthGraph[i]->Draw("P");
    else
      leadEtaWidthGraph[i]->Draw("P");
  }
  c1->SaveAs("tmp/leadetawidth.pdf");
  c1 = new TCanvas();
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiWidthGraph[i]->SetLineColor(i+1);
    subPhiWidthGraph[i]->SetMarkerStyle(29);
    subPhiWidthGraph[i]->SetMarkerSize(3);
    subPhiWidthGraph[i]->SetMarkerColor(i+1);
    subPhiWidthGraph[i]->SetTitle("Recoil Jet #Delta#phi - Gaussian Widths");
    subPhiWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    subPhiWidthGraph[i]->GetYaxis()->SetTitle("Width");
    subPhiWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      subPhiWidthGraph[i]->Draw("P");
    else
      subPhiWidthGraph[i]->Draw("P");
  }
  c1->SaveAs("tmp/subphiwidth.pdf");
  c1 = new TCanvas();
  for ( int i = 0; i < nFiles; ++i ) {
    subPhiDifWidthGraph[i]->SetLineColor(i+1);
    subPhiDifWidthGraph[i]->SetMarkerStyle(29);
    subPhiDifWidthGraph[i]->SetMarkerSize(3);
    subPhiDifWidthGraph[i]->SetMarkerColor(i+1);
    subPhiDifWidthGraph[i]->SetTitle("Recoil Jet #Delta#eta Subtracted #Delta#phi - Gaussian Widths");
    subPhiDifWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    subPhiDifWidthGraph[i]->GetYaxis()->SetTitle("Width");
    subPhiDifWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      subPhiDifWidthGraph[i]->Draw("P");
    else
      subPhiDifWidthGraph[i]->Draw("P");
  }
  c1->SaveAs("tmp/subphidifwidth.pdf");
  c1 = new TCanvas();
  for ( int i = 0; i < nFiles; ++i ) {
    subEtaWidthGraph[i]->SetLineColor(i+1);
    subEtaWidthGraph[i]->SetMarkerStyle(29);
    subEtaWidthGraph[i]->SetMarkerSize(3);
    subEtaWidthGraph[i]->SetMarkerColor(i+1);
    subEtaWidthGraph[i]->SetTitle("Recoil Jet #Delta#eta - Gaussian Widths");
    subEtaWidthGraph[i]->GetXaxis()->SetTitle("p_{T}");
    subEtaWidthGraph[i]->GetYaxis()->SetTitle("Width");
    subEtaWidthGraph[i]->GetYaxis()->SetRangeUser( 0, 1);
    if ( i == 0)
      subEtaWidthGraph[i]->Draw("P");
    else
      subEtaWidthGraph[i]->Draw("P");
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
    
    leadPhiGraph[i]->Draw("P");
    leadPhiDifGraph[i]->Draw("P");
    leadEtaGraph[i]->Draw("P");
    
    TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);
    leg->SetHeader("Trigger Jet Yields"); // option "C" allows to center the header
    leg->AddEntry(leadPhiGraph[i],"#Delta#phi","lep");
    leg->AddEntry(leadPhiDifGraph[i],"#Delta#phi #Delta#eta subtracted","lep");
    leg->AddEntry(leadEtaGraph[i],"#Delta#eta","lep");
    leg->Draw();
    
    TString graphOutName = "tmp/graph_out_"+analysisNames[i]+"_lead.pdf";
    
    c1->SaveAs( graphOutName );
    
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
    
    subPhiGraph[i]->Draw("P");
    subPhiDifGraph[i]->Draw("P");
    subEtaGraph[i]->Draw("P");
    
    leg = new TLegend(0.65,0.7,0.9,0.9);
    leg->SetHeader("Recoil Jet Yields"); // option "C" allows to center the header
    leg->AddEntry(subPhiGraph[i],"#Delta#phi","lep");
    leg->AddEntry(subPhiDifGraph[i],"#Delta#phi #Delta#eta subtracted","lep");
    leg->AddEntry(subEtaGraph[i],"#Delta#eta","lep");
    leg->Draw();

    graphOutName = "tmp/graph_out_"+analysisNames[i]+"_sub.pdf";

    c1->SaveAs( graphOutName );
    
    // joern asked for dphi without subtraction overlayed
    graphOutName = "tmp/dphi_trig_recoil_yield_"+analysisNames[i]+".pdf";
    c1 = new TCanvas();
    leadPhiGraph[i]->SetTitle("");
    leadPhiGraph[i]->SetLineColor(1);
    leadPhiGraph[i]->SetMarkerColor(1);
    subPhiGraph[i]->SetLineColor(7);
    subPhiGraph[i]->SetMarkerColor(7);
    subPhiGraph[i]->SetMarkerStyle(20);
    
    leadPhiGraph[i]->Draw("P");
    subPhiGraph[i]->Draw("P");
    
    leg = new TLegend(0.65,0.7,0.9,0.9);
    leg->SetHeader("#Delta#phi Jet Yields");
    leg->AddEntry(leadPhiGraph[i],"Trigger","lep");
    leg->AddEntry(subPhiGraph[i],"Recoil","lep");
    leg->Draw();
    
    c1->SaveAs( graphOutName );
  }
  
  
  return 0;
}
