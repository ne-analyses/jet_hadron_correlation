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
  TString ptBinString[nPtBins] = { "0.5-1.0", "1.0-2.0", "2.0-3.0", "3.0-4.0", "4.0-6.0" };
  
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
  corrCentVz.resize( nFiles );
  mixCentVz.resize( nFiles );
  for ( int i = 0; i < nFiles; ++i ) {
    corrCentVz[i].resize( corrAnalysis::binsCentrality );
    mixCentVz[i].resize( corrAnalysis::binsCentrality );
    
    for ( int j = 0; j < corrAnalysis::binsCentrality; ++j ) {
      corrCentVz[i][j].resize( corrAnalysis::binsVz );
      mixCentVz[i][j].resize( corrAnalysis::binsVz );
    }
  }
  
  for ( int i = 0; i < nFiles; ++i ) {
    
    nEvents[i] = corrFiles[i]->Get( "nevents" );
    hVz[i] = corrFiles[i]->Get( "vzdist" );
    corrHist[i] = corrFiles[i]->Get( "leadjetcorr" );
    mixHist[i] = mixFiles[i]->Get( "leadjetcorr" );
  }
  
  return 0;
}









//int jet_dijet() {
//  
//  // Define Pi to float precision
//  const float pi = 3.141592;
//  gStyle->SetOptStat(0);
//  
//  TFile* dijetFile = new TFile( "histograms/auaudijet.root", "READ");
//  TFile* jetLoFile = new TFile( "histograms/auaujetlo.root", "READ");
//  TFile* jetHiFile = new TFile( "histograms/auaujethi.root", "READ");
//  TFile* ppFile = new TFile( "histograms/ppdijet.root", "READ");
//  TFile* ppJetLoFile = new TFile( "histograms/ppjetlo.root", "READ");
//  TFile* ppJetHiFile = new TFile( "histograms/ppjethi.root", "READ");
//  
//  
//  TH3D* dijet_corr = (TH3D*) dijetFile->Get("leadjetcorr");
//  TH3D* dijetEvents = (TH3D*) dijetFile->Get("nevents");
//  dijet_corr->SetName("dijet_3d_corr");
//  dijetEvents->SetName("dijetEvents");
//  TH2D* auau_dijet = (TH2D*) dijet_corr->Project3D("ZY");
//  
//  TH3D* jetLo_corr = (TH3D*) jetLoFile->Get("triggerjetcorr");
//  TH3D* jetLoEvents = (TH3D*) jetLoFile->Get("nevents");
//  jetLo_corr->SetName("jetLo_3d_corr");
//  jetLoEvents->SetName("jetLoEvents");
//  TH2D* auau_jetlo = (TH2D*) jetLo_corr->Project3D("ZY");
//  
//  TH3D* jetHi_corr = (TH3D*) jetHiFile->Get("triggerjetcorr");
//  TH3D* jetHiEvents = (TH3D*) jetHiFile->Get("nevents");
//  jetHi_corr->SetName("jetHi_3d_corr");
//  jetHiEvents->SetName("jetHiEvents");
//  TH2D* auau_jethi = (TH2D*) jetHi_corr->Project3D("ZY");
//  
//  TH3D* pp_corr = (TH3D*) ppFile->Get("ppleadjetcorr");
//  TH3D* ppEvents = (TH3D*) ppFile->Get("binvzdist");
//  pp_corr->SetName("pp_3d_corr");
//  ppEvents->SetName("ppEvents");
//  TH2D* pp_dijet = (TH2D*) pp_corr->Project3D("ZY");
//  
//  TH3D* pp_lo_corr = (TH3D*) ppJetLoFile->Get("pptriggerjetcorr");
//  TH3D* ppLoEvents = (TH3D*) ppJetLoFile->Get("binvzdist");
//  pp_lo_corr->SetName("pp_lo_3d_corr");
//  ppLoEvents->SetName("ppLoEvents");
//  TH2D* pp_jetlo = (TH2D*) pp_lo_corr->Project3D("ZY");
//  
//  TH3D* pp_hi_corr = (TH3D*) ppJetHiFile->Get("pptriggerjetcorr");
//  TH3D* ppHiEvents = (TH3D*) ppJetHiFile->Get("binvzdist");
//  pp_hi_corr->SetName("pp_hi_3d_corr");
//  ppHiEvents->SetName("ppHiEvents");
//  TH2D* pp_jethi = (TH2D*) pp_hi_corr->Project3D("ZY");
//  
//  cout<<auau_dijet<<std::endl;
//  cout<< auau_jethi<<std::endl;
//  cout<< auau_jetlo<<std::endl;
//  cout<< pp_dijet<<std::endl;
//  cout<< pp_jetlo <<std::endl;
//  cout<< pp_jethi <<std::endl;
//  
//  // Pt bins
//  const int nPtBins = 5;
//  double ptBinLo[nPtBins] = { 3, 5, 9, 13, 17 };
//  double ptBinHi[nPtBins] = { 4, 8, 12, 16, 24 };
//  TString ptBinString[nPtBins] = { "0.5-1.0", "1.0-2.0", "2.0-3.0", "3.0-4.0", "4.0-6.0" };
//  
//  // Project
//  TH1D* dijet[nPtBins];
//  TH1D* jetlo[nPtBins];
//  TH1D* jethi[nPtBins];
//  TH1D* ppdijet[nPtBins];
//  TH1D* ppjetlo[nPtBins];
//  TH1D* ppjethi[nPtBins];
//  TF1* dijet_fit[nPtBins];
//  TF1* jetlo_fit[nPtBins];
//  TF1* jethi_fit[nPtBins];
//  TF1* ppdijet_fit[nPtBins];
//  TF1* ppjetlo_fit[nPtBins];
//  TF1* ppjethi_fit[nPtBins];
//  
//  // Fit definitions
//  TString phi_form = "[0]+[1]*exp(-0.5*((x-[2])/[3])**2)+[4]*exp(-0.5*((x-[5])/[6])**2)";
//  double phiMin = -pi + pi/2.0;
//  double phiMax = pi + pi/2.0;
//  
//  // Now for the yields
//  double auaudijet_yield[nPtBins];
//  double auaujetlo_yield[nPtBins];
//  double auaujethi_yield[nPtBins];
//  double ppdijet_yield[nPtBins];
//  double ppjetlo_yield[nPtBins];
//  double ppjethi_yield[nPtBins];
//  
//  double auaudijet_yield_away[nPtBins];
//  double auaujetlo_yield_away[nPtBins];
//  double auaujethi_yield_away[nPtBins];
//  double ppdijet_yield_away[nPtBins];
//  double ppjetlo_yield_away[nPtBins];
//  double ppjethi_yield_away[nPtBins];
//  
//  for (int i = 0; i < nPtBins; ++i ) {
//    
//    std::cout << std::endl;
//    std::cout << "NEW PT BIN: " << i << std::endl;
//    std::cout << std::endl;
//    
//    TString dijet_name = "dijet_1d_";
//    TString jetlo_name = "jetlo_1d_";
//    TString jethi_name = "jethi_id_";
//    TString ppdijet_name = "ppdijet_1d_";
//    TString ppjetlo_name = "ppjetlo_1d_";
//    TString ppjethi_name = "ppjethi_1d_";
//    
//    dijet_name += i;
//    jetlo_name += i;
//    jethi_name += i;
//    ppdijet_name += i;
//    ppjetlo_name += i;
//    ppjethi_name += i;
//    
//    std::cout<< "BIN WIDTHS" << std::endl;
//    std::cout<< "Low edge: " << auau_dijet->GetYaxis()->GetBinLowEdge( ptBinLo[i] ) << std::endl;
//    std::cout<< "High edge: " << auau_dijet->GetYaxis()->GetBinUpEdge( ptBinHi[i] ) << std::endl;
//    
//    dijet[i] = (TH1D*) auau_dijet->ProjectionX( dijet_name, ptBinLo[i], ptBinHi[i] );
//    jetlo[i] = (TH1D*) auau_jetlo->ProjectionX( jetlo_name, ptBinLo[i], ptBinHi[i] );
//    jethi[i] = (TH1D*) auau_jethi->ProjectionX( jethi_name, ptBinLo[i], ptBinHi[i] );
//    ppdijet[i] = (TH1D*) pp_dijet->ProjectionX( ppdijet_name, ptBinLo[i], ptBinHi[i] );
//    ppjetlo[i] = (TH1D*) pp_jetlo->ProjectionX( ppjetlo_name, ptBinLo[i], ptBinHi[i] );
//    ppjethi[i] = (TH1D*) pp_jethi->ProjectionX( ppjethi_name, ptBinLo[i], ptBinHi[i] );
//    
//    // Scale
//    dijet[i]->Scale( 1/double(dijetEvents->GetEntries()) );
//    dijet[i]->Scale( 1/dijet[i]->GetBinWidth(1) );
//    
//    jetlo[i]->Scale( 1/double(jetLoEvents->GetEntries()) );
//    jetlo[i]->Scale( 1/jetlo[i]->GetBinWidth(1) );
//    
//    jethi[i]->Scale( 1/double(jetHiEvents->GetEntries()) );
//    jethi[i]->Scale( 1/jethi[i]->GetBinWidth(1) );
//    
//    ppdijet[i]->Scale( 1/double(ppEvents->GetEntries()) );
//    ppdijet[i]->Scale( 1/ppdijet[i]->GetBinWidth(1) );
//    
//    ppjetlo[i]->Scale( 1/double(ppLoEvents->GetEntries()) );
//    ppjetlo[i]->Scale( 1/ppjetlo[i]->GetBinWidth(1) );
//    
//    ppjethi[i]->Scale( 1/double(ppHiEvents->GetEntries()) );
//    ppjethi[i]->Scale( 1/ppjethi[i]->GetBinWidth(1) );
//    
//    
//    if ( i == 0 ) {
//      dijet[i]->Scale( 1.0/0.65 );
//      jetlo[i]->Scale( 1.0/0.65 );
//      jethi[i]->Scale( 1.0/0.65 );
//      ppdijet[i]->Scale( 1.0/0.8 );
//      ppjetlo[i]->Scale( 1.0/0.8 );
//      ppjethi[i]->Scale( 1.0/0.8 );
//    }
//    if ( i == 1 ) {
//      dijet[i]->Scale( 1.0/0.67 );
//      jetlo[i]->Scale( 1.0/0.67 );
//      jethi[i]->Scale( 1.0/0.67 );
//      ppdijet[i]->Scale( 1.0/0.8 );
//      ppjetlo[i]->Scale( 1.0/0.8 );
//      ppjethi[i]->Scale( 1.0/0.8 );
//    }
//    if ( i == 2 ) {
//      dijet[i]->Scale( 1.0/0.68 );
//      jetlo[i]->Scale( 1.0/0.68 );
//      jethi[i]->Scale( 1.0/0.68 );
//      ppdijet[i]->Scale( 1.0/0.8 );
//      ppjetlo[i]->Scale( 1.0/0.8 );
//      ppjethi[i]->Scale( 1.0/0.8 );
//    }
//    if ( i == 3 ) {
//      dijet[i]->Scale( 1.0/0.7 );
//      jetlo[i]->Scale( 1.0/0.7 );
//      jethi[i]->Scale( 1.0/0.7 );
//      ppdijet[i]->Scale( 1.0/0.8 );
//      ppjetlo[i]->Scale( 1.0/0.8 );
//      ppjethi[i]->Scale( 1.0/0.8 );
//    }
//    if ( i == 4 ) {
//      dijet[i]->Scale( 1.0/0.71 );
//      jetlo[i]->Scale( 1.0/0.71 );
//      jethi[i]->Scale( 1.0/0.71 );
//      ppdijet[i]->Scale( 1.0/0.8 );
//      ppjetlo[i]->Scale( 1.0/0.8 );
//      ppjethi[i]->Scale( 1.0/0.8 );
//    }
//    
//    
//    
//    // Now subtract the background!
//    TString dijet_name_fit = "dijet_fit_"; dijet_name_fit += i;
//    TString jetlo_name_fit = "jetlo_fit_"; jetlo_name_fit += i;
//    TString jethi_name_fit = "jethi_fit_"; jethi_name_fit += i;
//    TString ppdijet_name_fit = "ppdijet_fit_"; ppdijet_name_fit += i;
//    TString ppjetlo_name_fit = "ppjetlo_fit_"; ppjetlo_name_fit += i;
//    TString ppjethi_name_fit = "ppjethi_fit_"; ppjethi_name_fit += i;
//    
//    // define each fit, set the basics
//    TF1* dijet_tmp = new TF1( dijet_name_fit, phi_form, phiMin, phiMax );
//    dijet_tmp->FixParameter( 2, 0 );
//    dijet_tmp->FixParameter( 5, pi );
//    dijet_tmp->SetParameter( 3, 0.2 );
//    dijet_tmp->SetParameter( 6, 0.2 );
//    TF1* jetlo_tmp = new TF1( jetlo_name_fit, phi_form, phiMin, phiMax );
//    jetlo_tmp->FixParameter( 2, 0 );
//    jetlo_tmp->FixParameter( 5, pi );
//    jetlo_tmp->SetParameter( 3, 0.2 );
//    jetlo_tmp->SetParameter( 6, 0.2 );
//    TF1* jethi_tmp = new TF1( jethi_name_fit, phi_form, phiMin, phiMax );
//    jethi_tmp->FixParameter( 2, 0 );
//    jethi_tmp->FixParameter( 5, pi );
//    jethi_tmp->SetParameter( 3, 0.2 );
//    jethi_tmp->SetParameter( 6, 0.2 );
//    TF1* ppdijet_tmp = new TF1( ppdijet_name_fit, phi_form, phiMin, phiMax );
//    ppdijet_tmp->FixParameter( 2, 0 );
//    ppdijet_tmp->FixParameter( 5, pi );
//    ppdijet_tmp->SetParameter( 3, 0.2 );
//    ppdijet_tmp->SetParameter( 6, 0.2 );
//    TF1* ppjetlo_tmp = new TF1( ppjetlo_name_fit, phi_form, phiMin, phiMax );
//    ppjetlo_tmp->FixParameter( 2, 0 );
//    ppjetlo_tmp->FixParameter( 5, pi );
//    ppjetlo_tmp->SetParameter( 3, 0.2 );
//    ppjetlo_tmp->SetParameter( 6, 0.2 );
//    TF1* ppjethi_tmp = new TF1( ppjethi_name_fit, phi_form, phiMin, phiMax );
//    ppjethi_tmp->FixParameter( 2, 0 );
//    ppjethi_tmp->FixParameter( 5, pi );
//    ppjethi_tmp->SetParameter( 3, 0.2 );
//    ppjethi_tmp->SetParameter( 6, 0.2 );
//    
//    
//    dijet[i]->Fit( dijet_tmp );
//    jetlo[i]->Fit( jetlo_tmp );
//    jethi[i]->Fit( jethi_tmp );
//    ppdijet[i]->Fit( ppdijet_tmp );
//    ppjetlo[i]->Fit( ppjetlo_tmp );
//    ppjethi[i]->Fit( ppjethi_tmp );
//    
//    
//    // Now subtract the background
//    dijet_name_fit += "_sub";
//    jetlo_name_fit += "_sub";
//    jethi_name_fit += "_sub";
//    ppdijet_name_fit += "_sub";
//    ppjetlo_name_fit += "_sub";
//    ppjethi_name_fit += "_sub";
//    
//    TF1* dijet_sub = new TF1( dijet_name_fit, "[0]", phiMin, phiMax );
//    dijet_sub->SetParameter( 0, dijet_tmp->GetParameter(0) );
//    TF1* jetlo_sub = new TF1( jetlo_name_fit, "[0]", phiMin, phiMax );
//    jetlo_sub->SetParameter( 0, jetlo_tmp->GetParameter(0) );
//    TF1* jethi_sub = new TF1( jethi_name_fit, "[0]", phiMin, phiMax );
//    jethi_sub->SetParameter( 0, jethi_tmp->GetParameter(0) );
//    
//    TF1* ppdijet_sub = new TF1( ppdijet_name_fit, "[0]", phiMin, phiMax );
//    ppdijet_sub->SetParameter( 0, ppdijet_tmp->GetParameter(0) );
//    TF1* ppjetlo_sub = new TF1( ppjetlo_name_fit, "[0]", phiMin, phiMax );
//    ppjetlo_sub->SetParameter( 0, ppjetlo_tmp->GetParameter(0) );
//    TF1* ppjethi_sub = new TF1( ppjethi_name_fit, "[0]", phiMin, phiMax );
//    ppjethi_sub->SetParameter( 0, ppjethi_tmp->GetParameter(0) );
//    
//    //if ( i == 1 )
//    //jethi_sub->SetParameter( 0, jethi_tmp->GetParameter(0) - 0.3 );
//    
//    dijet[i]->Add( dijet_sub, -1 );
//    jetlo[i]->Add( jetlo_sub, -1 );
//    jethi[i]->Add( jethi_sub, -1 );
//    ppdijet[i]->Add( ppdijet_sub, -1 );
//    ppjetlo[i]->Add( ppjetlo_sub, -1 );
//    ppjethi[i]->Add( ppjethi_sub, -1 );
//    
//    // Actual Fit
//    dijet_name += "_fit";
//    jetlo_name += "_fit";
//    jethi_name += "_fit";
//    ppdijet_name += "_fit";
//    ppjetlo_name += "_fit";
//    ppjethi_name += "_fit";
//    
//    dijet_fit[i] = new TF1( dijet_name, phi_form, phiMin, phiMax );
//    dijet_fit[i]->SetLineColor(kBlack);
//    dijet_fit[i]->FixParameter( 2, 0 );
//    dijet_fit[i]->FixParameter( 5, pi );
//    dijet_fit[i]->SetParameter( 3, 0.2 );
//    dijet_fit[i]->SetParameter( 6, 0.2 );
//    jetlo_fit[i] = new TF1( jetlo_name, phi_form, phiMin, phiMax );
//    jetlo_fit[i]->SetLineColor(kRed);
//    jetlo_fit[i]->FixParameter( 2, 0 );
//    jetlo_fit[i]->FixParameter( 5, pi );
//    jetlo_fit[i]->SetParameter( 3, 0.2 );
//    jetlo_fit[i]->SetParameter( 6, 0.2 );
//    jethi_fit[i] = new TF1( jethi_name, phi_form, phiMin, phiMax );
//    jethi_fit[i]->SetLineColor(kBlue);
//    jethi_fit[i]->FixParameter( 2, 0 );
//    jethi_fit[i]->FixParameter( 5, pi );
//    jethi_fit[i]->SetParameter( 3, 0.2 );
//    jethi_fit[i]->SetParameter( 6, 0.2 );
//    ppdijet_fit[i] = new TF1( ppdijet_name, phi_form, phiMin, phiMax );
//    ppdijet_fit[i]->FixParameter( 2, 0 );
//    ppdijet_fit[i]->FixParameter( 5, pi );
//    ppdijet_fit[i]->SetParameter( 3, 0.2 );
//    ppdijet_fit[i]->SetParameter( 6, 0.2 );
//    ppjetlo_fit[i] = new TF1( ppjetlo_name, phi_form, phiMin, phiMax );
//    ppjetlo_fit[i]->FixParameter( 2, 0 );
//    ppjetlo_fit[i]->FixParameter( 5, pi );
//    ppjetlo_fit[i]->SetParameter( 3, 0.2 );
//    ppjetlo_fit[i]->SetParameter( 6, 0.2 );
//    ppjethi_fit[i] = new TF1( ppjethi_name, phi_form, phiMin, phiMax );
//    ppjethi_fit[i]->FixParameter( 2, 0 );
//    ppjethi_fit[i]->FixParameter( 5, pi );
//    ppjethi_fit[i]->SetParameter( 3, 0.2 );
//    ppjethi_fit[i]->SetParameter( 6, 0.2 );
//    
//    std::cout << std::endl;
//    std::cout << "------------------------------------------------------" << std::endl;
//    std::cout << " FINAL FITS FOR BIN: "<< i << std::endl;
//    std::cout << std::endl;
//    
//    
//    std::cout<<"DIJET"<<std::endl;
//    dijet[i]->SetLineColor(kBlack);
//    dijet_fit[i]->SetLineColor(kBlack);
//    dijet[i]->Fit( dijet_fit[i] );
//    std::cout<<std::endl;
//    std::cout<<"LOW JET" << std::endl;
//    jetlo[i]->SetLineColor(kRed);
//    jetlo_fit[i]->SetLineColor(kRed);
//    jetlo[i]->Fit( jetlo_fit[i] );
//    std::cout<<std::endl;
//    std::cout<<" HIGH JET" << std::endl;
//    jethi[i]->SetLineColor(kBlue);
//    jethi_fit[i]->SetLineColor(kBlue);
//    jethi[i]->Fit( jethi_fit[i] );
//    std::cout<<std::endl;
//    std::cout<<" PP DIJET "<< std::endl;
//    ppdijet[i]->SetLineColor(8);
//    ppdijet_fit[i]->SetLineColor(8);
//    ppdijet[i]->Fit( ppdijet_fit[i] );
//    std::cout<<" PP LOW JET"<< std::endl;
//    ppjetlo[i]->SetLineColor(46);
//    ppjetlo_fit[i]->SetLineColor(46);
//    ppjetlo[i]->Fit( ppjetlo_fit[i] );
//    std::cout<<" PP HIGH JET"<< std::endl;
//    ppjethi[i]->SetLineColor(36);
//    ppjethi_fit[i]->SetLineColor(36);
//    ppjethi[i]->Fit( ppjethi_fit[i] );
//    
//    std::cout << std::endl;
//    std::cout << " FINAL FITS FOR BIN: "<< i << " COMPLETE" << std::endl;
//    std::cout << "------------------------------------------------------" << std::endl;
//    std::cout << std::endl;
//    
//    // now pull out the integrals...
//    auaudijet_yield[i] =  sqrt(2.0*pi)*dijet_fit[i]->GetParameter(1)*dijet_fit[i]->GetParameter(3);
//    if ( auaudijet_yield[i] < 0.0 )
//      auaudijet_yield[i] = -1.0*auaudijet_yield[i];
//    
//    auaujetlo_yield[i] =  sqrt(2.0*pi)*jetlo_fit[i]->GetParameter(1)*jetlo_fit[i]->GetParameter(3);
//    if ( auaujetlo_yield[i] < 0.0 )
//      auaujetlo_yield[i] = -1.0*auaujetlo_yield[i];
//    
//    auaujethi_yield[i] =  sqrt(2.0*pi)*jethi_fit[i]->GetParameter(1)*jethi_fit[i]->GetParameter(3);
//    if ( auaujethi_yield[i] < 0.0 )
//      auaujethi_yield[i] = -1.0*auaujethi_yield[i];
//    
//    ppdijet_yield[i] = sqrt(2.0*pi)*ppdijet_fit[i]->GetParameter(1)*ppdijet_fit[i]->GetParameter(3);
//    if ( ppdijet_yield[i] < 0.0 )
//      ppdijet_yield[i] = -1.0*ppdijet_yield[i];
//    ppjetlo_yield[i] = sqrt(2.0*pi)*ppjetlo_fit[i]->GetParameter(1)*ppjetlo_fit[i]->GetParameter(3);
//    if ( ppjetlo_yield[i] < 0.0 )
//      ppjetlo_yield[i] = -1.0*ppjetlo_yield[i];
//    ppjethi_yield[i] = sqrt(2.0*pi)*ppjethi_fit[i]->GetParameter(1)*ppjethi_fit[i]->GetParameter(3);
//    if ( ppjethi_yield[i] < 0.0 )
//      ppjethi_yield[i] = -1.0*ppjethi_yield[i];
//    
//    auaudijet_yield_away[i] = sqrt(2.0*pi)*dijet_fit[i]->GetParameter(4)*dijet_fit[i]->GetParameter(6);
//    auaujetlo_yield_away[i] = sqrt(2.0*pi)*jetlo_fit[i]->GetParameter(4)*jetlo_fit[i]->GetParameter(6);
//    auaujethi_yield_away[i] = sqrt(2.0*pi)*jethi_fit[i]->GetParameter(4)*jethi_fit[i]->GetParameter(6);
//    ppdijet_yield_away[i] = sqrt(2.0*pi)*ppdijet_fit[i]->GetParameter(4)*ppdijet_fit[i]->GetParameter(6);
//    ppjetlo_yield_away[i] = sqrt(2.0*pi)*ppjetlo_fit[i]->GetParameter(4)*ppjetlo_fit[i]->GetParameter(6);
//    ppjethi_yield_away[i] = sqrt(2.0*pi)*ppjethi_fit[i]->GetParameter(4)*ppjethi_fit[i]->GetParameter(6);
//    
//  }
//  
//  
//  TCanvas* canvas;
//  TLegend* leg;
//  for ( int i = 0; i < nPtBins; ++i ) {
//    canvas = new TCanvas();
//    leg = new TLegend(0.1,0.7,0.48,0.9);
//    
//    dijet[i]->SetLineColor(kBlack);
//    dijet[i]->GetXaxis()->SetTitle("#Delta #phi");
//    dijet[i]->GetYaxis()->SetTitle("1/N_{dijets} dN/d#Delta#phi");
//    if ( i == 0 ) {
//      dijet[i]->SetTitle(" 0.5 < pt_{assoc} < 1.0");
//      dijet[i]->GetYaxis()->SetRangeUser( -0.5, 4);
//    }
//    else if ( i == 1 ) {
//      dijet[i]->SetTitle(" 1.0 < pt_{assoc} < 2.0 ");
//    }
//    else if ( i == 2 ) {
//      dijet[i]->SetTitle(" 2.0 < pt_{assoc} < 3.0 ");
//    }
//    else if ( i == 3 ) {
//      dijet[i]->SetTitle(" 3.0 < pt_{assoc} < 4.0 ");
//      dijet[i]->GetYaxis()->SetRangeUser( -0.2, 2.0 );
//    }
//    else if ( i == 4 ) {
//      dijet[i]->SetTitle(" 4.0 < pt_{assoc} < 6.0 ");
//      dijet[i]->GetYaxis()->SetRangeUser( -0.2, 3.0);
//    }
//    
//    dijet_fit[i]->SetLineColor(kBlack);
//    jetlo[i]->SetLineColor(kRed);
//    jetlo_fit[i]->SetLineColor(kRed);
//    jethi[i]->SetLineColor(kBlue);
//    jethi_fit[i]->SetLineColor(kBlue);
//    
//    //leg->SetHeader("Correlation comparisons","C"); // option "C" allows to center the header
//    leg->AddEntry(dijet_fit[i],"dijet","lep");
//    leg->AddEntry(jetlo_fit[i],"jet 10 < pt < 15", "lep");
//    leg->AddEntry(jethi_fit[i], "jet 20 < pt< 40", "lep");
//    
//    //dijet[i]->SetTitle("");
//    dijet[i]->Draw();
//    jetlo[i]->Draw("SAME");
//    jethi[i]->Draw("SAME");
//    //leg->Draw();
//    
//    TString outname = "tmp/auau_pt_";
//    outname += i;
//    outname += ".pdf";
//    canvas->SaveAs(outname);
//  }
//  
//  for ( int i = 0; i < nPtBins; ++i ) {
//    canvas = new TCanvas();
//    leg = new TLegend(0.1,0.7,0.48,0.9);
//    
//    ppdijet[i]->SetLineColor(kBlack);
//    ppdijet[i]->GetXaxis()->SetTitle("#Delta #phi");
//    ppdijet[i]->GetYaxis()->SetTitle("1/N_{dijets} dN/d#Delta#phi");
//    if ( i == 0 ) {
//      ppdijet[i]->SetTitle(" 0.5 < pt_{assoc} < 1.0");
//      ppdijet[i]->GetYaxis()->SetRangeUser( -0.5, 4);
//    }
//    else if ( i == 1 ) {
//      ppdijet[i]->SetTitle(" 1.0 < pt_{assoc} < 2.0 ");
//    }
//    else if ( i == 2 ) {
//      ppdijet[i]->SetTitle(" 2.0 < pt_{assoc} < 3.0 ");
//    }
//    else if ( i == 3 ) {
//      ppdijet[i]->SetTitle(" 3.0 < pt_{assoc} < 4.0 ");
//      ppdijet[i]->GetYaxis()->SetRangeUser( -0.2, 2.0 );
//    }
//    else if ( i == 4 ) {
//      ppdijet[i]->SetTitle(" 4.0 < pt_{assoc} < 6.0 ");
//      ppdijet[i]->GetYaxis()->SetRangeUser( -0.2, 3.0);
//    }
//    
//    ppdijet_fit[i]->SetLineColor(kBlack);
//    ppjetlo[i]->SetLineColor(kRed);
//    ppjetlo_fit[i]->SetLineColor(kRed);
//    ppjethi[i]->SetLineColor(kBlue);
//    ppjethi_fit[i]->SetLineColor(kBlue);
//    
//    //leg->SetHeader("Correlation comparisons","C"); // option "C" allows to center the header
//    //leg->AddEntry(dijet_fit[i],"dijet","lep");
//    //leg->AddEntry(jetlo_fit[i],"jet 10 < pt < 15", "lep");
//    //leg->AddEntry(jethi_fit[i], "jet 20 < pt< 40", "lep");
//    
//    //ppdijet[i]->SetTitle("");
//    ppdijet[i]->Draw();
//    ppjetlo[i]->Draw("SAME");
//    ppjethi[i]->Draw("SAME");
//    //leg->Draw();
//    
//    TString outname = "tmp/pp_pt_";
//    outname += i;
//    outname += ".pdf";
//    canvas->SaveAs(outname);
//  }
//  
//  
//  for ( int i = 0; i < nPtBins; ++i ) {
//    
//    dijet[i]->Draw();
//    ppdijet[i]->Draw("SAME");
//    TString tmpName = "tmp/dijetpp_";
//    tmpName += i; tmpName += ".pdf";
//    canvas->SaveAs( tmpName );
//  }
//  
//  
//  // Now we want to subtract the hi jet from the dijet
//  for ( int i = 0; i < nPtBins; ++i ) {
//    TH1D* tmpAuAuDijet = dijet[i]->Clone("auau_dijet_dif");
//    TH1D* tmpPPDijet = ppdijet[i]->Clone("pp_dijet_dif");
//    TH1D* tmpAuAuHiJet = jethi[i]->Clone("auau_jethi_dif");
//    TH1D* tmpPPHiJet = ppjethi[i]->Clone("pp_jethi_dif");
//    
//    tmpAuAuDijet->Add(tmpAuAuHiJet, -1 );
//    tmpPPDijet->Add( tmpPPHiJet, -1 );
//    
//    TString tmpName = "tmp/diff_";
//    tmpName += i;
//    tmpName += ".pdf";
//    TString tmpFitName = "dif_fit_";
//    TString ppTmpFitName = "pp_dif_fit_";
//    tmpFitName += i;
//    ppTmpFitName += i;
//    
//    TF1* auau_dif_tmp = new TF1( tmpFitName, phi_form, phiMin, phiMax );
//    auau_dif_tmp->SetLineColor(kBlack);
//    auau_dif_tmp->FixParameter( 2, 0 );
//    auau_dif_tmp->FixParameter( 5, pi );
//    auau_dif_tmp->SetParameter( 3, 0.2 );
//    auau_dif_tmp->SetParameter( 6, 0.2 );
//    TF1* pp_dif_tmp = new TF1( ppTmpFitName, phi_form, phiMin, phiMax );
//    tmpPPDijet->SetLineColor(kBlue);
//    pp_dif_tmp->SetLineColor(kBlue);
//    pp_dif_tmp->FixParameter( 2, 0 );
//    pp_dif_tmp->FixParameter( 5, pi );
//    pp_dif_tmp->SetParameter( 3, 0.2 );
//    pp_dif_tmp->SetParameter( 6, 0.2 );
//    
//    tmpAuAuDijet->Fit( auau_dif_tmp );
//    tmpPPDijet->Fit(pp_dif_tmp );
//    
//    //tmpAuAuDijet->SetTitle("");
//    tmpAuAuDijet->Draw();
//    tmpPPDijet->Draw("SAME");
//    canvas->SaveAs(tmpName);
//    
//  }
//  
//  for (int i = 0; i < nPtBins; ++i ) {
//    cout<<"bin: "<< i << " dijet: "<< auaudijet_yield[i] << std::endl;
//    cout<<"bin: "<< i << " dijet: "<< auaudijet_yield_away[i] << std::endl;
//    cout<<"bin: "<< i << " jetlo: "<< auaujetlo_yield[i] << std::endl;
//    cout<<"bin: "<< i << " jetlo: "<< auaujetlo_yield_away[i] << std::endl;
//    cout<<"bin: "<< i << " jethi: "<< auaujethi_yield[i] << std::endl;
//    cout<<"bin: "<< i << " jethi: "<< auaujethi_yield_away[i] << std::endl;
//  }
//  
//  
//  // for the widths
//  double binCenters[nPtBins] = { 0.75, 1.5, 2.5, 3.5, 5 };
//  TGraph* dijetGraph = new TGraph( nPtBins, binCenters, auaudijet_yield );
//  TGraph* dijetGraphAway = new TGraph( nPtBins, binCenters, auaudijet_yield_away );
//  
//  TGraph* jetloGraph = new TGraph( nPtBins, binCenters, auaujetlo_yield );
//  TGraph* jetloGraphAway = new TGraph( nPtBins, binCenters, auaujetlo_yield_away );
//  
//  TGraph* jethiGraph = new TGraph( nPtBins, binCenters, auaujethi_yield );
//  TGraph* jethiGraphAway = new TGraph( nPtBins, binCenters, auaujethi_yield_away );
//  
//  TGraph* ppdijetGraph = new TGraph( nPtBins, binCenters, ppdijet_yield );
//  TGraph* ppdijetGraphAway = new TGraph( nPtBins, binCenters, ppdijet_yield_away );
//  
//  TGraph* ppjetloGraph = new TGraph( nPtBins, binCenters, ppjetlo_yield );
//  TGraph* ppjetloGraphAway = new TGraph( nPtBins, binCenters, ppjetlo_yield_away );
//  
//  TGraph* ppjethiGraph = new TGraph( nPtBins, binCenters, ppjethi_yield );
//  TGraph* ppjethiGraphAway = new TGraph( nPtBins, binCenters, ppjethi_yield_away );
//  
//  TCanvas* AuAuNearCanvas = new TCanvas();
//  
//  dijetGraph->GetYaxis()->SetRangeUser(0, 7);
//  dijetGraph->SetTitle("");
//  dijetGraph->GetXaxis()->SetTitle("p_{T}");
//  dijetGraph->GetYaxis()->SetTitle("Yield");
//  dijetGraph->Draw();
//  jetloGraph->SetLineColor(kRed);
//  jetloGraph->Draw("SAME");
//  jethiGraph->SetLineColor(kBlue);
//  jethiGraph->Draw("SAME");
//  
//  AuAuNearCanvas->SaveAs("tmp/auau_near_yields.pdf");
//  
//  TCanvas* AuAuAwayCanvas = new TCanvas();
//  dijetGraphAway->GetXaxis()->SetTitle("p_{T}");
//  dijetGraphAway->GetYaxis()->SetTitle("Yield");
//  dijetGraphAway->GetYaxis()->SetRangeUser(0, 7);
//  dijetGraphAway->SetTitle("");
//  dijetGraphAway->Draw();
//  jetloGraphAway->SetLineColor(kRed);
//  jetloGraphAway->Draw("SAME");
//  jethiGraphAway->SetLineColor(kBlue);
//  jethiGraphAway->Draw("SAME");
//  
//  AuAuAwayCanvas->SaveAs("tmp/auau_away_yields.pdf");
//  
//  TCanvas* PPNearCanvas = new TCanvas();
//  ppdijetGraph->GetYaxis()->SetRangeUser( 0.0, 2.0 );
//  ppdijetGraph->GetXaxis()->SetTitle("p_{T}");
//  ppdijetGraph->GetYaxis()->SetTitle("Yield");
//  ppdijetGraph->SetTitle("");
//  ppdijetGraph->SetLineColor(8);
//  ppdijetGraph->Draw();
//  ppjetloGraph->SetLineColor(46);
//  ppjetloGraph->Draw("SAME");
//  ppjethiGraph->SetLineColor(36);
//  ppjethiGraph->Draw("SAME");
//  
//  PPNearCanvas->SaveAs("tmp/pp_near_yields.pdf");
//  
//  TCanvas* PPAwayCanvas = new TCanvas();
//  ppdijetGraphAway->GetYaxis()->SetRangeUser( 0.0, 2.0 );
//  ppdijetGraphAway->GetXaxis()->SetTitle("p_{T}");
//  ppdijetGraphAway->GetYaxis()->SetTitle("Yield");
//  ppdijetGraphAway->SetTitle("");
//  ppdijetGraphAway->SetLineColor(8);
//  ppdijetGraphAway->Draw();
//  ppjetloGraphAway->SetLineColor(46);
//  ppjetloGraphAway->Draw("SAME");
//  ppjethiGraphAway->SetLineColor(36);
//  ppjethiGraphAway->Draw("SAME");
//  
//  PPAwayCanvas->SaveAs("tmp/pp_away_yields.pdf");
//  
//  return 0;
//}
