// testing for whatever

#include "corrParameters.hh"
#include "corrFunctions.hh"
#include "outputFunctions.hh"

#include <cmath>
#include <iostream>
#include <random>

// ROOT
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
#include "TCanvas.h"
#include "TStopwatch.h"

// Data is read in by TStarJetPico
// Library, we convert to FastJet::PseudoJet
// TStarJetPico headers
#include "TStarJetPicoReader.h"
#include "TStarJetPicoEvent.h"
#include "TStarJetPicoEventHeader.h"
#include "TStarJetPicoEventCuts.h"
#include "TStarJetPicoPrimaryTrack.h"
#include "TStarJetPicoTower.h"
#include "TStarJetPicoTrackCuts.h"
#include "TStarJetPicoTowerCuts.h"
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "TStarJetPicoTriggerInfo.h"
#include "TStarJetPicoUtils.h"


int main() {
  
  TFile in ( "out/added/pythia/pythia.root", "READ" );
  
  TH3D* corrLead = (TH3D*) in.Get("correlationsLead");
  TH3D* subLead = (TH3D*) in.Get("correlationsSub");
  TH1D* counts = (TH1D*) in.Get("events");
  TH1D* ptCount = (TH1D*) in.Get("ptcount");
  TH1D* ptCountSub = (TH1D*) in.Get("ptcountSub");
  TH1D* ptCountAll = (TH1D*) in.Get("ptcountAll");

  ptCount->Scale( 1.0 / counts->GetEntries() );
  ptCountSub->Scale( 1.0 / counts->GetEntries() );
  ptCountAll->Scale( 1.0 / counts->GetEntries() );
  
  std::vector<double> yieldsLead;
  std::vector<double> yieldsSub;
  
  jetHadron::binSelector selector;
  
  double area = jetHadron::pi*0.4*0.4 / (4*jetHadron::pi);
  
  for ( int i = 0; i < selector.nPtBins; ++i  ) {
    yieldsLead.push_back( (ptCount->Integral( selector.ptBinLowEdge(i), selector.ptBinHighEdge(i) )  - ptCountAll->Integral( selector.ptBinLowEdge(i), selector.ptBinHighEdge(i) ) )/selector.GetPtBinWidth(i)  );
    yieldsSub.push_back(  ( ptCountSub->Integral( selector.ptBinLowEdge(i), selector.ptBinHighEdge(i) ) - ptCountAll->Integral( selector.ptBinLowEdge(i), selector.ptBinHighEdge(i) ))/selector.GetPtBinWidth(i)  );
  }
  
  
  for ( int i = 0; i < yieldsLead.size(); ++i ) {
    std::cout<<"bin: "<< i<< std::endl;
    std::cout<<"lead int: "<< yieldsLead[i]<<std::endl;
    std::cout<<"sub int: "<< yieldsSub[i]<<std::endl;
  }
  
	return 0;
}
