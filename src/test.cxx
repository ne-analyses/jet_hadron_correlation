// testing nick's new correlation implementations

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
  
  // test read in
  std::vector<TFile*> inFile;
  std::vector<TFile*> inFileMix;
  inFile.resize(1);
  inFile[0] = new TFile("out/tmp/auau_corr_20_10.root", "READ");
  inFileMix.resize(1);
  inFileMix[0] = new TFile("out/tmp/auau_mix_20_10.root", "READ");
  
  // testing
  jetHadron::binSelector selector;
  
  std::vector<TH3F*> nEvents;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leading;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > sub;
  std::vector<TH3F*> nEventsMix;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leadingMix;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > subMix;
  
  jetHadron::ReadInFiles( inFile, leading, sub, nEvents, selector );
  jetHadron::ReadInFilesMix( inFileMix, leadingMix, subMix, nEventsMix, selector );
  
  
  std::vector<std::vector<TH2F*> > mixedEvents = jetHadron::RecombineMixedEvents( leadingMix, selector );
  jetHadron::ScaleMixedEvents( mixedEvents );
  
  // Find the pt bin center for future use
  std::vector<TH1F*> ptSpectra;
  std::vector<std::vector<double> > ptBinCenters = jetHadron::FindPtBinCenter( leading, ptSpectra, selector );
  
  
  // Now testing redoing the correlations to 2d in pt
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > leadingCorr;
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > subCorr;
  
  jetHadron::BuildSingleCorrelation( leading, leadingCorr, selector );
  
  std::vector<std::vector<TH2F*> > reducedCorr = jetHadron::AverageCorrelations( leadingCorr, selector );
  
  jetHadron::Print2DHistograms( mixedEvents[0], "tmp/test", "mixing", selector );
  jetHadron::Print2DHistograms( reducedCorr[0], "tmp/corr", "corr", selector);
  
  TH3F* leadingTmp = (TH3F*) inFile[0]->Get("leadjetcorr");
  
  leadingTmp->GetZaxis()->SetRange( selector.ptBinLowEdge(4), selector.ptBinHighEdge(4) );
  
  TCanvas c1;
  ((TH2F*)leadingTmp->Project3D("YX"))->Draw("surf1");
  c1.SaveAs("tmp/leading.pdf");
  
  TH2F* test = 0;
  
  for ( int i = 0; i < leadingCorr.size(); ++i ) {
    for ( int j = 0; j < leadingCorr[i].size(); ++j ) {
      for ( int k = 0; k < leadingCorr[i][j].size(); ++k ) {
        if ( test == 0 ) {
          test = (TH2F*) leadingCorr[i][j][k][4]->Clone();
        }
        else {
          test->Add( leadingCorr[i][j][k][4] );
        }
      }
    }
  }
  
  test->Draw("surf1");
  c1.SaveAs("tmp/compare.pdf");
  
	return 0;
}
