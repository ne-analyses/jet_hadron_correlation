// testing nick's new correlation implementations

#include "corrParameters.hh"
#include "corrFunctions.hh"

#include <cmath>
#include <iostream>

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
  
  TH1D* test = new TH1D("test", "Effective Acceptance", 500, -2, 2);
  
  // set up generator
  std::random_device rd1;
  std::mt19937 gen1(rd1());
  std::uniform_real_distribution<> random1( -0.6, 0.6 );
  
  // set up generator
  std::random_device rd2;
  std::mt19937 gen2(rd2());
  std::uniform_real_distribution<> random2( -1, 1 );
  
  
  for ( int i = 0; i < 10000; ++i ) {
    double trigger = random1(gen1);
    for ( int j = 0; j < 10000; ++j ) {
      double assoc = random2(gen1);
      test->Fill( trigger - assoc );
    }
  }
  
  TCanvas c1;
  test->Draw();
  c1.SaveAs("test.pdf");
	return 0;
}
