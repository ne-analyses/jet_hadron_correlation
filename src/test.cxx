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
  inFile.resize(1);
  inFile[0] = new TFile("out/tmp/dijet_20_10.root");
  
  // testing
  jetHadron::binSelector selector;
  
  std::vector<TH3F*> nEvents;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leading;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > sub;
  
  jetHadron::ReadInFiles( inFile, leading, sub, nEvents, selector );
  
	return 0;
}
