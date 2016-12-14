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
  inFile[0] = new TFile("out/tmp/dijet_corr_20_10.root", "READ");
  inFileMix.resize(1);
  inFile[0] = new TFile("out/tmp/dijet_mix_20_10.root", "READ");
  
  // testing
  jetHadron::binSelector selector;
  
  std::vector<TH3F*> nEvents;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leading;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > sub;
  std::vector<TH3F*> nEventsMix;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > leadingMix;
  std::vector<std::vector<std::vector<std::vector<TH3F*> > > > subMix;
  
  jetHadron::ReadInFiles( inFile, leading, sub, nEvents, selector );
  jetHadron::ReadInFiles( inFileMix, leadingMix, subMix, nEventsMix, selector );
  
  
  std::vector<std::vector<TH2F*> > mixedEvents = jetHadron::RecombineMixedEvents( leadingMix, selector );
  jetHadron::ScaleMixedEvents( mixedEvents );
  
  // Find the pt bin center for future use
  std::vector<TH1F*> ptSpectra;
  std::vector<std::vector<double> > ptBinCenters = jetHadron::FindPtBinCenter( leading, ptSpectra, selector );
  
  jetHadron::Print2DHistograms( mixedEvents[0], "tmp/test", "mixing", selector );
  
  
  
	return 0;
}
