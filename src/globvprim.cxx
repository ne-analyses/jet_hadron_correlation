// testing nick's new correlation implementations

#include "corrParameters.hh"
#include "corrFunctions.hh"

#include <cmath>
#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TRandom3.h"

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
  
  
  TH2D* events = new TH2D("globalvprime", "Global vs Prime", 1500, -0.5, 1499.5, 6000, -0.5, 5999.5 );
  
  TChain* chain = new TChain("JetTree");
  chain = TStarJetPicoUtils::BuildChainFromFileList( "auau_list/grid_AuAuy7HT.list" );
  
  TStarJetPicoReader reader;
  jetHadron::InitReader( reader, chain, "auau", jetHadron::triggerAll, jetHadron::allEvents );
  
  TStarJetPicoEventHeader* header;
  
  try{
    while ( reader.NextEvent() ) {
      
      header = reader.GetEvent()->GetHeader();
      
      events->Fill( header->GetNOfPrimaryTracks(), header->GetNGlobalTracks() );
      
      
    }
  } catch ( std::exception& e) {
    std::cerr << "Caught " << e.what() << std::endl;
    return -1;
  }
  
  TCanvas c1;
  events->Draw("colz");
  c1.SaveAs("test.pdf");
  
  TFile outFile("src/globvprime.root", "RECREATE");
  events->Write();
  outFile.Close();
  
	return 0;
}
