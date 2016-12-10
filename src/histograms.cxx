// ____________________________________________________________________________________
// Class implementation
// corrAnalysis::histograms
// Nick Elsey

#include "corrParameters.hh"
#include "corrFunctions.hh"
#include "histograms.hh"

namespace corrAnalysis {
  
  // These are used by fill functions
  // to check for consistency
  bool histograms::IsPP() {
    
    if ( analysisType.find("pp") != std::string::npos )
    return true;
    
    return false;
  }
  
  bool histograms::IsAuAu() {
    
    if ( analysisType.find("pp") == std::string::npos )
    return true;
    
    return false;
  }
  
  bool histograms::IsDijet() {
    
    if ( analysisType.find("dijet") != std::string::npos )
    return true;
    
    return false;
  }
  
  bool histograms::IsJet() {
    
    if ( analysisType.find("dijet") == std::string::npos && analysisType.find("jet") != std::string::npos )
    return true;
    
    return false;
  }
  
  bool histograms::IsMix() {
    
    if ( analysisType.find("mix") != std::string::npos )
    return true;
    
    return false;
  }
  
  // Used to check if the histograms have been initialized
  // before allowing any filling to stop seg faults
  bool histograms::IsInitialized() {
    if ( initialized ) {
      return true;
    }
    else {
      __ERR("histograms instance not initialized")
      return false;
    }
  }
  
  // Used to build the cent/vz(/aj) arrays used to
  // hold correlations
  void histograms::BuildArrays() {
    
    // split by analysis type
    
    if ( analysisType == "dijet" || analysisType == "ppdijet" || analysisType == "dijetmix" || analysisType == "ppdijetmix" ) {
      
      //now build the full 3D vz/centrality binned histograms
      TH3F* tmpHistLead, * tmpHistSub;
      
      leadingArrays = new TObjArray**[binsAj];
      subleadingArrays = new TObjArray**[binsAj];
      
      for ( int i = 0; i < binsAj; ++i ) {
        leadingArrays[i] = new TObjArray*[binsCentrality];
        subleadingArrays[i] = new TObjArray*[binsCentrality];
        
        
        for ( int j = 0; j < binsCentrality; ++j ) {
          leadingArrays[i][j] = new TObjArray();
          subleadingArrays[i][j] = new TObjArray();
          
          for ( int k = 0; k < binsVz; ++k ) {
            // create unique identifiers for each histogram
            std::stringstream s1, s2, s3;
            s1 << i;
            s2 << j;
            s3 << k;
            
            TString leadName = "lead_aj_";
            TString subName = "sub_aj_";
            if ( analysisType == "dijetmix" || analysisType == "ppdijetmix" ) {
              leadName = "mix_lead_cent_";
              subName = "mix_sub_cent_";
            }
            
            leadName += s1.str() + "_cent_" + s2.str() + "_vz_" + s3.str();
            subName += s1.str() + "_cent_" + s2.str() + "_vz_" + s3.str();
            
            // make each histogram
            tmpHistLead = new TH3F(leadName, leadName+";eta;phi;centrality", binsEta, dEtaLowEdge+etaBinShift, dEtaHighEdge+etaBinShift, binsPhi, phiLowEdge+phiBinShift, phiHighEdge+phiBinShift, binsPt, ptLowEdge, ptHighEdge );
            tmpHistSub = new TH3F(subName, subName+";eta;phi;centrality", binsEta, dEtaLowEdge+etaBinShift, dEtaHighEdge+etaBinShift, binsPhi, phiLowEdge+phiBinShift, phiHighEdge+phiBinShift, binsPt, ptLowEdge, ptHighEdge );
            
            // add to the correct bin
            leadingArrays[i][j]->AddLast( tmpHistLead );
            subleadingArrays[i][j]->AddLast( tmpHistSub );
          }
        }
        std::cout<<"finished i: "<<i<<std::endl;
      }
    }
    
    if ( analysisType == "jet" || analysisType == "ppjet" || analysisType == "jetmix" || analysisType == "ppjetmix" ) {
      
      
      TH3F* tmpHistTrig;
      leadingArrays = new TObjArray**[binsAj];
      
      for ( int i = 0; i < binsAj; ++i ) {
        leadingArrays[i] = new TObjArray*[binsCentrality];
        
        for ( int j = 0; j < binsCentrality; ++j ) {
          leadingArrays[i][j] = new TObjArray();
          
          for ( int k = 0; k < binsVz; ++k ) {
            // create unique identifiers for each histogram
            std::stringstream s1, s2, s3;
            s1 << i;
            s2 << j;
            s3 << k;
            TString leadName = "lead_aj_";
            if ( analysisType == "jetmix" || analysisType == "ppjetmix" ) {
              leadName = "mix_lead_aj_";
            }
            
            leadName += s1.str() + "_cent_" + s2.str() + "_vz_" + s3.str();
            
            tmpHistTrig = new TH3F(leadName, leadName+";eta;phi;centrality", binsEta, dEtaLowEdge+etaBinShift, dEtaHighEdge+etaBinShift, binsPhi, phiLowEdge+phiBinShift, phiHighEdge+phiBinShift, binsPt, ptLowEdge, ptHighEdge );
            
            // add to the correct bin
            leadingArrays[i][j]->AddLast( tmpHistTrig );
          }
        }
      }
    }
  }
  
  // Used internally to pick histogram edges
  // that give a bin centered at zero in correlation plots
  void histograms::FindBinShift() {
    
    // first get the initial variables -
    // bin edges from corrParameters.hh
    // and number of bins from
    // histograms::binsEta, histograms::binsPhi
    
    double phiMin = phiLowEdge;
    double phiMax = phiHighEdge;
    double etaMin = dEtaLowEdge;
    double etaMax = dEtaHighEdge;
    
    double phiBinWidth = ( phiMax - phiMin ) / binsPhi;
    double etaBinWidth = ( etaMax - etaMin ) / binsEta;
    
    phiBinShift = phiMin + ( phiBinWidth / 2.0 );
    etaBinShift = etaMin + ( etaBinWidth / 2.0 );
    
    for ( int i = 0; i < binsPhi; ++i ) {
      if ( fabs(phiBinShift) < fabs(phiBinShift + phiBinWidth) )
      break;
      phiBinShift += phiBinWidth;
    }
    for ( int i = 0; i < binsEta; ++i ) {
      if ( fabs(etaBinShift) < fabs(etaBinShift + etaBinWidth) )
      break;
      etaBinShift += etaBinWidth;
      
    }
    
    if (phiBinShift < 0.001)
    phiBinShift = 0.0;
    if (etaBinShift < 0.001)
    etaBinShift = 0.0;
    
  }

  // Used to find the respective Aj bin
  int histograms::FindAjBin(double aj) {
    double ajBinWidth = ( ajHighEdge - ajLowEdge ) / binsAj;

    for ( int i = 0; i < binsAj; ++i ) {
      if ( aj < ajBinWidth*( i + 1 ) )
        return i; 
    }

  }
  

  histograms::histograms() {
    analysisType = "none";
    initialized = false;
    binsEta = 0;
    binsPhi = 0;
    
    phiBinShift = 0;
    etaBinShift = 0;
    
    hLeadJetPt 	= 0;
    hLeadEtaPhi = 0;
    hSubJetPt 	= 0;
    hSubEtaPhi 	= 0;
    hAssocPt 		= 0;
    hAssocEtaPhi= 0;
    hAjHigh 		= 0;
    hAjLow 			= 0;
    hAjDif      = 0;
    hCentVz			= 0;
    hBinVz			= 0;
    hGRefMult   = 0;
    hVz 				= 0;
    h3DimCorrLead	= 0;
    h3DimCorrSub 	= 0;
    leadingArrays = 0;
    subleadingArrays = 0;
    hAjStruct    = 0;
  }
  
  histograms::histograms( std::string anaType, unsigned tmpBinsEta, unsigned tmpBinsPhi ) {
    analysisType = anaType;
    initialized = false;
    
    binsEta = tmpBinsEta;
    binsPhi = tmpBinsPhi;
    
    phiBinShift = 0;
    etaBinShift = 0;
    
    hLeadJetPt 	= 0;
    hLeadEtaPhi = 0;
    hSubJetPt 	= 0;
    hSubEtaPhi 	= 0;
    hAssocPt 		= 0;
    hAssocEtaPhi= 0;
    hAjHigh 		= 0;
    hAjLow 			= 0;
    hAjDif      = 0;
    hCentVz			= 0;
    hBinVz			= 0;
    hGRefMult   = 0;
    hVz 				= 0;
    h3DimCorrLead	= 0;
    h3DimCorrSub 	= 0;
    leadingArrays = 0;
    subleadingArrays = 0;
    hAjStruct = 0;
  }
  
  histograms::~histograms() {
    Clear();
  }
  
  void histograms::Clear() {
    if ( hLeadJetPt )
    delete hLeadJetPt;
    if ( hLeadEtaPhi )
    delete hLeadEtaPhi;
    if ( hSubJetPt )
    delete hSubJetPt;
    if ( hSubEtaPhi )
    delete hSubEtaPhi;
    if ( hAssocPt )
    delete hAssocPt;
    if ( hAssocEtaPhi )
    delete hAssocEtaPhi;
    if ( hAjLow )
    delete hAjLow;
    if ( hAjHigh )
    delete hAjHigh;
    if ( hAjDif )
    delete hAjDif;
    if ( hCentVz )
    delete hCentVz;
    if ( hBinVz )
    delete hBinVz;
    if ( hGRefMult )
    delete hGRefMult;
    if ( hVz )
    delete hVz;
    if ( h3DimCorrLead )
    delete h3DimCorrLead;
    if ( h3DimCorrSub )
    delete h3DimCorrSub;
    if ( hAjStruct )
    delete hAjStruct;
    
    if ( leadingArrays ) {
      for ( int i = 0; i < binsAj; ++i ) {
        for ( int j = 0; j < binsCentrality; ++j ) {
          leadingArrays[i][j]->Delete();
        }
      }
    }
    if ( subleadingArrays ) {
      for ( int i = 0; i < binsAj; ++i ) {
        for ( int j = 0; j < binsCentrality; ++j ) {
          subleadingArrays[i][j]->Delete();
        }
      }
    }
  }
  
  
  bool histograms::SetAnalysisType( std::string type ) {
    if ( type == analysisType )
    return true;
    
    else if ( type == "dijet" || type == "jet" || type == "ppdijet" || type == "ppjet" || type == "dijetmix" || type == "jetmix" || type == "ppdijetmix" || type == "ppjetmix" ) {
      initialized = false;
      Clear();
      analysisType = type;
      return true;
    }
    
    __ERR("Unknown type")
    return false;
    
  }
  
  int histograms::Init() {
    if ( initialized )
    return 0;
    
    // now set up the phi and eta bin shifts
    FindBinShift();
    
    hCentVz 		= new TH3F( "nevents","Event Count;aj;Centrality;VzBin", binsAj, ajLowEdge, ajHighEdge, binsCentrality, centLowEdge, centHighEdge, binsVz, -0.5, (double) binsVz - 0.5 );
    hGRefMult 	= new TH1D( "grefmultdist", "grefmultdist", 1000, -0.5, 999.5 );
    hVz					 = new TH1D("vzdist", "vzdist", 100, -30, 30);
    
    hLeadJetPt 	= new TH1D( "leadjetpt", "Leading Jet Pt;p_{T}", 80, 0, 80 );
    hLeadEtaPhi = new TH2D( "leadjetetaphi", "Leading Jet Eta Phi;eta;phi", 40, -1, 1, 40, -pi, pi );
    hSubJetPt 	= new TH1D( "subjetpt", "Subleading Jet Pt;p_{T}", 80, 0, 80 );
    hSubEtaPhi 	= new TH2D( "subjetetaphi", "Subleading Jet Eta Phi;eta;phi", 40, -1, 1, 40, -pi, pi );
    
    hAssocPt 		= new TH1D("assocpt", "Associated Track Pt;p_{T}", 80, 0, 12 );
    hAssocEtaPhi= new TH2D("assocetaphi", "Associated Track Eta Phi;#eta;#phi", 40, -1, 1, 40, -pi, pi );
    
    if ( analysisType == "dijet" || analysisType == "ppdijet" ) {
      hAjHigh 		= new TH1D( "ajhigh", "A_{J} High P_{T} Constituents;A_{J};fraction", 30, 0, 0.9 );
      hAjLow 			= new TH1D( "ajlow", "A_{J} Low P_{T} Constituents;A_{J};fraction", 30, 0, 0.9 );
      hAjDif      = new TH3F( "ajdif", "A_{J} difference by A_{J} hard and soft", 30, 0, 1, 30, 0, 1, 30, 0, 1 );
      hAjStruct = new TH3F("hAjStruct", "Aj Centrality Pt", 25, 0, 1, binsCentrality, centLowEdge, centHighEdge, 30, 0, 60 );
    }
    h3DimCorrLead		= new TH3F("leadjetcorr", "Lead Jet - Hadron Correlation;#eta;#phi;p_{T}", binsEta, dEtaLowEdge+etaBinShift, dEtaHighEdge+etaBinShift, binsPhi, phiLowEdge+phiBinShift, phiHighEdge+phiBinShift, binsPt, ptLowEdge, ptHighEdge );
    h3DimCorrSub		= new TH3F("subjetcorr", "Sub Jet - Hadron Correlation;#eta;#phi;p_{T}", binsEta, dEtaLowEdge+etaBinShift, dEtaHighEdge+etaBinShift, binsPhi, phiLowEdge+phiBinShift, phiHighEdge+phiBinShift, binsPt, ptLowEdge, ptHighEdge );
    
    BuildArrays();
    
    initialized = true;

    return 1;
  }
  
  void histograms::Write() {
    
    if ( hCentVz )
    hCentVz->Write();
    if ( hBinVz )
    hBinVz->Write();
    if ( hGRefMult )
    hGRefMult->Write();
    if ( hVz )
    hVz->Write();
    
    if ( hLeadJetPt )
    hLeadJetPt->Write();
    if ( hLeadEtaPhi )
    hLeadEtaPhi->Write();
    if ( hSubJetPt )
    hSubJetPt->Write();
    if ( hSubEtaPhi )
    hSubEtaPhi->Write();
    
    if ( hAssocEtaPhi )
    hAssocEtaPhi->Write();
    if ( hAssocPt )
    hAssocPt->Write();
    
    if ( hAjHigh )
    hAjHigh->Write();
    if ( hAjLow )
    hAjLow->Write();
    if ( hAjDif )
    hAjDif->Write();
    
    if ( h3DimCorrLead )
    h3DimCorrLead->Write();
    if ( h3DimCorrSub )
    h3DimCorrSub->Write();
    if ( hAjStruct )
    hAjStruct->Write();
    
    for ( int i = 0; i < binsAj; ++i ) {
      if (!leadingArrays || !leadingArrays[i])
        continue;
      for ( int j = 0; j < binsCentrality; ++j )
        if ( !leadingArrays[i][j] ) {
          continue;
        for ( int k = 0; k < binsVz; ++k ) {
          leadingArrays[i][j][k].Write();
          if ( IsDijet() )
            subleadingArrays[i][j][k].Write();
        }
      }
    }
  }
  
  
  // --------------------------- Histogram Filling Functions ------------------------------- //
  bool histograms::CountEvent( int vzbin, int centrality, double aj ) {
    if ( !IsInitialized() ) { return false; }
    
    hCentVz->Fill( aj, centrality,(double) vzbin );
    return true;
  }
  
  bool histograms::FillGRefMult( int gRefMult ) {
    if ( !IsInitialized() ) { return false; }
    
    hGRefMult->Fill( gRefMult );
    return true;
  }
  
  bool histograms::FillVz( double vz ) {
    if ( !IsInitialized() ) { return false; }
    
    hVz->Fill( vz );
    return true;
  }
  
  bool histograms::FillAjHigh( double aj ) {
    if ( !IsInitialized() ) { return false; }
    
    hAjHigh->Fill( aj );
    return true;
  }
  
  bool histograms::FillAjLow( double aj ) {
    if ( !IsInitialized() ) { return false; }
    
    hAjLow->Fill( aj );
    return true;
  }
  
  bool histograms::FillAjDif( double high, double low ) {
    if ( !IsInitialized() ) { return false; }
    
    hAjDif->Fill( high, low, fabs(high - low) );
    return true;
  }
  
  bool histograms::FillJetPt( double pt ) {
    if ( !IsInitialized() ) { return false; }
    
    hLeadJetPt->Fill( pt );
    return true;
  }
  
  bool histograms::FillJetEtaPhi( double eta, double phi ) {
    if ( !IsInitialized() ) { return false; }
    
    hLeadEtaPhi->Fill( eta, phi );
    return true;
  }
  
  bool histograms::FillCorrelation( double dEta,  double dPhi, double assocPt, double weight, int vzBin, int centBin ) {
    if ( !IsInitialized() ) { return false; }
    
    if ( dPhi < phiLowEdge+phiBinShift)
    dPhi += 2.0*pi;
    
    h3DimCorrLead->Fill( dEta, dPhi, assocPt, weight );
      
    // now do the bin-divided fill
    TH3F* tmpHist = (TH3F*) leadingArrays[0][centBin]->At(vzBin);
    tmpHist->Fill( dEta, dPhi, assocPt, weight );
    return true;
  }
  
  bool histograms::FillLeadJetPt( double pt ) {
    if ( !IsInitialized() ) { return false; }
    
    hLeadJetPt->Fill( pt );
    return true;
  }
  
  bool histograms::FillSubJetPt( double pt ) {
    if ( !IsInitialized() ) { return false; }
    
    hSubJetPt->Fill( pt );
    return true;
  }
  
  bool histograms::FillLeadEtaPhi( double eta, double phi ) {
    if ( !IsInitialized() ) { return false; }
    
    hLeadEtaPhi->Fill( eta, phi );
    return true;
  }
  
  bool histograms::FillSubEtaPhi( double eta, double phi ) {
    if ( !IsInitialized() )
      return false;
    
    hSubEtaPhi->Fill( eta, phi );
    return true;
  }
  
  bool histograms::FillCorrelationLead( double dEta, double dPhi, double assocPt, double weight, double aj, int vzBin, int centBin ) {
    if ( !IsInitialized() ) { return false; }
    
    if ( dPhi < phiLowEdge+phiBinShift )
    dPhi += 2.0*pi;
    
    // find aj bin, if applicable
    int binAj = FindAjBin( aj );
    
    h3DimCorrLead->Fill( dEta, dPhi, assocPt, weight );
    
    // now do the bin-divided fill
    TH3F* tmpHist = (TH3F*) leadingArrays[binAj][centBin]->At(vzBin);
    tmpHist->Fill( dEta, dPhi, assocPt, weight );
      
    return true;
  }
  
  bool histograms::FillCorrelationSub( double dEta, double dPhi, double assocPt, double weight, double aj, int vzBin, int centBin ) {
    if ( !IsInitialized() ) { return false; }
    
    if ( dPhi < phiLowEdge+phiBinShift )
    dPhi += 2.0*pi;
    
    // find aj bin, if applicable
    int binAj = FindAjBin( aj );
    
    h3DimCorrSub->Fill( dEta, dPhi, assocPt, weight );
      
    // now do the bin-divided fill
    TH3F* tmpHist = (TH3F*) subleadingArrays[binAj][centBin]->At(vzBin);
    tmpHist->Fill( dEta, dPhi, assocPt, weight );
      
    return true;
  }
  
  bool histograms::FillAssocPt( double pt ) {
    if ( !IsInitialized() ) { return false; }
    
    hAssocPt->Fill( pt );
    return true;
  }
  
  bool histograms::FillAssocEtaPhi( double eta, double phi ) {
    if ( !IsInitialized() ) { return false; }
    
    hAssocEtaPhi->Fill( eta, phi );
    return true;
  }
  
  // Looking for correlations between Aj, Cent, and Pt
  bool histograms::FillAjStruct( double aj, double pt, int centrality ) {
    if ( !IsInitialized() ) { return false; }
    
      
    hAjStruct->Fill( aj, centrality, pt );
    return true;
  }

} // end namespace
