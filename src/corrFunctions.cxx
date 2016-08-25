// Implementation of basic functions used
// in the dijet-hadron and jet-hadron
// correlation analysis
// Nick Elsey

#include "corrFunctions.hh"
#include "corrParameters.hh"

namespace corrAnalysis {
	
	// -------------------------
	// IO/OS Manip functionality
	// -------------------------
	// Used to understand which format of input file is being used
	// ( .root file, .txt, .list, etc )
	// ---------------------------------------------------------------------
	bool HasEnding (std::string const &full_string, std::string const &ending) {
		if (full_string.length() >= ending.length()) {
			return (0 == full_string.compare (full_string.length() - ending.length(), ending.length(), ending) );
		} else {
			return false;
		}
	}

	
	// Used to find current working directory
	// Trying to make this relatively machine independent
	std::string getPWD() {
		char buffer[255];
		char *answer = getcwd( buffer, sizeof(buffer) );
		std::string s_cwd;
		if (answer) { s_cwd = answer; }
		return s_cwd;
	}
	
	// -------------------------
	// Analysis functionaliy
	// -------------------------

	// Function to return DPhi
	// Not necessary in current version
	// Analysis uses fastjet::PseudoJet::delta_phi_to
	// ------------------------------------------------
	double GetdPhi ( double trigPhi, double assocPhi ){
		if ( trigPhi < -1*pi ) trigPhi += ( 2 * pi );
		else if ( trigPhi > pi ) trigPhi -= ( 2 * pi );
		if ( assocPhi < -1 * pi ) assocPhi += ( 2 * pi );
		else if ( assocPhi > pi ) assocPhi -= ( 2 * pi );
		double dphi= assocPhi - trigPhi;
		if ( dphi < -1 * pi ) dphi += ( 2 * pi );
		else if ( dphi > pi ) dphi -= ( 2 * pi );
		return dphi;
	}

	// Calculate Aj fraction
	// Used in dijet analysis to compare
	// To the Aj results for consistency
	// Checks That there are the correct # of Jets, then returns A_j
	// -------------------------------------------------------------
	double CalcAj ( std::vector<fastjet::PseudoJet>& jets ){
		if ( jets.size()!=2 ){
			throw ( -1 );
			return -1e10;
		}
		return fabs (( jets.at(0).pt()-jets.at(1).pt() ) / ( jets.at(0).pt()+jets.at(1).pt() ));
	}
  
  // Used to get reference centrality from gRefMult
  // Uses definitions in corrParameters.hh
  // ----------------------------------------------
  int GetReferenceCentrality( int gRefMult ) {
    for ( int i = 8; i >= 0;  --i )
      if ( gRefMult >= y7RefMultCent[i] ) {
        return i;
      }
    // in case there is strange input
    __ERR("Reference Centrality not found")
    return -1;
  }
  
  // Used to get the inverse definition of reference centality
  // ( 8->0, 7->1, etc)
  // ----------------------------------------------
  int GetReferenceCentralityAlt( int RefCent ) {
    return abs( 8 - RefCent );
  }
  
  // Used to find the Vz Bin - definitions for Vz
  // Binning found in corrParameters
  // ----------------------------------------------
  int GetVzBin( double Vz ) {
    int VzBin = -1;
    // Check the boundaries
    if ( Vz > vzHighEdge || Vz <= vzLowEdge )
      return -1;
    
    // Loop and check each bin
    for ( int i = binsVz-1; i >= 0; --i ) {
      if ( Vz > ( vzLowEdge + dVz*i ) ) {
        return i;
      }
    }
    __ERR("There is a problem with VzBin finding")
    return VzBin;
  }
  
  // Fills our working container after converting TStarJetVectors into PseudoJets
  // Also makes sure to empty the containers if they are full
  //
  // ------------------------------------------------------------------------------
  void ConvertTStarJetVector( TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles, bool ClearVector ) {
    // Empty the container
    // if called for
    if ( ClearVector )
    	particles.clear();
    
    // Transform TStarJetVectors into (FastJet) PseudoJets
    // ---------------------------------------------------
    TStarJetVector* sv;
    for ( int i=0; i < container->GetEntries() ; ++i ){
      sv = container->Get(i);
      
      fastjet::PseudoJet tmpPJ = fastjet::PseudoJet( *sv );
      tmpPJ.set_user_index( sv->GetCharge() );
      particles.push_back( tmpPJ );
      
      
    }
  }
  
  // finds the triggers and saves them, if requireTrigger == True
  void GetTriggers( bool requireTrigger, TClonesArray* triggerObjs, std::vector<fastjet::PseudoJet> & triggers ) {
    
    // empty the container
    triggers.clear();
    // If we use triggers, pull every trigger from the event
    // And convert the high tower triggers
    if ( requireTrigger ) {
      TIter nextTrigger(triggerObjs);
      TStarJetPicoTriggerInfo* trigger = 0;
      while ( ( trigger = (TStarJetPicoTriggerInfo*) nextTrigger() ) ) {
        if ( trigger->GetTriggerFlag() == 1 ) {
          fastjet::PseudoJet tmpTrig;
          tmpTrig.reset_PtYPhiM(0.1, trigger->GetEta(), trigger->GetPhi(), 0 );
          triggers.push_back( tmpTrig );
        }
      }
    }
  }
  
  // for the pp data where the trigger objects dont seem to be working
  void GetTriggersPP( bool requireTrigger, std::vector<fastjet::PseudoJet> ppParticles, std::vector<fastjet::PseudoJet>& triggers ) {
    // empty the container
    triggers.clear();
    
    // if we're using triggers, run over all towers and get any with E > triggerThreshold
    if ( requireTrigger ) {
      for ( int i = 0; i < ppParticles.size(); ++i ) {
        fastjet::PseudoJet tmpParticle = ppParticles[i];
        if ( tmpParticle.pt() > triggerThreshold )
          triggers.push_back( tmpParticle );
      
      }
    }
  }

	// ----------------------
	// Verbose output scripts
	// ----------------------

	// Settings summary, called before dijet event loop
	// -------------------------------------------------------
	void BeginSummaryDijet ( double jetRadius, double leadJetPtMin, double subLeadJetPtMin, double jetMaxPt, double hardJetConstPt, double softJetConstPt, int nVzBins, double vzRange, std::string dijetFile, std::string corrFile ) {
		std::cout<<" ------- SUMMARY OF SETTINGS ------"<<std::endl;
		std::cout<<" Jetfinding Algorithm: Anti-kt"<<std::endl;
		std::cout<<" Resolution Parameter: "<<jetRadius<<std::endl;
		std::cout<<" Leading Jet minimum pt: "<<leadJetPtMin<<std::endl;
		std::cout<<" SubLeading Jet minimum pt: "<<subLeadJetPtMin<<std::endl;
		std::cout<<" Jet Maximum Pt: "<<jetMaxPt<<std::endl;
		std::cout<<" Hard Jet Constituent minimum pt: "<< hardJetConstPt<<std::endl;
		std::cout<<" Soft Jet Constituent minimum pt: "<< softJetConstPt<<std::endl;
		std::cout<<" Number of bins in Vertex Z position: "<<nVzBins<<std::endl;
		std::cout<<" Vertex Z range (symmetric about 0): "<<vzRange<<std::endl;
		std::cout<<" Outputting dijet tree to: "<<dijetFile<<std::endl;
		std::cout<<" Outputting histograms to: "<<corrFile<<std::endl;
		std::cout<<" ---------- END SUMMARY ----------"<<std::endl;
	}

	// Settings summary, called before dijet event loop
	// -------------------------------------------------------
	void BeginSummaryJet ( double jetRadius, double jetPtMin, double jetPtMax, double jetConstPt, int nVzBins, double vzRange, std::string jetFile, std::string corrFile ) {
	std::cout<<" ------- SUMMARY OF SETTINGS ------"<<std::endl;
		std::cout<<" Jetfinding Algorithm: Anti-kt"<<std::endl;
		std::cout<<" Resolution Parameter: "<<jetRadius<<std::endl;
		std::cout<<" Leading Jet minimum pt: "<<jetPtMin<<std::endl;
		std::cout<<" Jet Constituent minimum pt: "<< jetConstPt<<std::endl;
		std::cout<<" Number of bins in Vertex Z position: "<<nVzBins<<std::endl;
		std::cout<<" Vertex Z range (symmetric about 0): "<<vzRange<<std::endl;
		std::cout<<" Outputting jet tree to: "<<jetFile<<std::endl;
		std::cout<<" Outputting histograms to: "<<corrFile<<std::endl;
		std::cout<<" ---------- END SUMMARY ----------"<<std::endl;
	}

	// called after dijet correlation event loop is complete
	// ---------------------------------------------------------------------
	void EndSummaryDijet ( int ntotal, int nviable, int nused, double time ) {
		std::cout<<"  ----------------- SUMMARY ----------------- "<<std::endl;
		std::cout<<"  Processed "<< ntotal <<" events in "<< time << " seconds."<<std::endl;
		std::cout<<"  Of these, "<< nviable << " produced hard dijet pairs,"<<std::endl;
		std::cout<<"  Which corresponds to "<< 100.0*(double) nviable / (double) ntotal <<"%"<<std::endl;
		std::cout<<"  Chance per event to find a hard dijet pair"<<std::endl;
		std::cout<<"  Of these "<< nviable <<" hard dijets, "<< nused <<" produced full dijets that were used"<<std::endl;
		std::cout<<"  for correlation."<<std::endl;
		std::cout<<"  Overall Efficiency: "<< time/ (double) nused <<" seconds per dijet"<<std::endl;
	}

	// Called after jet correlation event loop is complete
	// ---------------------------------------------------------------------
	void EndSummaryJet ( int ntotal, int nused, double time ) {
		std::cout<<"  ----------------- SUMMARY ----------------- "<<std::endl;
		std::cout<<"  Processed "<< ntotal <<" events in "<< time << " seconds."<<std::endl;
		std::cout<<"  Of these, "<< nused << " produced useable leading jets"<<std::endl;
		std::cout<<"  Which corresponds to "<< 100.0*(double) nused / (double) ntotal <<"%"<<std::endl;
		std::cout<<"  Chance per event to find a leading jet"<<std::endl;
		std::cout<<"  for correlation."<<std::endl;
		std::cout<<"  Overall Efficiency: "<< time/ (double) nused <<" seconds per jet"<<std::endl;
	}
	
  // Used to initialized the reader - will set the event cuts,
  // Tower cuts, track cuts and hadronic correction
  // ---------------------------------------------------------------------
  void InitReader( TStarJetPicoReader & reader, TChain* chain, std::string collisionType, std::string triggerString, int nEvents ) {
    
    // First tolower() on the analysisType
    // shouldnt be necessary....
    std::transform(collisionType.begin(), collisionType.end(), collisionType.begin(), ::tolower);
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader.SetApplyFractionHadronicCorrection( hadronicCorrection );
    reader.SetFractionHadronicCorrection( hadronicCorrectionFraction );
    reader.SetRejectTowerElectrons( kFALSE );
    
    // Event and track selection
    // -------------------------
    
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetTriggerSelection( triggerString.c_str() );
    evCuts->SetVertexZCut ( vertexZCut );
    evCuts->SetMaxEventPtCut( eventPtCut );
    evCuts->SetMaxEventEtCut( eventEtCut );
    evCuts->SetVertexZDiffCut( vertexZDiffCut );
    if ( collisionType == "auau" ) {
    	evCuts->SetRefMultCut ( y7RefMultCut );
    }
    else if ( collisionType == "pp" ) {
      evCuts->SetRefMultCut( 0 );
    }
    else
      __ERR("unknown collision system")
    
    // Tracks cuts
    TStarJetPicoTrackCuts* trackCuts = reader.GetTrackCuts();
    trackCuts->SetDCACut( DCACut );
    trackCuts->SetMinNFitPointsCut( minFitPoints );
    trackCuts->SetFitOverMaxPointsCut( minFitFrac );
    trackCuts->SetMaxPtCut ( trackPtCut );
    
    
    std::cout << "Using these track cuts:" << std::endl;
    std::cout << " dca : " << trackCuts->GetDCACut(  ) << std::endl;
    std::cout << " nfit : " <<   trackCuts->GetMinNFitPointsCut( ) << std::endl;
    std::cout << " nfitratio : " <<   trackCuts->GetFitOverMaxPointsCut( ) << std::endl;
    std::cout << " maxpt : " << trackCuts->GetMaxPtCut (  ) << std::endl;
    
    // Towers
    TStarJetPicoTowerCuts* towerCuts = reader.GetTowerCuts();
    towerCuts->SetMaxEtCut( towerEtCut );
    if ( collisionType == "auau" ) {
    	towerCuts->AddBadTowers( y7AuAuTowerList.c_str() );
    }
    else if ( collisionType == "pp" ) {
      towerCuts->AddBadTowers( y6PPTowerList.c_str() );
    }
    else
      __ERR("unknown collision system")
    std::cout << "Using these tower cuts:" << std::endl;
    std::cout << "  GetMaxEtCut = " << towerCuts->GetMaxEtCut() << std::endl;
    std::cout << "  Gety8PythiaCut = " << towerCuts->Gety8PythiaCut() << std::endl;
    
    // V0s: Turn off
    reader.SetProcessV0s(false);
    
    // Initialize the reader
    reader.Init( nEvents ); //runs through all events with -1

  }
  
  // Use this to decide if there are 2 dijets for dijet analysis
  // in the proper pt ranges, and if they're back to back
  // Or for jet analysis if there is a single jet
  bool CheckHardCandidateJets( std::string analysisType, std::vector<fastjet::PseudoJet> & HiResult, double leadJetPtMin, double subJetPtMin ) {
    if ( analysisType == "dijet" || analysisType == "ppdijet" ) {
      if ( HiResult.size() < 2 ) 									{ return false; }
      if ( HiResult.at(0).pt() < leadJetPtMin )   { return false; }
      if ( HiResult.at(1).pt() < subJetPtMin ) 		{ return false; }
      if ( fabs( fabs( HiResult.at(0).delta_phi_to( HiResult.at(1) ) ) - TMath::Pi() ) > jetDPhiCut  ) { return false; }
    }
    else if ( analysisType == "jet" || analysisType == "ppjet" ) {
      if ( HiResult.size() == 0 )									{ return false;	}
    }
    else {
      __ERR("Unrecognized analysis type")
      throw(-1);
    }
    return true;
  }
  
  // Use this to select either one or two jets
  // Depending on the analysis type
  // -----------------------------------------
  std::vector<fastjet::PseudoJet> BuildHardJets( std::string analysisType, std::vector<fastjet::PseudoJet> & HiResult ) {
    std::vector<fastjet::PseudoJet> tmpJets;
    if ( analysisType == "dijet" || analysisType == "ppdijet" ) {
      if ( HiResult.size() < 2 ) {
        __ERR("Dijet analysis but less than two candidate jets")
        throw(-1);
      }
      tmpJets.push_back( HiResult.at(0) );
      tmpJets.push_back( HiResult.at(1) );
    }
    else if ( analysisType == "jet" || analysisType == "ppjet" ) {
      if ( HiResult.size() < 1 ) {
        __ERR("Jet analysis but less than one candidate jet")
        throw(-1);
      }
      return HiResult;
    }
    else {
      __ERR("Undefined analysis type")
      throw(-1);
    }
    
    
    return tmpJets;
  }
  
  //
  std::vector<fastjet::PseudoJet> BuildMatchedJets( std::string analysisType, std::vector<fastjet::PseudoJet> & hardJets, std::vector<fastjet::PseudoJet> & LoResult, bool requireTrigger, std::vector<fastjet::PseudoJet> & triggers, double jetRadius ) {
    if ( analysisType == "dijet" || analysisType == "ppdijet" ) {
      
      // make sure the input makes sense
      if ( hardJets.size() != 2 ) {
        __ERR("needs two hard jets for dijet analysis")
        throw( -1 );
      }
      
      // Match hard jets to full jets
      // match the leading jet
      fastjet::Selector selectMatchedLead = fastjet::SelectorCircle( jetRadius );
      selectMatchedLead.set_reference( hardJets.at(0) );
      std::vector<fastjet::PseudoJet> matchedToLead = sorted_by_pt( selectMatchedLead( LoResult ));
      
      // match the subleading jet
      fastjet::Selector selectMatchedSub = fastjet::SelectorCircle( jetRadius );
      selectMatchedSub.set_reference( hardJets.at(1) );
      std::vector<fastjet::PseudoJet> matchedToSub = sorted_by_pt( selectMatchedSub( LoResult ));
      
      if ( matchedToLead.size() == 0 || matchedToSub.size() == 0 ) {
        __OUT("Couldn't match hard and soft jets")
        return std::vector<fastjet::PseudoJet>();
      }
      
      std::vector<fastjet::PseudoJet> matchedToDijet;
      matchedToDijet.push_back( matchedToLead.at(0) );
      matchedToDijet.push_back( matchedToSub.at(0) );
      
      // if we need a trigger, match jets
      // if the trigger is in the subleading jet,
      // make sub jet the leading jet
      // otherwise, return the dijets without matching
      if ( requireTrigger ) {
        bool matchedLeadTrigger = false;
        bool matchedSubTrigger = false;
        for ( int i = 0; i < triggers.size(); ++i ) {
          if ( matchedToDijet.at(0).delta_R( triggers.at(i) ) < jetRadius )
            matchedLeadTrigger = true;
          else if ( matchedToDijet.at(1).delta_R( triggers.at(i) ) < jetRadius )
            matchedSubTrigger = true;
        }
        
        // now return in proper order
        // if subjet was matched but not the leading jet,
        // reverse
        if ( matchedSubTrigger && !matchedLeadTrigger ) {
          std::vector<fastjet::PseudoJet> reversedJets;
          reversedJets.push_back( matchedToDijet.at(1) );
          reversedJets.push_back( matchedToDijet.at(0) );
          return reversedJets;
        }
        else if ( matchedLeadTrigger )
          return matchedToDijet;
       else
         return std::vector<fastjet::PseudoJet>();
      }
      
      return matchedToDijet;
      
    }
    else if ( analysisType == "jet" || analysisType == "ppjet" ) {
      if ( hardJets.size() == 0 ) {
        __ERR("need at least one jet for jet analysis")
        throw(-1);
      }
			
      // Jet analysis uses hard jet so no need to match
      // Check if it has to be matched to HT trigger
      if ( requireTrigger ) {
        // build selector to match triggers to the jet
        fastjet::Selector selectMatchedTrigger = fastjet::SelectorCircle( jetRadius );
        
        for ( int j = 0; j < triggers.size(); ++ j ) {
	        selectMatchedTrigger.set_reference( triggers.at(j) );
					
          std::vector<fastjet::PseudoJet> matchedToJet = sorted_by_pt(selectMatchedTrigger( hardJets ) );
          if ( matchedToJet.size() > 0 ) {
            matchedToJet.resize(1);
            return matchedToJet;
          }
        	
        }
        // if we didnt find a matched trigger, return empty
        return std::vector<fastjet::PseudoJet>();
      }
      else {
        // if we dont match to triggers, simply return
        // highest pt
        std::vector<fastjet::PseudoJet> matchedToJet;
        matchedToJet.push_back( hardJets.at(0) );
        return matchedToJet;
      }
    }
    else {
    	__ERR("Undefined analysis type")
    	throw( -1 );
    }
    
    __ERR("Something weird happened")
    throw( -1 );
  }
  
  // Used to correlate jets and their charged associated particles
  // Checks to make sure the efficiency is sane
  // First, to check that the track makes all cuts
  bool useTrack( fastjet::PseudoJet& assocTrack, double efficiency ) {
    // Check to make sure the its a charged track within our eta acceptance
    if ( assocTrack.eta() > maxTrackRap )			{ return false; }
    if ( assocTrack.user_index() == 0 )  			{ return false; }
    
    // Check to make sure the efficiency is not crazy
    // ( the parameterization isnt perfect, about 3.5% of tracks return nonsense efficiencies )
    if ( efficiency <= 0.0 )        { return false;  }
    if ( efficiency > 1.0 ) 				{ return false;  }
    
    return true;
  }
  
  bool correlateLeading( std::string analysisType, int vzBin, int centBin, histograms* histogram, fastjet::PseudoJet& leadJet, fastjet::PseudoJet& assocTrack, double efficiency ) {
    
    // check if track is ok
    if ( !useTrack( assocTrack, efficiency ) )
      return false;
    
    // track can be used, so get dPhi and dEta
    double deltaEta = leadJet.eta() - assocTrack.eta();
    double deltaPhi = leadJet.delta_phi_to( assocTrack );
    double assocPt =	assocTrack.pt();
    double weight = 1.0/efficiency;
    
    // Fill some debug info
    histogram->FillAssocEtaPhi( assocTrack.eta(), assocTrack.phi_std() );
    histogram->FillAssocPt( assocPt );
    
    // now fill the histograms
    histogram->FillCorrelationLead( deltaEta, deltaPhi, assocPt, weight, vzBin, centBin );
    
    return true;
  }
  
  bool correlateSubleading( std::string analysisType, int vzBin, int centBin, histograms* histogram, fastjet::PseudoJet& subJet, fastjet::PseudoJet& assocTrack, double efficiency ) {
    
    // check if track is ok
    if ( !useTrack( assocTrack, efficiency ) )
      return false;
    
    double deltaEta = subJet.eta() - assocTrack.eta();
    double deltaPhi = subJet.delta_phi_to( assocTrack );
    double assocPt =  assocTrack.pt();
    double weight = 1.0/efficiency;
    
    // Fill some debug info
    histogram->FillAssocEtaPhi( assocTrack.eta(), assocTrack.phi_std() );
    histogram->FillAssocPt( assocPt );
    
    // now fill the histograms
    histogram->FillCorrelationSub( deltaEta, deltaPhi, assocPt, weight, vzBin, centBin );
    
    return true;
  }
  
  bool correlateTrigger( std::string analysisType, int vzBin, int centBin, histograms* histogram, fastjet::PseudoJet& triggerJet, fastjet::PseudoJet& assocTrack, double efficiency ) {
    
    // check if track is ok
    if ( !useTrack( assocTrack, efficiency ) )
      return false;
    
    double deltaEta = triggerJet.eta() - assocTrack.eta();
    double deltaPhi = triggerJet.delta_phi_to( assocTrack );
    double assocPt =  assocTrack.pt();
    double weight = 1.0/efficiency;
    
    // Fill some debug info
    histogram->FillAssocEtaPhi( assocTrack.eta(), assocTrack.phi_std() );
    histogram->FillAssocPt( assocPt );
    
    // now fill the histograms
    histogram->FillCorrelation( deltaEta, deltaPhi, assocPt, weight, vzBin, centBin );
    
    return true;
  }

  
	
	// JetFinding and Selector functionality
	
  // Build the analysis jet definition using anti-kt
	// -----------------------------------------------
  fastjet::JetDefinition AnalysisJetDefinition( double R ) {
    return fastjet::JetDefinition( fastjet::antikt_algorithm, R );
  }
  
	// Build the background jet definition using kt
  // --------------------------------------------
  fastjet::JetDefinition BackgroundJetDefinition( double R ) {
    return fastjet::JetDefinition( fastjet::kt_algorithm, R );
  }
	
  // Selector for low Pt constituents -
  // clustering with background subtraction
  // --------------------------------------------
  fastjet::Selector	SelectLowPtConstituents( double trackMaxEta, double trackMinPtLow ) {
    return fastjet::SelectorAbsRapMax( trackMaxEta ) * fastjet::SelectorPtMin( trackMinPtLow );
  }
  
  // Selector for high Pt constituents -
  // clustering without background subtraction initially
  // --------------------------------------------
  fastjet::Selector SelectHighPtConstituents( double trackMaxEta, double trackMinPtHigh ) {
    return fastjet::SelectorAbsRapMax( trackMaxEta ) * fastjet::SelectorPtMin( trackMinPtHigh );
  }
	
  // Selects candidate jets with | eta | < trackMaxEta - jetRadius
  // to keep jets away from the edges of our acceptance
  // --------------------------------------------
  fastjet::Selector SelectJetCandidates( double trackMaxEta, double jetRadius, double jetMinPt, double jetMaxPt ) {
    return fastjet::SelectorPtMin( jetMinPt ) * fastjet::SelectorPtMax( jetMaxPt ) * fastjet::SelectorAbsRapMax( trackMaxEta - jetRadius );
  }
  
  // Ghosted area spec, used for background
  // estimation and subtraction
  // --------------------------------------
  fastjet::GhostedAreaSpec GhostedArea( double trackMaxEta, double jetRadius ) {
    double ghostMaxRap = (trackMaxEta - jetRadius) + 2.0*jetRadius;
    return fastjet::GhostedAreaSpec( ghostMaxRap, ghostRepeat, ghostArea );
  }
  
  // Area definition used for background estimation
  // and subtraction
  // ----------------------------------------------
  fastjet::AreaDefinition  AreaDefinition( fastjet::GhostedAreaSpec ghostAreaSpec ) {
    return fastjet::AreaDefinition( fastjet::active_area_explicit_ghosts, ghostAreaSpec );
  }
  
  // Selector used to estimate background
  fastjet::Selector SelectBkgEstimator( double maxTrackRap, double jetRadius ) {
    return fastjet::SelectorAbsRapMax( maxTrackRap - jetRadius ) * (!fastjet::SelectorNHardest(2));
  }
	
	// ____________________________________________________________________________________
	// Class implementation
	// corrAnalysis::histograms
	
	histograms::histograms() {
		analysisType = "none";
		initialized = false;
		
		hLeadJetPt 	= 0;
		hLeadEtaPhi = 0;
		hSubJetPt 	= 0;
		hSubEtaPhi 	= 0;
    hAssocPt 		= 0;
    hAssocEtaPhi= 0;
		hAjHigh 		= 0;
		hAjLow 			= 0;
		hCentVz			= 0;
		hBinVz			= 0;
		hGRefMult   = 0;
		hVz 				= 0;
		h3DimCorrLead	= 0;
		h3DimCorrSub 	= 0;
    leadingArrays = 0;
    subleadingArrays = 0;
	}
	
	histograms::histograms( std::string anaType ) {
		analysisType = anaType;
		initialized = false;
		hLeadJetPt 	= 0;
		hLeadEtaPhi = 0;
		hSubJetPt 	= 0;
		hSubEtaPhi 	= 0;
    hAssocPt 		= 0;
    hAssocEtaPhi= 0;
		hAjHigh 		= 0;
		hAjLow 			= 0;
		hCentVz			= 0;
		hBinVz			= 0;
		hGRefMult   = 0;
		hVz 				= 0;
		h3DimCorrLead	= 0;
		h3DimCorrSub 	= 0;
    leadingArrays = 0;
    subleadingArrays = 0;
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
    
    if ( leadingArrays ) {
      for ( int i = 0; i < binsCentrality; ++i ) {
        leadingArrays[i]->Delete();
      }
    }
    if ( subleadingArrays ) {
      for ( int i = 0; i < binsCentrality; ++i ) {
        subleadingArrays[i]->Delete();
      }
    }
	}
	
	
	bool histograms::SetAnalysisType( std::string type ) {
		if ( type == analysisType )
			return true;
		
		else if ( type == "dijet" || type == "jet" || type == "ppdijet" || type == "ppjet" ) {
			initialized = false;
			Clear();
			analysisType = type;
			return true;
		}
		
		std::cerr << "histograms::SetAnalysisType(): Error: Unknown type" << std::endl;
		return false;
		
	}
	
	int histograms::Init() {
		if ( initialized )
			return 0;
		
		if ( analysisType == "dijet" ) {
			hCentVz 		= new TH2D( "nevents","Event Count;Centrality;VzBin", binsCentrality, centLowEdge, centHighEdge, binsVz, -0.5, (double) binsVz - 0.5 );
			hGRefMult 	= new TH1D( "grefmultdist", "grefmultdist", 1000, -0.5, 999.5 );
			hVz					 = new TH1D("vzdist", "vzdist", 100, -30, 30);
			
			hLeadJetPt 	= new TH1D( "leadjetpt", "Leading Jet Pt;p_{T}", 80, 0, 80 );
			hLeadEtaPhi = new TH2D( "leadjetetaphi", "Leading Jet Eta Phi;eta;phi", 40, -1, 1, 40, -pi, pi );
			hSubJetPt 	= new TH1D( "subjetpt", "Subleading Jet Pt;p_{T}", 80, 0, 80 );
			hSubEtaPhi 	= new TH2D( "subjetetaphi", "Subleading Jet Eta Phi;eta;phi", 40, -1, 1, 40, -pi, pi );
      
      hAssocPt 		= new TH1D("assocpt", "Associated Track Pt;p_{T}", 80, 0, 12 );
      hAssocEtaPhi= new TH2D("assocetaphi", "Associated Track Eta Phi;#eta;#phi", 40, -1, 1, 40, -pi, pi );
			
			hAjHigh 		= new TH1D( "ajhigh", "A_{J} High P_{T} Constituents;A_{J};fraction", 30, 0, 0.9 );
			hAjLow 			= new TH1D( "ajlow", "A_{J} Low P_{T} Constituents;A_{J};fraction", 30, 0, 0.9 );
			
			h3DimCorrLead		= new TH3D("leadjetcorr", "Lead Jet - Hadron Correlation;#eta;#phi;p_{T}", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
			h3DimCorrSub		= new TH3D("subjetcorr", "Sub Jet - Hadron Correlation;#eta;#phi;p_{T}", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
      
      // now build the full 3D vz/centrality binned histograms
      TH3D* tmpHistLead, * tmpHistSub;
      leadingArrays = new TObjArray*[binsCentrality];
      subleadingArrays = new TObjArray*[binsCentrality];
      for ( int i = 0; i < binsCentrality; ++i ) {
        leadingArrays[i] = new TObjArray();
        subleadingArrays[i] = new TObjArray;
        for ( int j = 0; j < binsVz; ++j ) {
          // Unique name for each TH3D
          std::stringstream s1, s2;
          s1 << i;
          s2 << j;
          TString leadName = "lead_cent_";
          TString subName = "sub_cent_";
          leadName += s1.str() + "_vz_" + s2.str();
          subName += s1.str() + "_vz_" + s2.str();
          
          // make each histogram
          tmpHistLead = new TH3D(leadName, leadName+";eta;phi;centrality", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
          tmpHistSub = new TH3D(subName, subName+";eta;phi;centrality", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
          
          // add to the correct bin
          leadingArrays[i]->AddLast( tmpHistLead );
          subleadingArrays[i]->AddLast( tmpHistSub );

        }
      }
			
			initialized = true;
			return 0;
		}
		else if ( analysisType == "jet" ) {
			hCentVz 		= new TH2D( "nevents","Event Count;Centrality;VzBin", binsCentrality, centLowEdge, centHighEdge, binsVz, -0.5, (double) binsVz - 0.5 );
			hGRefMult 	= new TH1D( "grefmultdist", "grefmultdist", 1000, -0.5, 999.5 );
			hVz					 = new TH1D("vzdist", "vzdist", 100, -30, 30);
			
			hLeadJetPt 	= new TH1D( "triggerjetpt", "Trigger Jet Pt;p_{T}", 80, 0, 80 );
			hLeadEtaPhi = new TH2D( "triggerjetetaphi", "Trigger Jet Eta Phi;eta;phi", 40, -1, 1, 40, -pi, pi );
      
      hAssocPt 		= new TH1D("assocpt", "Associated Track Pt;p_{T}", 80, 0, 12 );
      hAssocEtaPhi= new TH2D("assocetaphi", "Associated Track Eta Phi;#eta;#phi", 40, -1, 1, 40, -pi, pi );
			
			h3DimCorrLead		= new TH3D("triggerjetcorr", "Jet - Hadron Correlation;#eta;#phi;p_{T}", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
      
      // now build the full 3D vz/centrality binned histograms
      TH3D* tmpHistTrig;
      leadingArrays = new TObjArray*[binsCentrality];
      for ( int i = 0; i < binsCentrality; ++i ) {
        leadingArrays[i] = new TObjArray();
        for ( int j = 0; j < binsVz; ++j ) {
          // Unique name for each TH3D
          std::stringstream s1, s2;
          s1 << i;
          s2 << j;
          TString leadName = "jet_cent_";
          leadName += s1.str() + "_vz_" + s2.str();
          
          // make each histogram
          tmpHistTrig = new TH3D(leadName, leadName+";eta;phi;centrality", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
          
          // add to the correct bin
          leadingArrays[i]->AddLast( tmpHistTrig );
          
        }
      }
			
			initialized = true;
			return 0;
		}
		else if ( analysisType == "ppdijet" ) {
			hBinVz			= new TH1D( "binvzdist", "Vz Bin Distribution", binsVz, -0.5, (double) binsVz - 0.5 );
			hVz					= new TH1D( "vzdist", "Vz Distribution", 100, -30, 30);
			
			hLeadJetPt 	= new TH1D( "ppleadjetpt", "PP Leading Jet Pt;p_{T}", 80, 0, 80 );
			hLeadEtaPhi = new TH2D( "ppleadjetetaphi", "PP Leading Jet Eta Phi;eta;phi", 40, -1, 1, 40, -pi, pi );
			hSubJetPt 	= new TH1D( "ppsubjetpt", "PP Subleading Jet Pt;p_{T}", 80, 0, 80 );
			hSubEtaPhi 	= new TH2D( "ppsubjetetaphi", "PP Subleading Jet Eta Phi;eta;phi", 40, -1, 1, 40, -pi, pi );
      
      hAssocPt 		= new TH1D("assocpt", "Associated Track Pt;p_{T}", 80, 0, 12 );
      hAssocEtaPhi= new TH2D("assocetaphi", "Associated Track Eta Phi;#eta;#phi", 40, -1, 1, 40, -pi, pi );
			
			hAjHigh 		= new TH1D( "ppajhigh", "PP A_{J} High P_{T} Constituents;A_{J};fraction", 30, 0, 0.9 );
			hAjLow 			= new TH1D( "ppajlow", "PP A_{J} Low P_{T} Constituents;A_{J};fraction", 30, 0, 0.9 );
			
			h3DimCorrLead		= new TH3D("ppleadjetcorr", "PP Lead Jet - Hadron Correlation;#eta;#phi;p_{T}", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
			h3DimCorrSub		= new TH3D("ppsubjetcorr", "PP Sub Jet - Hadron Correlation;#eta;#phi;p_{T}", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
      
      // now build the full 3D vz/centrality binned histograms
      TH3D* tmpHistLead, * tmpHistSub;
      leadingArrays = new TObjArray*[binsCentrality];
      subleadingArrays = new TObjArray*[binsCentrality];
      for ( int i = 0; i < binsCentrality; ++i ) {
        leadingArrays[i] = new TObjArray();
        subleadingArrays[i] = new TObjArray;
        for ( int j = 0; j < binsVz; ++j ) {
          // Unique name for each TH3D
          std::stringstream s1, s2;
          s1 << i;
          s2 << j;
          TString leadName = "pp_lead_cent_";
          TString subName = "pp_sub_cent_";
          leadName += s1.str() + "_vz_" + s2.str();
          subName += s1.str() + "_vz_" + s2.str();
          
          // make each histogram
          tmpHistLead = new TH3D(leadName, leadName+";eta;phi;centrality", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
          tmpHistSub = new TH3D(subName, subName+";eta;phi;centrality", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
          
          // add to the correct bin
          leadingArrays[i]->AddLast( tmpHistLead );
          subleadingArrays[i]->AddLast( tmpHistSub );
          
        }
      }

			
			initialized = true;
			return 0;
		}
		else if ( analysisType == "ppjet" ) {
			hBinVz			= new TH1D( "binvzdist", "Vz Bin Distribution", binsVz, -0.5, (double) binsVz - 0.5 );
			hVz					 = new TH1D("vzdist", "vzdist", 100, -30, 30);
			
			hLeadJetPt 	= new TH1D( "pptriggerjetpt", "PP Trigger Jet Pt;p_{T}", 80, 0, 80 );
			hLeadEtaPhi = new TH2D( "pptriggerjetetaphi", "PP Trigger Jet Eta Phi;eta;phi", 40, -1, 1, 40, -pi, pi );
      
      hAssocPt 		= new TH1D("assocpt", "Associated Track Pt;p_{T}", 80, 0, 12 );
      hAssocEtaPhi= new TH2D("assocetaphi", "Associated Track Eta Phi;#eta;#phi", 40, -1, 1, 40, -pi, pi );
			
			h3DimCorrLead		= new TH3D("pptriggerjetcorr", "PP Jet - Hadron Correlation;#eta;#phi;p_{T}", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
      
      // now build the full 3D vz/centrality binned histograms
      TH3D* tmpHistTrig;
      leadingArrays = new TObjArray*[binsCentrality];
      for ( int i = 0; i < binsCentrality; ++i ) {
        leadingArrays[i] = new TObjArray();
        for ( int j = 0; j < binsVz; ++j ) {
          // Unique name for each TH3D
          std::stringstream s1, s2;
          s1 << i;
          s2 << j;
          TString leadName = "pp_jet_cent_";
          leadName += s1.str() + "_vz_" + s2.str();
          
          // make each histogram
          tmpHistTrig = new TH3D(leadName, leadName+";eta;phi;centrality", binsEta, dEtaLowEdge, dEtaHighEdge, binsPhi, phiLowEdge, phiHighEdge, binsPt, ptLowEdge, ptHighEdge );
          
          // add to the correct bin
          leadingArrays[i]->AddLast( tmpHistTrig );
          
        }
      }
			
			initialized = true;
			return 0;
		}

		else {
			__ERR("Unrecognized analysis type")
			
			initialized = false;
			return -1;
		}
		
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
		
		if ( h3DimCorrLead )
			h3DimCorrLead->Write();
		if ( h3DimCorrSub )
			h3DimCorrSub->Write();
    
    for ( int i = 0; i < binsCentrality; ++i ) {
      if ( leadingArrays )
        if ( leadingArrays[i] )
          leadingArrays[i]->Write();
      if ( subleadingArrays )
        if ( subleadingArrays[i] )
          subleadingArrays[i]->Write();
    }
	}
	
	
	// --------------------------- Histogram Filling Functions ------------------------------- //
	bool histograms::CountEvent( int centrality, int vzbin ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "dijet" || analysisType == "jet" ) {
			hCentVz->Fill( centrality, vzbin );
			return true;
		}
		
		__ERR("hCentVz not initialized for pp")
		return false;
	}
	
	bool histograms::CountEvent( int vzbin ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "ppdijet" || analysisType == "ppjet" ) {
			hBinVz->Fill(  vzbin );
			return true;
		}
		
		__ERR("hBinVz not initialized for auau")
		return false;
	}
	
	bool histograms::FillGRefMult( int gRefMult ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "dijet" || analysisType == "jet" ) {
			hGRefMult->Fill( gRefMult );
			return true;
		}
		
		__ERR("hGRefMult not initialized for pp")
		return false;
	}
	
	bool histograms::FillVz( double vz ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		hVz->Fill( vz );
		return true;
	}
	
	bool histograms::FillAjHigh( double aj ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "dijet" || analysisType == "ppdijet" ) {
			hAjHigh->Fill( aj );
			return true;
		}
		
		__ERR("hAjHigh not initialized for jet-hadron")
		return false;
	}
	
	bool histograms::FillAjLow( double aj ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "dijet" || analysisType == "ppdijet" ) {
			hAjLow->Fill( aj );
			return true;
		}
		
		__ERR("hAjLow not initialized for jet-hadron")
		return false;
	}
	
	bool histograms::FillJetPt( double pt ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "jet" || analysisType == "ppjet" ) {
			hLeadJetPt->Fill( pt );
			return true;
		}
		
		__ERR(" FillJetPt not used for dijet-hadron")
		return false;
	}
	
	bool histograms::FillJetEtaPhi( double eta, double phi ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "jet" || analysisType == "ppjet" ) {
			hLeadEtaPhi->Fill( eta, phi );
			return true;
		}

		__ERR("FillJetEtaPhi() not used for dijet-hadron")
		return false;
	}
	
	bool histograms::FillCorrelation( double dEta,  double dPhi, double assocPt, double weight, int vzBin, int centBin ) {
		
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( dPhi < phiLowEdge )
			dPhi += 2.0*pi;
		
		if ( analysisType == "jet" ) {
			h3DimCorrLead->Fill( dEta, dPhi, assocPt, weight );
      
      // now do the bin-divided fill
      TH3D* tmpHist = (TH3D*) leadingArrays[centBin]->At(vzBin);
      tmpHist->Fill( dEta, dPhi, assocPt, weight );
			return true;
		}
    else if ( analysisType == "ppjet" ) {
      h3DimCorrLead->Fill( dEta, dPhi, assocPt, weight );
      
      // now do the bin-divided fill
      TH3D* tmpHist = (TH3D*) leadingArrays[0]->At(vzBin);
      tmpHist->Fill( dEta, dPhi, assocPt, weight );
      return true;
    }

    
		__ERR("FillCorrelation() not used for dijet-hadron")
		return false;
	}
	
	bool histograms::FillLeadJetPt( double pt ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "dijet" || analysisType == "ppdijet" ) {
			hLeadJetPt->Fill( pt );
			return true;
		}
		
		__ERR("FillLeadJetPt() not used for jet-hadron")
		return false;
	}
	
	bool histograms::FillSubJetPt( double pt ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "dijet" || analysisType == "ppdijet" ) {
			hSubJetPt->Fill( pt );
			return true;
		}
		
		__ERR("FillSubJetPt() not used for jet-hadron")
		return false;
	}
	
	bool histograms::FillLeadEtaPhi( double eta, double phi ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "dijet" || analysisType == "ppdijet" ) {
			hLeadEtaPhi->Fill( eta, phi );
			return true;
		}
		
		__ERR("FillLeadJetEtaPhi() not used for jet-hadron")
		return false;
	}
	
	bool histograms::FillSubEtaPhi( double eta, double phi ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( analysisType == "dijet" || analysisType == "ppdijet" ) {
			hSubEtaPhi->Fill( eta, phi );
			return true;
		}
		
		__ERR("FillSubJetEtaPhi() not used for jet-hadron")
		return false;
	}
	
	bool histograms::FillCorrelationLead( double dEta, double dPhi, double assocPt, double weight, int vzBin, int centBin ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( dPhi < phiLowEdge )
			dPhi += 2.0*pi;
		
		if ( analysisType == "dijet" ) {
			h3DimCorrLead->Fill( dEta, dPhi, assocPt, weight );
      
      // now do the bin-divided fill
      TH3D* tmpHist = (TH3D*) leadingArrays[centBin]->At(vzBin);
      tmpHist->Fill( dEta, dPhi, assocPt, weight );
      
			return true;
		}
    else if ( analysisType == "ppdijet" )  {
      h3DimCorrLead->Fill( dEta, dPhi, assocPt, weight );
      
      // now do the bin-divided fill
      TH3D* tmpHist = (TH3D*) leadingArrays[0]->At(vzBin);
      tmpHist->Fill( dEta, dPhi, assocPt, weight );
      
      return true;

      
    }
    
		
		__ERR("FillCorrelationLead() not used for jet-hadron")
		return false;

	}
	
	bool histograms::FillCorrelationSub( double dEta, double dPhi, double assocPt, double weight, int vzBin, int centBin ) {
		if ( !initialized ) {
			__ERR("histograms instance not initialized")
			return false;
		}
		
		if ( dPhi < phiLowEdge )
			dPhi += 2.0*pi;
		
		if ( analysisType == "dijet" ) {
			h3DimCorrSub->Fill( dEta, dPhi, assocPt, weight );
      
      // now do the bin-divided fill
      TH3D* tmpHist = (TH3D*) subleadingArrays[centBin]->At(vzBin);
      tmpHist->Fill( dEta, dPhi, assocPt, weight );
      
			return true;
		}
    else if ( analysisType == "ppdijet" ) {
      h3DimCorrSub->Fill( dEta, dPhi, assocPt, weight );
      
      // now do the bin-divided fill
      TH3D* tmpHist = (TH3D*) subleadingArrays[0]->At(vzBin);
      tmpHist->Fill( dEta, dPhi, assocPt, weight );
      
      return true;
    }

    
		__ERR("FillCorrelationSub() not used for jet-hadron")
		return false;
	}
	
  bool histograms::FillAssocPt( double pt ) {
    if ( !initialized ) {
      __ERR("histograms instance not initialized")
      return false;
    }
    
    hAssocPt->Fill( pt );
    return true;
  }
  
  bool histograms::FillAssocEtaPhi( double eta, double phi ) {
    if ( !initialized ) {
      __ERR("histograms instance not initialized")
      return false;
    }
    
    hAssocEtaPhi->Fill( eta, phi );
    return true;
  }
  
}

