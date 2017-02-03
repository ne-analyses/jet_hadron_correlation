// Implementation of basic functions used
// in the dijet-hadron and jet-hadron
// correlation analysis
// Nick Elsey

#include "corrFunctions.hh"
#include "corrParameters.hh"
#include "histograms.hh"

#include <time.h>
#include <random>

namespace jetHadron {
	
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

  // Other string manipulations
  // Checks if a string begins with a certain substring
  bool BeginsWith (std::string const &full_string, std::string const &beginning) {
    if (full_string.length() >= beginning.length()) {
      return (0 == full_string.compare (0, beginning.length(), beginning) );
    } else {
      return false;
    }
  }
  
  // Used to pull the current directory from its absolute path
  // Used in mixing to get mixing parameters
  std::string GetDirFromPath( std::string path ) {
    return path.substr( path.find_last_of("\\/")+1, path.size() );
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
  void ConvertTStarJetVector( TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles, bool ClearVector, double towerScale ) {
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
      if ( sv->GetCharge() == 0 )
        tmpPJ *= towerScale;
      tmpPJ.set_user_index( sv->GetCharge() );
      particles.push_back( tmpPJ );
      
      
    }
  }
  
  // applies an effective 90% relative efficiency compared to auau
  void ConvertTStarJetVectorPP( TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles, ktTrackEff& eff, int64_t seed, bool ClearVector, double towerScale ) {
    // Empty the container
    // if called for
    if ( ClearVector )
      particles.clear();
    
    // create a RNG for shuffling events
    std::random_device rd;
    std::mt19937 g(rd());
    g.seed( seed );
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    // Transform TStarJetVectors into (FastJet) PseudoJets
    // ---------------------------------------------------
    TStarJetVector* sv;
    for ( int i=0; i < container->GetEntries() ; ++i ){
      sv = container->Get(i);
      
      if ( sv->GetCharge() ) {
        double ratio = eff.EffRatio_20(sv->Eta(),sv->Pt());
        double random_ = dis(g);
        if ( random_ > ratio ) {
          continue;
        }
      }
      fastjet::PseudoJet tmpPJ = fastjet::PseudoJet( *sv );
      if ( sv->GetCharge() == 0 )
        tmpPJ *= towerScale;
      tmpPJ.set_user_index( sv->GetCharge() );
      particles.push_back( tmpPJ );
    }
    
  }
  
  // For AuAu being embedded into PP
  void ConvertTStarJetVectorPPEmbedded( TStarJetVectorContainer<TStarJetVector>* container, std::vector<fastjet::PseudoJet> & particles, bool allTracks, double towerScale ) {
    
    // Transform TStarJetVectors into (FastJet) PseudoJets
    // ---------------------------------------------------
    TStarJetVector* sv;
    for ( int i=0; i < container->GetEntries() ; ++i ){
      sv = container->Get(i);
      
      fastjet::PseudoJet tmpPJ = fastjet::PseudoJet( *sv );
      if ( sv->GetCharge() == 0 )
        tmpPJ *= towerScale;
      tmpPJ.set_user_index( sv->GetCharge() );
      double pt = tmpPJ.pt();
      
      // only add if allTracks is selected, or if the track has pt > 2.0
      if ( allTracks || pt > 2.0 )
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
        if ( trigger->GetTriggerFlag() == 1 || trigger->GetTriggerFlag() == 118 || trigger->GetTriggerFlag() == 125 ) {
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
    
    // if we're using triggers, run over all towers and get any with Et > triggerThreshold
    if ( requireTrigger ) {
      for ( int i = 0; i < ppParticles.size(); ++i ) {
        fastjet::PseudoJet tmpParticle = ppParticles[i];
        if ( tmpParticle.user_index() == 0 && tmpParticle.pt() > triggerThreshold )
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
  void InitReader( TStarJetPicoReader & reader, TChain* chain, std::string collisionType, std::string triggerString, double softwareTrigger, int nEvents ) {
    
    // First tolower() on the analysisType
    // shouldnt be necessary....
    std::transform(collisionType.begin(), collisionType.end(), collisionType.begin(), ::tolower);
    
    // set the chain
    reader.SetInputChain( chain );
    // apply hadronic correction - subtract 100% of charged track energy from towers
    reader.SetApplyFractionHadronicCorrection( hadronicCorrection );
    reader.SetFractionHadronicCorrection( hadronicCorrectionFraction);
    reader.SetRejectTowerElectrons( kFALSE );
    
    // Event and track selection
    // -------------------------
    
    TStarJetPicoEventCuts* evCuts = reader.GetEventCuts();
    evCuts->SetTriggerSelection( triggerString.c_str() );
    evCuts->SetVertexZCut ( vertexZCut );
    evCuts->SetMaxEventPtCut( eventPtCut );
    evCuts->SetMaxEventEtCut( eventEtCut );
    
    if( fabs(softwareTrigger) < 0.0001 ) {
      evCuts->SetMinEventEtCut( softwareTrigger );
      evCuts->SetUseRawForMinEventEtCut( true );
    }
    
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
    if ( collisionType == "auau" || collisionType == "pp" ) {
    	towerCuts->AddBadTowers( y7AuAuTowerList.c_str() );
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
      
      // SWITCH TO HARD CORE
      //matchedToDijet.push_back( hardJets.at(0) );
      //matchedToDijet.push_back( hardJets.at(1) );
      
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
        
        // check to make sure the matched jets are within the
        // accepted eta range
        for ( int i = 0; i < matchedToDijet.size(); ++i )
          if ( fabs( matchedToDijet[i].eta() ) > jetHadron::maxTrackRap - jetRadius )
            return std::vector<fastjet::PseudoJet>();
        
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
    if ( fabs( assocTrack.eta() ) > maxTrackRap )			{ return false; }
    if ( assocTrack.user_index() == 0 )  			{ return false; }
    
    // Check to make sure the efficiency is not crazy
    // ( the parameterization isnt perfect, about 3.5% of tracks return nonsense efficiencies )
    if ( efficiency <= 0.01 )        { return false;  }
    if ( efficiency > 1.0 ) 				{ return false;  }
    
    return true;
  }
  
  bool correlateLeading( std::string analysisType, int vzBin, int centBin, histograms* histogram, fastjet::PseudoJet& leadJet, fastjet::PseudoJet& assocTrack, double efficiency, double aj ) {
    
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
    histogram->FillCorrelationLead( deltaEta, deltaPhi, assocPt, weight, aj, vzBin, centBin );
    
    return true;
  }
  
  bool correlateSubleading( std::string analysisType, int vzBin, int centBin, histograms* histogram, fastjet::PseudoJet& subJet, fastjet::PseudoJet& assocTrack, double efficiency, double aj ) {
    
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
    histogram->FillCorrelationSub( deltaEta, deltaPhi, assocPt, weight, aj, vzBin, centBin );
    
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
  
  // --------------------------
  // ------ Event Mixing ------
  // --------------------------
  
  // This function will check if it can parse the analysis String
  // if its unrecognized it will use defaults
  // if it recognizes the string but can't parse, it returns -2, which will exit the mixing
  int GetVarsFromString( std::string& analysisType, std::string analysisString, double& leadPt, double& subPt, double& maxPt, double& jetRadius, double& hardPt, bool& useEff, bool& reqTrigger, unsigned& binsEta, unsigned& binsPhi ) {
    
    // First, pick out the analysis type from the analysis string
    if ( BeginsWith( analysisString, "dijet" ) )
      analysisType = "dijetmix";
    else if ( BeginsWith( analysisString, "ppdijet" ) )
      analysisType = "ppdijetmix";
    else if ( BeginsWith( analysisString, "jet" ) )
      analysisType = "jetmix";
    else if ( BeginsWith( analysisString, "ppjet" ) )
      analysisType = "ppjetmix";
    
    // everything matches: we now split the string
    std::stringstream parser(analysisString);
    std::string tmp;
    std::vector<std::string> varHolder;
    while ( std::getline( parser, tmp, '_' ) ) {
      varHolder.push_back(tmp);
    }
    
    // convert and set the analysis variables
    for (int i = 0; i < varHolder.size()-1; ++i ) {
      if ( varHolder[i] == "lead" )
        leadPt = atof( varHolder[i+1].c_str() );
      else if ( varHolder[i] == "sub" )
        subPt = atof( varHolder[i+1].c_str() );
      else if ( varHolder[i] == "max" )
        maxPt = atof( varHolder[i+1].c_str() );
      else if ( varHolder[i] == "rad" )
        jetRadius = atof( varHolder[i+1].c_str() );
      else if ( varHolder[i] == "hardpt" )
        hardPt = atof( varHolder[i+1].c_str() );
      else if ( varHolder[i] == "eta" )
        binsEta = atoi( varHolder[i+1].c_str() );
      else if ( varHolder[i] == "phi" )
        binsPhi = atoi( varHolder[i+1].c_str() );
      else if ( varHolder[i] == "trigger" ) {
        if ( varHolder[i+1] == "true" )
          reqTrigger = true;
        else
          reqTrigger = false;
      }
      else if ( varHolder[i] == "eff" ) {
        if ( varHolder[i+1] == "true" )
          useEff = true;
        else
          useEff = false;
      }
    }
    
    if ( leadPt == -999 || subPt == -999 || jetRadius == -999 || maxPt == -999 || binsEta == 1000 || binsPhi == 1000 )
      return -2;
    
    return 1;
  }
  
  // Used to decide, if using HT events
  // What the max jet pt can be before the event must be thrown away
  double GetMixEventJetPtMax( bool useMB, std::string analysisType, double leadJetPtMin ) {
    // Check to make sure the analysis type is 'mix'
    // If not, return -999
    if ( !HasEnding( analysisType, "mix" ) ) {
      __ERR("This function is only used in event mixing")
      return -999;
    }
    
    
    // If using MB data, return 0, unneccessary
    // To jetfind, use all events
    if ( useMB )
      return 0;
    // If the jet pt < 10 GeV ( getting close to trigger threshhold )
    // We will not be able to use HT data
    else if ( leadJetPtMin < 10.0 )
      return -1;
    // Otherwise return 80% of the leading jet pt min
    // ( Or trigger jet for jet-hadron )
    else
      return 0.8*leadJetPtMin;
    
  }

  // Decides whether an event should be used
  // In mixing or not - logic depends on analysis type
  // And on Data set being used
  bool UseEventInMixing( std::string analysisType, bool isMB, std::vector<fastjet::PseudoJet>& highPtConsJets, int refMult, int vzBin ) {
    
    // If it is HT data and a jet was found above the threshold, discard the event
    if ( !isMB && highPtConsJets.size() > 0 )
      return false;
   
    //If it is AuAu data, check reference centrality boundaries
    if ( analysisType == "dijetmix" || analysisType == "jetmix" ) {
      int refCentrality = GetReferenceCentrality( refMult );
      // Check to see if we use those centralities
      if ( refCentrality < 0 )                      								 	{ return false; }
      if ( refCentrality < jetHadron::y7EfficiencyRefCentLower )   { return false; }
      if ( refCentrality > jetHadron::y7EfficiencyRefCentUpper )   { return false; }
    }
    
    // Make sure Vz bin is accepted
    if ( vzBin < 0 )
      return false;
    
    // it made the cuts - we'll use the event
    return true;
  }
  
  
} // end namespace



