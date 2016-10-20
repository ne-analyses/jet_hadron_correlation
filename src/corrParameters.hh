// Parameters for all correlation analyses
// Nick Elsey

#include <string>

// Define a namespace for the variables

#ifndef CORRPARAMETERS_HH
#define CORRPARAMETERS_HH

#define __ERR(message) {std::cerr << "[" << __FILE__ << "::" << __func__ << "()] -- ERR: " << message << std::endl;}
#define __OUT(message) {std::cout << "[" << __FILE__ << "::" << __func__ << "()] -- OUT: " << message << std::endl;}

namespace corrAnalysis {
	
	// Math
	// ----------------------------
	// Define pi to float precision
	const float pi = 3.141592;
	
	// Reader settings
	// ----------------------------
  
  // default nEvents for all events in chain
  const int		allEvents = -1;									// selects all events in the chain
  
	// TriggerStrings
	const std::string triggerAll = "All";				// accept all events
	const std::string	triggerHT = "HT";					// accept events triggered by high E hit in calorimeter
	const std::string triggerMB = "MB";					// accept minimum bias events
	const std::string triggerPP = "pp";					// accept pp events
	const std::string triggerPPHT = "ppHT";			// accept pp high tower events
  const std::string triggerPPJP = "ppJP";			// accept pp jet patch events
  
  const double triggerThreshold = 5.0;				// the required energy for a tower to be considered a trigger
	
	// Event
  const int 		y7RefMultCut = 269;										// refmult cut for 0-20% centrality
	const bool 		hadronicCorrection = true;						// will towers have matched tracks energy subtracted?
	const double	hadronicCorrectionFraction = 0.9999;	// fraction of energy subtracted
	const double  vertexZCut = 30;											// | Vz | < vertexZCut - vertex Z position requirement
	const double	vertexZDiffCut = 9999;								// | Vz - VPDVz | < vertexZDiffCut
  
  // Max Et and Pt
  const double  eventPtCut = 30.0;					// if a track has higher than 30 GeV Pt, reject event
	const double  eventEtCut = 30.0;					// if a tower has higher than 30 GeV E, reject event
  
	// TPC tracks
	const double 	DCACut = 1.0;								// distance of closest approach to primary vertex
	const int 		minFitPoints = 20;					// required number of fit points used in track reconstruction
	const double 	minFitFrac = 0.52;					// required fraction of fit points out of the total available
	const double  trackPtCut = 9999.0;				// we want this to be nonexistent so that the event is rejected by eventPtCut
	
	// Calorimeter towers
	const double  towerEtCut = 9999.0;				// we want this to be nonexistent so that the event is rejected by eventEtCut
	// Trying to make everything machine independent
	// Define a path to bad tower list
	const std::string y7AuAuTowerList = "src/y7_AuAu_HT_hot_list.txt";
  const std::string y6PPTowerList = "src/Combined_y7_PP_Nick.txt";
	
	// If that doesnt work, here are machine dependent paths
	// bad tower list location on RHIC121
	const std::string y7AuAuTowerRHIC121 = "/Users/nickelsey/physics/software/eventStructuredAu/y7_AuAu_HT_hot_list.txt";
  const std::string y7PPTowerRHIC121 = "/Users/nickelsey/physics/software/eventStructuredAu/Combined_y7_PP_Nick.txt";
	// bad tower list location on WSU grid
	const std::string y7AuAuTowerWSUGRID = "NA";
  const std::string y7PPTowerWSUGRID = "NA";
	// bad tower list location on Gauss
	const std::string y7AuAuTowerGAUSS   = "/Users/nick/physics/software/eventStructuredAu/y7_AuAu_HT_hot_list.txt";
  const std::string y7PPTowerGAUSS   = "/Users/nick/physics/software/eventStructuredAu/Combined_y7_PP_Nick.txt";
	// --------------------------------------------------------
	
	// Binning for histograms,
	// and centrality/Vz bins
	// ---------------------------
	// Vz binning
	const double 	vzRange = 60.0; 					// total accepted vz range (symmetric about vz = 0 )
	const int		 	binsVz = 40;							// range is split into vzBins number of bins
	const double 	dVz = vzRange/binsVz;			// each bin has a range of dVz
	const double 	vzLowEdge = -vzRange/2.0;	// lower edge for accepted Vz range
	const double 	vzHighEdge = vzRange/2.0;	// upper edge for accepted Vz range
	
	// Centrality binning
	const int 		binsCentrality = 9;				// using the 9 centrality bins definition
	// The centrality bins
	const std::string centBin9String[9] = { "0-5", "5-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80" };
	const double 	centLowEdge = -0.5;				// lower edge for centrality ( to keep 0 in first bin, 8 in 9th bin )
	const double 	centHighEdge = (double) binsCentrality - 0.5; // upper edge for centrality
	// GRefMult centrality definitions for y7
	const int 		y7RefMultCent[9] = { 10, 21, 39, 69, 114, 178, 269, 399, 485 };
	
	// Pt binning
	const int 	 	binsPt = 120;							// number of pt bins
	const double 	ptLowEdge = 0.0;					// lower edge for accepted pt range
	const double 	ptHighEdge = 30.0;				// upper edge for accepted pt range
	
	// Phi binning ( full 2*pi coverage, from -pi/4 to 7*pi/4 )
	const int 	 	binsPhi = 25;							// number of phi bins
	const double 	phiLowEdge = -pi/2.0;			// lower edge for accepted phi range
	const double 	phiHighEdge = 3.0*pi/2.0;	// upper edge for accepted phi range
	
	// Eta binning
	const int		 binsEta = 25;							// number of eta bins
	const double etaLowEdge = -1;						// lower edge for accepted eta range
	const double etaHighEdge = 1;						// upper edge for accepted phi range
	const double dEtaLowEdge = 2*etaLowEdge;		// lower edge for accepted dEta range
	const double dEtaHighEdge = 2*etaHighEdge; 	// upper edge for accepted dEta range
	
	// Jetfinding Settings
	// ---------------------------
	// (Jet Pt min/max are set in the analysis)
	const double jetPtHardMax = 1000.0;   	// for consistency, used when no max is needed
	
	const double jetRadius = 0.4;
	
	// Constituent parameters used in jet selection
	const double maxTrackRap 	= 1.0;				// accept tracks with rapidity [ -1, 1 ]
	const double trackMinPt  	= 0.2;				// only accept tracks with Pt > 0.2
	const double hardTrackMinPt = 2.0;			// hard jets clustered with pt > 2.0
  
  // dijet analysis requires jets to be back to back
  const double jetDPhiCut 	= 0.4;				// require jets to be pi - 0.4 away in phi
	
	// Variables for the area definition
	// Ghosts, etc
	const int ghostRepeat = 1;
	const double ghostArea = 0.01;
	
	
	// Associated efficiency information
	// ----------------------------
	// Trying to make this all machine independent
	// Path to efficiency file from jet_hadron_corr/
	const std::string y7EfficiencyFile = "src/run7eff.root";
	
	// If that doesnt work, machine dependent stuff
	// Efficiency file location on RHIC121
	const std::string y7EfficiencyRHIC121 = "/Users/nickelsey/physics/analysis/jet_hadron_corr/src/run7eff.root";
	// Efficiency file location on WSU grid
	const std::string y7EfficiencyWSUGRID = "NA";
	// Efficiency file location on Gauss
	const std::string y7EfficiencyGAUSS   = "/Users/nick/physics/analysis/jet_hadron_corr/src/run7eff.root";
	// --------------------------------------------------------
	
	// We only have efficiency for 0-20%
	// For year 7
	// These can't be changed
	const int y7EfficiencyRefCentLower = 6;
	const int y7EfficiencyRefCentUpper = 8;
	
}

#endif
