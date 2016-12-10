// histogram class used in my correlation analysis
// holds correlations in bins of Aj, Vz and Centrality
// Nick Elsey

#include "corrParameters.hh"
#include "corrFunctions.hh"

#ifndef HISTOGRAMS_HH
#define HISTOGRAMS_HH

namespace corrAnalysis {
 
  // Histogram holder
  // ----------------
  class histograms {
    
  private:
    
    std::string analysisType;			// Used by Init() to create proper histograms
    bool initialized;							// Used for control flow - must be true before filling
    
    unsigned binsEta;             // Binning for histograms
    unsigned binsPhi;             // Binning for histograms
    
    float phiBinShift;            // value for shifting histogram edges to force a bin centered at zero
    float etaBinShift;            // value for shifting histogram edges to force a bin centered at zero
    
    // Event statistics
    TH3D* hCentVz;
    TH2D* hBinVz;
    TH1D* hGRefMult;
    TH1D* hVz;
    
    // Jet kinematics
    TH1D* hLeadJetPt;
    TH2D* hLeadEtaPhi;
    TH1D* hSubJetPt;
    TH2D* hSubEtaPhi;
    
    // Associated kinematics
    TH1D* hAssocPt;
    TH2D* hAssocEtaPhi;
    
    // A_j if doing dijets
    TH1D* hAjHigh;
    TH1D* hAjLow;
    TH3D* hAjDif;
    
    // Correlations that are not differentiated by Vz and Centrality
    TH3D* h3DimCorrLead;
    TH3D* h3DimCorrSub;
    
    // Holders for the vz/cent binned histograms
    TObjArray*** leadingArrays;
    TObjArray*** subleadingArrays;
    
    // a histogram to see if there is some event or jet structure
    // relating to Aj
    TH3D* hAjStruct;
    
    // Used internally when filling histograms
    bool IsPP();
    bool IsAuAu();
    bool IsDijet();
    bool IsJet();
    bool IsMix();
    
    // Checked to make sure the histogram class
    // has been initialized before filling
    bool IsInitialized();
    
    // Used internally during initialization
    // to generate the histogram arrays properly
    // depending on analysis settings
    void BuildArrays();
    
    // Used internally to pick histogram edges
    // that give a bin centered at zero in correlation plots
    void FindBinShift();

    // Used to find the respective Aj bin
    int FindAjBin(double aj);
    
  public:
    histograms( );
    histograms( std::string type, unsigned binsEta = 24, unsigned binsPhi = 24 ); // In general, this should be used, passing "dijet" or "jet" for analysis
    ~histograms();
    
    // Deletes all histograms
    void Clear();
    
    // Get the analysis type
    std::string GetAnalysisType()  { return analysisType; }
    
    // Can set analysisType or Aj splitting - careful, if it changes after Init() is called
    // Reinitialization will be needed
    bool SetAnalysisType( std::string type );
    
    // Must be called before filling
    // Checks analysisType and creates histograms
    int	Init();
    
    // Writes histograms to current root directory
    void Write();
    
    // Get Histograms
    TH3D* GetCentVz()				{ return hCentVz; }
    TH2D* GetBinVz()				{ return hBinVz; }
    TH1D*	GetGRefMult()			{ return hGRefMult; }
    TH1D* GetVz()						{ return hVz; }
    TH1D* GetLeadPt()				{ return hLeadJetPt; }
    TH2D* GetLeadEtaPhi() 	{ return hLeadEtaPhi; }
    TH1D* GetSubPt() 				{ return hSubJetPt; }
    TH2D* GetSubEtaPhi()		{ return hSubEtaPhi; }
    TH1D* GetAjHigh()				{ return hAjHigh; }
    TH1D* GetAjLow()				{ return hAjLow; }
    TH3D* Get3DLeadCorr()		{ return h3DimCorrLead; }
    TH3D* Get3DSubCorr()		{ return h3DimCorrSub; }
    
    
    // Fill histogram functions
    bool CountEvent( int vzBin, int centrality = 0, double aj = 0.01 ); 	// Used to count events
    
    bool FillGRefMult( int gRefMult );							// For AuAu events, records gRefMult
    bool FillVz( double vz );												// Records Vz distribution
    
    bool FillAjHigh( double aj );										// Records Aj for initial hard jets
    bool FillAjLow( double aj );										// Records Aj for jets with soft constituents
    bool FillAjDif( double ajHigh, double ajLow );         // Records Aj high, Ajlow, Ajdif
    
    // For jet-hadron
    bool FillJetPt( double pt );										// For Jet-hadron: records accepted trigger jet pt
    bool FillJetEtaPhi( double eta, double phi );		// Records accepted trigger jet eta-phi
    // records trigger-associated correlations
    bool FillCorrelation( double dEta, double dPhi, double assocPt, double weight, int vzBin, int centBin = 0 );
    
    // For dijet-hadron
    bool FillLeadJetPt( double pt );								// For dijet-hadron: records lead jet pt
    bool FillSubJetPt( double pt );									// Records sub jet pt
    bool FillLeadEtaPhi( double eta, double phi );	// Records lead jet eta-phi
    bool FillSubEtaPhi( double eta, double phi );		// Records sub jet eta-phi
    // Records trigger-associated correlations with trigger = leading/subleading
    bool FillCorrelationLead( double dEta, double dPhi, double assocPt, double weight, double aj, int vzBin, int centBin = 0 );
    bool FillCorrelationSub( double dEta, double dPhi, double assocPt, double weight, double aj, int vzBin, int centBin = 0 );
    
    // Associated track info
    bool FillAssocPt( double pt );
    bool FillAssocEtaPhi( double eta, double phi );
    
    // Looking for correlations between Aj, Cent, and Pt
    bool FillAjStruct( double aj, double pt, int centrality );
    
  };

  
}

#endif
