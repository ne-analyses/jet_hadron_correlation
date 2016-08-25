// Making averaged efficiencies for each centrality
// AuAu - STAR - year 7
// Nick Elsey

// ROOT headers
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
#include "TStopWatch.h"

// standard library headers
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>

#include "ktTrackEff.hh"


int main () {
	TString ktEff_hist_name = "src/run7eff.root";
	
	TFile ktEff_hist_file( ktEff_hist_name, "READ" );
	
	TH2D* h2D_eta_phi[3];
	TH1D* h1D_eta_corr[3];
	
	ktTrackEff efficiency_engine( ktEff_hist_name );
	// Get the centrality bin histograms
	// And makes the output histograms
	// Plus one for pp y7
	TH2D* out_hist[3];
	int nBinsPt = 30;
	double ptLow = 0.0;
	double ptHigh = 5.0;
	int nBinsEta = 30;
	double etaLow = -1;
	double etaHigh = 1;
	TH2D* out_hist_pp = new TH2D("pp_efficiency_pt_eta", "pp_efficiency_pt_eta;pt;eta;efficiency", nBinsPt, ptLow, ptHigh, nBinsEta, etaLow, etaHigh );
	// And the final averaged result
	TH1D* out_pt[3];
	TH1D* out_pt_pp;
	for (int i = 0; i < 3; ++i ) {
		// get the histogram name
		TString h2D_name_base = "ptetaScale_";
		TString h1D_name_base = "etaScale_";
		std::stringstream ss;
		ss << i;
		h2D_name_base += ss.str();
		h1D_name_base += ss.str();
		
		h2D_eta_phi[i] = (TH2D*) ktEff_hist_file.Get(h2D_name_base);
		h1D_eta_corr[i] = (TH1D*) ktEff_hist_file.Get(h1D_name_base);
		
		// make a name for the output
		TString out_name = "efficiency_pt_eta_cent_";
		out_name += ss.str();
		out_hist[i] = new TH2D(out_name, out_name+";pt;eta;efficiency", nBinsPt, ptLow, ptHigh, nBinsEta, etaLow, etaHigh );
	}
	
	// Make an output file
	TFile janky_eff_file("src/janky_eff.root", "RECREATE");
	
	// Now to calculate the corrections... Following the steps in ktTrackEff
	double pt_bin_width = (ptHigh - ptLow)/(double)nBinsPt;
	double eta_bin_width = (etaHigh - etaLow ) / (double)nBinsEta;
	double pt_half_bin = pt_bin_width/2.0;
	double eta_half_bin = eta_bin_width/2.0;
	
	std::cout<<" pt bin width: "<<pt_bin_width;
	std::cout<<" pt half bin: "<<pt_half_bin;
	std::cout<<" eta bin width: "<<eta_bin_width;
	std::cout<<" eta half bin: "<<eta_half_bin;
	for ( int i = 0; i < 3; ++i )
		for ( int j = 0; j < nBinsPt; ++j )
			for ( int k = 0; k < nBinsEta; ++k ) {
				double bin_pt = ptLow + pt_half_bin + j * pt_bin_width;
				double bin_eta = etaLow + eta_half_bin + k * eta_bin_width;
				
				double efficiency = efficiency_engine.EffAAY07( bin_eta, bin_pt, i);
				double efficiency_pp = efficiency_engine.EffPPY06(bin_eta, bin_pt );
				out_hist[i]->Fill( bin_pt, bin_eta, efficiency );
				
				if ( i == 0 )
					out_hist_pp->Fill( bin_pt, bin_eta, efficiency_pp );
				
			}
				
	
	// Now to average over eta... Write to file
	for (int i = 0; i < 3; ++i ) {
		out_hist[i]->Write();
		out_pt[i] = (TH1D*) out_hist[i]->ProjectionX();
		std::stringstream ss;
		ss << i;
		TString h1d_name = "efficiency_pt_cent_";
		out_pt[i]->SetName(h1d_name + ss.str() );
		out_pt[i]->SetTitle(h1d_name + ss.str() + ";pt;efficiency");
		out_pt[i]->Scale(1.0/(double)nBinsEta);
		out_pt[i]->Rebin(3);
		out_pt[i]->Scale(1.0/3.0);
		out_pt[i]->Write();
	}
	
	std::cout<<std::endl;
	
	// Now do the same for pp...
	out_hist_pp->Write();
	out_pt_pp = (TH1D*) out_hist_pp->ProjectionX();
	out_pt_pp->SetName( "pp_efficiency_pt" );
	out_pt_pp->SetTitle( "pp_efficiency_pt;pt;efficiency" );
	out_pt_pp->Scale(1.0/(double)nBinsEta);
	out_pt_pp->Rebin(3);
	out_pt_pp->Scale(1.0/3.0);
	out_pt_pp->Write();
	
	
	return 0;
}