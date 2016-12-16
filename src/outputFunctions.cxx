// Implementation of basic functions
// used in the output workflow of
// the jet hadron correlation analysis

#include "outputFunctions.hh"
#include "corrParameters.hh"
#include "histograms.hh"

// to build directories we use boost
#include "boost/filesystem.hpp"

// the grid does not have std::to_string() for some ungodly reason
// replacing it here. Simply ostringstream
namespace patch {
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}

namespace jetHadron {
  
  // binSelector functionality to set bin information to match
  // current correlations
  void binSelector::SetHistogramBins( TH2F* h ) {
    
    bindEta = h->GetXaxis()->GetNbins();
    dEtaLow = h->GetXaxis()->GetBinLowEdge(1);
    dEtaHigh = h->GetXaxis()->GetBinUpEdge( h->GetXaxis()->GetNbins() );
    bindPhi = h->GetYaxis()->GetNbins();
    dPhiLow = h->GetYaxis()->GetBinLowEdge(1);
    dPhiHigh = h->GetYaxis()->GetBinUpEdge( h->GetYaxis()->GetNbins() );
    
  }
  
  // setting ranges for the histogram
  // if the radius is NOT 0.4, this needs to be used...
  void binSelector::ChangeRadius( double R ) {
    dEtaAcceptanceLow = R - 2.0;
    dEtaAcceptanceHigh = 2.0 - R;
  }
  
  
  // Function used to read in histograms from
  // the files passed in - it returns the correlations,
  // and the number of events, and selects using the centralities,
  // vz bin range, and aj ranges passed in via binSelector
  int ReadInFiles(std::vector<TFile*>& filesIn, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& leadingCorrelations, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& subLeadingCorrelations, std::vector<TH3F*>& nEvents, binSelector selector ) {
    
    // loop over all files and all aj, centrality, and vz bins
    // to return a 4D vector of histograms
    for ( int i = 0; i < filesIn.size(); ++i ) {
      // tell the user what is going on
      std::string outMsg = "Reading in file " + patch::to_string(i);
      __OUT(outMsg.c_str() )
      
      // for each file, get the number of events
      nEvents.push_back( (TH3F*) filesIn[i]->Get("nevents") );
      // rename because root can't handle simple crap
      std::string tmpName = "corr_nevents_" + patch::to_string(i);
      nEvents[i]->SetName( tmpName.c_str() );
      
      // push back the vectors
      leadingCorrelations.push_back( std::vector<std::vector<std::vector<TH3F*> > >() );
      subLeadingCorrelations.push_back( std::vector<std::vector<std::vector<TH3F*> > >() );
      
      for ( int j = selector.centLow; j <= selector.centHigh; ++j ) {
        
        int cent_index = j - selector.centLow;
        
        // push back the vectors
        leadingCorrelations[i].push_back( std::vector<std::vector<TH3F*> >() );
        subLeadingCorrelations[i].push_back( std::vector<std::vector<TH3F*> >() );
        
        for ( int k = selector.vzLow; k <= selector.vzHigh; ++k ) {
          
          int vz_index = k - selector.vzLow;
          
          // push back the vectors
          leadingCorrelations[i][cent_index].push_back( std::vector<TH3F*>() );
          subLeadingCorrelations[i][cent_index].push_back( std::vector<TH3F*>() );
          
          for ( int l = selector.ajLow; l <= selector.ajHigh; ++l ) {
            
            int aj_index = l - selector.ajLow;
            
            // build the in-file histogram names
            std::string leadName = "lead_aj_" + patch::to_string(l) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k);
            std::string subLeadName = "sub_aj_" + patch::to_string(l) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k);
            
            // check to make sure the file is properly formatted
            if ( !filesIn[i]->Get( leadName.c_str() ) )  {
              __ERR("Can't find histograms - maybe it has mixing correlations not signal?")
              return -1;
            }
            
            // get the correlation histograms
            leadingCorrelations[i][cent_index][vz_index].push_back( (TH3F*) filesIn[i]->Get( leadName.c_str() ) );
            
            // check to make sure it was successful
            if ( !leadingCorrelations[i][cent_index][vz_index][aj_index] ) {
              std::string errorMsg = "Couldn't read in leading correlation: " + patch::to_string(i) + " " + patch::to_string(j) + " " + patch::to_string(k) + " " + patch::to_string(l);
              __ERR( errorMsg.c_str() )
              continue;
            }
            
            
            // if subleading correlations are there, load them in
            if ( filesIn[i]->Get( subLeadName.c_str() ) ) {
              subLeadingCorrelations[i][cent_index][vz_index].push_back( (TH3F*) filesIn[i]->Get( subLeadName.c_str() ) );
            }
            else {
              subLeadingCorrelations[i][cent_index][vz_index].push_back(0x0);
            }
            
            //now differentiate the names by file
            leadName = "corr_file_" + patch::to_string(i) + leadName;
            subLeadName = "corr_file_" + patch::to_string(i) + subLeadName;
            
            leadingCorrelations[i][cent_index][vz_index][aj_index]->SetName( leadName.c_str() );
            if ( subLeadingCorrelations[i][cent_index][vz_index][aj_index] ) {
              subLeadingCorrelations[i][cent_index][vz_index][aj_index]->SetName( subLeadName.c_str() );
            }
          }
        }
      }
    }
    return 1;
  }
  
  int ReadInFilesMix(std::vector<TFile*>& filesIn, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& leadingCorrelations, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& subLeadingCorrelations, std::vector<TH3F*>& nEvents, binSelector selector ) {
    
    // loop over all files and all aj, centrality, and vz bins
    // to return a 4D vector of histograms
    for ( int i = 0; i < filesIn.size(); ++i ) {
      // tell the user what is going on
      std::string outMsg = "Reading in file " + patch::to_string(i);
      __OUT(outMsg.c_str() )
      
      // for each file, get the number of events
      nEvents.push_back( (TH3F*) filesIn[i]->Get("nevents") );
      // rename because root can't handle simple crap
      std::string tmpName = "mix_nevents_" + patch::to_string(i);
      nEvents[i]->SetName( tmpName.c_str() );
      
      // push back the vectors
      leadingCorrelations.push_back( std::vector<std::vector<std::vector<TH3F*> > >() );
      subLeadingCorrelations.push_back( std::vector<std::vector<std::vector<TH3F*> > >() );
      
      for ( int j = selector.centLow; j <= selector.centHigh; ++j ) {
        
        int cent_index = j - selector.centLow;
        
        // push back the vectors
        leadingCorrelations[i].push_back( std::vector<std::vector<TH3F*> >() );
        subLeadingCorrelations[i].push_back( std::vector<std::vector<TH3F*> >() );
        
        for ( int k = selector.vzLow; k <= selector.vzHigh; ++k ) {
          
          int vz_index = k - selector.vzLow;
          
          // push back the vectors
          leadingCorrelations[i][cent_index].push_back( std::vector<TH3F*>() );
          subLeadingCorrelations[i][cent_index].push_back( std::vector<TH3F*>() );
          
          for ( int l = selector.ajLow; l <= selector.ajHigh; ++l ) {
            
            int aj_index = l - selector.ajLow;
            
            // build the in-file histogram names
            std::string leadName = "mix_lead_aj_" + patch::to_string(l) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k);
            std::string subLeadName = "mix_sub_aj_" + patch::to_string(l) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k);
            
            // check to make sure the file is properly formatted
            if ( !filesIn[i]->Get( leadName.c_str() ) ) {
              __ERR("Can't find histograms - maybe it has signal correlations not event mixing?")
              return -1;
            }
            
            // get the correlation histograms
            leadingCorrelations[i][cent_index][vz_index].push_back( (TH3F*) filesIn[i]->Get( leadName.c_str() ) );
            
            // check to make sure it was successful
            if ( !leadingCorrelations[i][cent_index][vz_index][aj_index] ) {
              std::string errorMsg = "Couldn't read in leading correlation: " + patch::to_string(i) + " " + patch::to_string(j) + " " + patch::to_string(k) + " " + patch::to_string(l);
              __ERR( errorMsg.c_str() )
              continue;
            }
            
            // if this is a dijet analysis, get the subleading correlations
            if ( filesIn[i]->Get( subLeadName.c_str() ) ) {
              subLeadingCorrelations[i][cent_index][vz_index].push_back( (TH3F*) filesIn[i]->Get( subLeadName.c_str() ) );
            }
            else {
              subLeadingCorrelations[i][cent_index][vz_index].push_back(0x0);
            }
            
            //now differentiate the names by file
            leadName = "mix_corr_file_" + patch::to_string(i) + leadName;
            subLeadName = "mix_corr_file_" + patch::to_string(i) + subLeadName;
            
            leadingCorrelations[i][cent_index][vz_index][aj_index]->SetName( leadName.c_str() );
            if ( subLeadingCorrelations[i][cent_index][vz_index][aj_index] ) {
              subLeadingCorrelations[i][cent_index][vz_index][aj_index]->SetName( subLeadName.c_str() );
            }
          }
        }
      }
    }
    return 1;
  }

  
  // Function used to find the weighted center
  // for each pt bin for each file - vector<vector<double> >
  // and also creates pt spectra for each file
  std::vector<std::vector<double> > FindPtBinCenter( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& correlations, std::vector<TH1F*>& ptSpectra, binSelector selector ) {
    
    // make the returned object 
    std::vector<std::vector<double> > ptBinCenters;
    ptBinCenters.resize( correlations.size() );

    // first we loop over all histograms, project 
    // down to 1D and add together
    ptSpectra.resize( correlations.size() );
    std::vector<std::vector<TH1F*> > ptBinHolder;
    ptBinHolder.resize( correlations.size() );

    for ( int i = 0; i < correlations.size(); ++i ) {
      for  ( int j = 0; j < correlations[i].size(); ++j ) {
        for ( int k = 0; k < correlations[i][j].size(); ++k ) {
          for ( int l = 0; l < correlations[i][j][k].size(); ++l ) {
            for ( int m = 0; m < selector.nPtBins; ++m ) {
              if ( j == 0 && k == 0 && l == 0 ) {
                ptBinHolder[i].resize(selector.nPtBins);
                std::string tmp = "pt_file_" + patch::to_string(i);
                if ( m == 0 )
                  ptSpectra[i] = new TH1F( tmp.c_str(), tmp.c_str(), binsPt, ptLowEdge, ptHighEdge );
              }

              correlations[i][j][k][l]->GetZaxis()->SetRange();
              ptSpectra[i]->Add( (TH1F*) correlations[i][j][k][l]->ProjectionZ() );

              correlations[i][j][k][l]->GetZaxis()->SetRange( selector.ptBinLowEdge(m), selector.ptBinHighEdge(m) );
              if ( !ptBinHolder[i][m] ) {
                ptBinHolder[i][m] = ((TH1F*) ((TH1F*) correlations[i][j][k][l]->ProjectionZ())->Clone());
                std::string tmp = "pt_file_" + patch::to_string(i) + "_pt_" + patch::to_string(m);
                ptBinHolder[i][m]->SetName( tmp.c_str() );
              }
              else {
                ptBinHolder[i][m]->Add( (TH1F*) correlations[i][j][k][l]->ProjectionZ() );
              }
            }
          }
        }
      }
      
      // now extract the mean values
      ptBinCenters[i].resize(selector.nPtBins);
      for ( int j = 0; j < selector.nPtBins; ++j ) {
        ptBinCenters[i][j] = ptBinHolder[i][j]->GetMean();
        
      }
    }
    
    return ptBinCenters;
  }
  
  // For Correlations
  // Functions to project out the Aj dependence -
  // can either produce a single, Aj independent bin
  // or splits on an ajbin
  void BuildSingleCorrelation( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& correlations, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& reducedCorrelations, binSelector selector, std::string uniqueID ) {
    
    // expand the holder
    reducedCorrelations.resize( correlations.size() );
    
    for ( int i = 0; i < correlations.size(); ++i ) {
      reducedCorrelations[i].resize(correlations[i].size() );
      for ( int j = 0; j < correlations[i].size(); ++j ) {
        reducedCorrelations[i][j].resize(correlations[i][j].size() );
        for ( int k = 0; k < correlations[i][j].size(); ++k ) {
          reducedCorrelations[i][j][k].resize( selector.nPtBins );
          for ( int l = 0; l < correlations[i][j][k].size(); ++l ) {
            for ( int m = 0; m < selector.nPtBins; ++m ) {
              
              // select proper pt range
              correlations[i][j][k][l]->GetZaxis()->SetRange( selector.ptBinLowEdge(m), selector.ptBinHighEdge(m) );
              
              if ( !reducedCorrelations[i][j][k][m] ) {
                std::string tmp = uniqueID + "_corr_file_" + patch::to_string(i) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k);
                reducedCorrelations[i][j][k][m] = (TH2F*) ((TH2F*) correlations[i][j][k][l]->Project3D("YX"))->Clone();
                reducedCorrelations[i][j][k][m]->SetName( tmp.c_str() );
              }
              else {
                reducedCorrelations[i][j][k][m]->Add( (TH2F*) correlations[i][j][k][l]->Project3D("YX") );
              }
            }
          }
        }
      }
    }
  }
  
  void BuildAjSplitCorrelation( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& correlations, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& reducedCorrelationsHigh, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& reducedCorrelationsLow, binSelector selector, int ajBinSplit, std::string uniqueID ) {
    
    // expand the holder
    reducedCorrelationsHigh.resize( correlations.size() );
    reducedCorrelationsLow.resize( correlations.size() );

    for ( int i = 0; i < correlations.size(); ++i ) {
      reducedCorrelationsHigh[i].resize( correlations[i].size() );
      reducedCorrelationsLow[i].resize( correlations[i].size() );
      for ( int j = 0; j < correlations[i].size(); ++j ) {
        reducedCorrelationsHigh[i][j].resize( correlations[i][j].size() );
        reducedCorrelationsLow[i][j].resize( correlations[i][j].size() );
        for ( int k = 0; k < correlations[i][j].size(); ++k ) {
          reducedCorrelationsHigh[i][j][k].resize( selector.nPtBins );
          reducedCorrelationsLow[i][j][k].resize( selector.nPtBins );
          for ( int l = 0; l < correlations[i][j][k].size(); ++l ) {
            for ( int m = 0; m < selector.nPtBins; ++m ) {
              
              // select proper pt range
              correlations[i][j][k][l]->GetZaxis()->SetRange( selector.ptBinLowEdge(m), selector.ptBinHighEdge(m) );
              if ( l >= ajBinSplit ) {
                if ( !reducedCorrelationsHigh[i][j][k][m] ) {
                  std::string tmp = uniqueID +"_corr_aj_low_file_" + patch::to_string(i) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k);
                  reducedCorrelationsHigh[i][j][k][m] = (TH2F*) ((TH2F*) correlations[i][j][k][l]->Project3D("YX"))->Clone();
                  reducedCorrelationsHigh[i][j][k][m]->SetName( tmp.c_str() );
                }
                else {
                  reducedCorrelationsHigh[i][j][k][m]->Add( (TH2F*) correlations[i][j][k][l]->Project3D("YX") );
                }
              }
              else {
                if ( !reducedCorrelationsLow[i][j][k][m] ) {
                  std::string tmp = uniqueID + "_corr_aj_high_file_" + patch::to_string(i) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k);
                  reducedCorrelationsLow[i][j][k][m] = (TH2F*) ((TH2F*) correlations[i][j][k][l]->Project3D("YX"))->Clone();
                  reducedCorrelationsLow[i][j][k][m]->SetName( tmp.c_str() );
                }
                else {
                  reducedCorrelationsLow[i][j][k][m]->Add( (TH2F*) correlations[i][j][k][l]->Project3D("YX") );
                }
              }
              
            }
          }
        }
      }
    }
  }
  
  // Averages over all vz and centralities
  // to show uncorrected signals
  std::vector<std::vector<TH2F*> > AverageCorrelations( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& correlations, binSelector selector, std::string uniqueID ) {
    
    // build the initial holder
    std::vector<std::vector<TH2F*> > averagedCorrelations;
    averagedCorrelations.resize( correlations.size() );
    
    for ( int i = 0; i < correlations.size(); ++i ) {
      averagedCorrelations[i].resize(selector.nPtBins);
      for ( int j = 0; j < correlations[i].size(); ++j ) {
        for ( int k = 0; k < correlations[i][j].size(); ++k ) {
          for ( int l = 0; l < correlations[i][j][k].size(); ++l ) {
            
            if ( !averagedCorrelations[i][l] ) {
              std::string tmp = uniqueID + "_averaged_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              if ( TString(correlations[i][j][k][l]->GetName()).Contains("low") ) {
                tmp = uniqueID + "_averaged_aj_low_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              }
              if ( TString(correlations[i][j][k][l]->GetName()).Contains("high") ) {
                tmp = uniqueID + "_averaged_aj_high_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              }
              
              averagedCorrelations[i][l] = ((TH2F*) correlations[i][j][k][l]->Clone());
              averagedCorrelations[i][l]->SetName( tmp.c_str() );
            }
            else {
              
              averagedCorrelations[i][l]->Add( correlations[i][j][k][l] );
            }
          }
        }
      }
    }
    
    return averagedCorrelations;
  }
  
  
  // For Mixed Events
  // Used to recombine Aj and split in pt
  // to give 2D projections we can turn use
  // to correct the correlations
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > BuildMixedEvents( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& mixedEvents, binSelector selector, std::string uniqueID ) {
    
    std::vector<std::vector<std::vector<std::vector<TH2F*> > > > finalMixedEvents;
    finalMixedEvents.resize( mixedEvents.size() );
    
    for ( int i = 0; i < mixedEvents.size(); ++i ) {
      finalMixedEvents[i].resize( mixedEvents[i].size() );
      for ( int j = 0; j < mixedEvents[i].size(); ++j ) {
        finalMixedEvents[i][j].resize( mixedEvents[i][j].size() );
        for ( int k = 0; k < mixedEvents[i][j].size(); ++k ) {
          finalMixedEvents[i][j][k].resize( selector.nPtBins );
          for ( int l = 0; l < mixedEvents[i][j][k].size(); ++l ) {
            for ( int m = 0; m < selector.nPtBins; ++m ){
              
              // select the proper pt range for each histogram
              mixedEvents[i][j][k][l]->GetZaxis()->SetRange( selector.ptBinLowEdge(m), selector.ptBinHighEdge(m) );
              
              if ( !finalMixedEvents[i][j][k][m] ) {
                finalMixedEvents[i][j][k][m] = (TH2F*) ((TH2F*) mixedEvents[i][j][k][l]->Project3D("YX"))->Clone();
                std::string tmp = uniqueID + "_mix_file_" + patch::to_string(i) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k) + "_pt_" + patch::to_string(m);
                finalMixedEvents[i][j][k][m]->SetName( tmp.c_str() );
              }
              else {
                finalMixedEvents[i][j][k][m]->Add( (TH2F*) mixedEvents[i][j][k][l]->Project3D("YX") );
              }
            }
          }
        }
      }
    }
    
    return finalMixedEvents;
  }

  
  // Used to average the mixed event data to help
  // with the lower statistics
  std::vector<std::vector<TH2F*> > RecombineMixedEvents( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& mixedEvents, binSelector selector, std::string uniqueID ) {
    
    // create the return object
    std::vector<std::vector<TH2F*> > combinedMixedEvents;
    combinedMixedEvents.resize( mixedEvents.size() );
    
    for ( int i = 0; i < mixedEvents.size(); ++i ) {
      combinedMixedEvents[i].resize( 3 );
      
      for ( int j = 0; j < mixedEvents[i].size(); ++j ) {
        for ( int k = 0; k < mixedEvents[i][j].size(); ++k ) {
          for ( int l = 0; l < mixedEvents[i][j][k].size(); ++l ) {
            for ( int m = 0; m < selector.nPtBins; ++m ) {
              
              mixedEvents[i][j][k][l]->GetZaxis()->SetRange( selector.ptBinLowEdge(m), selector.ptBinHighEdge(m) );
              
              if ( m <= 2 ) {
                if ( !combinedMixedEvents[i][m] ) {
                  combinedMixedEvents[i][m] = (TH2F*) ((TH2F*) mixedEvents[i][j][k][l]->Project3D("YX"))->Clone();
                  std::string tmp = uniqueID + "_mix_file_" + patch::to_string(i) + "_pt_" + patch::to_string(m);
                  combinedMixedEvents[i][m]->SetName( tmp.c_str() );
                }
                else {
                  combinedMixedEvents[i][m]->Add( (TH2F*) mixedEvents[i][j][k][l]->Project3D("YX") );
                }
              }
              else {
                if ( !combinedMixedEvents[i][2] ) {
                  combinedMixedEvents[i][2] = (TH2F*) ((TH2F*) mixedEvents[i][j][k][l]->Project3D("YX"))->Clone();
                  std::string tmp = uniqueID + "mix_file_" + patch::to_string(i) + "_pt_" + patch::to_string(m);
                  combinedMixedEvents[i][2]->SetName( tmp.c_str() );
                }
                else {
                  combinedMixedEvents[i][2]->Add( (TH2F*) mixedEvents[i][j][k][l]->Project3D("YX") );
                }
              }
            }
          }
        }
      }
    }
    
    return combinedMixedEvents;
  }
  
  // Used to normalize mixed event histograms so
  // that the maximum bin content = 1
  // version for both the independent mixed events and the weighed averages
  void ScaleMixedEvents( std::vector<std::vector<TH2F*> >& mixedEvents ) {
    
    // scale each histogram
    for ( int i = 0; i < mixedEvents.size(); ++i ) {
      for ( int j = 0; j < mixedEvents[i].size(); ++j ) {
        if ( mixedEvents[i][j]->GetEntries() ) {
          mixedEvents[i][j]->Scale( 1.0 / mixedEvents[i][j]->GetMaximum() );
        }
      }
    }
  }
  
  void ScaleMixedEvents( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& mixedEvents ) {
    
    // scale each histogram
    for ( int i = 0; i < mixedEvents.size(); ++i ) {
      for ( int j = 0; j < mixedEvents[i].size(); ++j ) {
        for ( int k = 0; k < mixedEvents[i][j].size(); ++k ) {
          for ( int l = 0; l < mixedEvents[i][j][k].size(); ++l ) {
            if ( mixedEvents[i][j][k][l]->GetEntries() ) {
              mixedEvents[i][j][k][l]->Scale( 1.0 / mixedEvents[i][j][k][l]->GetMaximum() );
            }
          }
        }
      }
    }
  }
  
  // Used to perform the mixed event division
  // And add up everything into a 2D structure
  // only keeping differntial in file and Pt
  // Has a version for both the averaged and non
  // averaged event mixing
  std::vector<std::vector<TH2F*> > EventMixingCorrection( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& correlations, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& mixedEvents, binSelector selector, std::string uniqueID ) {
    
    // create holder for the output
    std::vector<std::vector<TH2F*> > correctedCorrelations;
    correctedCorrelations.resize(correlations.size() );
    
    for ( int i = 0; i < correlations.size(); ++i ) {
      correctedCorrelations[i].resize( selector.nPtBins );
      for ( int j = 0; j < correlations[i].size(); ++j ) {
        for ( int k = 0; k < correlations[i][j].size(); ++k ) {
          for ( int l = 0; l < correlations[i][j][k].size(); ++l ) {
            
            if ( !correctedCorrelations[i][l] ) {
              std::string tmp = uniqueID + "_corrected_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              if ( TString(correlations[i][j][k][l]->GetName()).Contains("low") ) {
                tmp = uniqueID + "_corrected_aj_low_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              }
              if ( TString(correlations[i][j][k][l]->GetName()).Contains("high") ) {
                tmp = uniqueID + "_corrected_aj_high_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              }
              
              if ( correlations[i][j][k][l]->GetEntries() && mixedEvents[i][j][k][l]->GetEntries() ) {
                TH2F* hTmp = ((TH2F*) correlations[i][j][k][l]->Clone());
                hTmp->Divide( mixedEvents[i][j][k][l] );
                correctedCorrelations[i][l] = (TH2F*) hTmp->Clone();
                correctedCorrelations[i][l]->SetName( tmp.c_str() );
              }
            }
            else {
              if ( correlations[i][j][k][l]->GetEntries() && mixedEvents[i][j][k][l]->GetEntries() ) {
                TH2F* hTmp = ((TH2F*) correlations[i][j][k][l]->Clone());
                hTmp->Divide( mixedEvents[i][j][k][l] );
                correctedCorrelations[i][l]->Add( hTmp );
              }
            }
          }
        }
      }
    }
    return correctedCorrelations;
  }
  
  std::vector<std::vector<TH2F*> > EventMixingCorrection( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& correlations, std::vector<std::vector<TH2F*> >& mixedEvents, binSelector selector, std::string uniqueID ) {
    
    // create holder for the output
    std::vector<std::vector<TH2F*> > correctedCorrelations;
    correctedCorrelations.resize(correlations.size() );
    
    for ( int i = 0; i < correlations.size(); ++i ) {
      correctedCorrelations[i].resize( selector.nPtBins );
      for ( int j = 0; j < correlations[i].size(); ++j ) {
        for ( int k = 0; k < correlations[i][j].size(); ++k ) {
          for ( int l = 0; l < correlations[i][j][k].size(); ++l ) {
            if ( !correctedCorrelations[i][l] ) {
              std::string tmp = uniqueID + "_corrected_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              if ( TString(correlations[i][j][k][l]->GetName()).Contains("low") ) {
                tmp = uniqueID + "_corrected_aj_low_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              }
              if ( TString(correlations[i][j][k][l]->GetName()).Contains("high") ) {
                tmp = uniqueID + "_corrected_aj_high_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              }
              
              if ( correlations[i][j][k][l]->GetEntries() ) {
                TH2F* hTmp = ((TH2F*) correlations[i][j][k][l]->Clone());
                if ( l <= 2 && mixedEvents[i][l]->GetEntries() ) {
                  hTmp->Divide( mixedEvents[i][l] );
                }
                else if ( mixedEvents[i][2]->GetEntries() )  {
                  hTmp->Divide( mixedEvents[i][2] );
                }
                else {
                  __ERR("Did not have any mixed event data to correct with")
                  continue;
                }
                correctedCorrelations[i][l] = (TH2F*) hTmp->Clone();
                correctedCorrelations[i][l]->SetName( tmp.c_str() );
              }
            }
            else {
              if ( correlations[i][j][k][l]->GetEntries() ) {
                TH2F* hTmp = ((TH2F*) correlations[i][j][k][l]->Clone());
                if ( l <= 2 && mixedEvents[i][l]->GetEntries() ) {
                  hTmp->Divide( mixedEvents[i][l] );
                }
                else if ( mixedEvents[i][2]->GetEntries() )  {
                  hTmp->Divide( mixedEvents[i][2] );
                }
                else {
                  __ERR("Did not have any mixed event data to correct with")
                  continue;
                }
                correctedCorrelations[i][l]->Add( hTmp );
              }
            }
          }
        }
      }
    }
    return correctedCorrelations;
  }
  
  // Used to extract 1D projections from
  // the 2D histograms - allows for setting
  // ranges for the projection ( e.g. projecting
  // only the near side of the dPhi range in a dEta
  // projection )
  std::vector<std::vector<TH1F*> > ProjectDphi( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, std::string uniqueID, bool restrictDeta ) {
    
    // build the return vector
    std::vector<std::vector<TH1F*> > projections;
    projections.resize( correlation2d.size() );
    
    // now loop over every 2d histogram and project
    for ( int i = 0; i < correlation2d.size(); ++i ) {
      projections[i].resize( correlation2d[i].size() );
      for ( int j = 0; j < correlation2d[i].size(); ++j ) {
        
        // new name for the projection
        std::string tmp = uniqueID + "_dphi_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
        
        if ( restrictDeta ) {
          correlation2d[i][j]->GetXaxis()->SetRangeUser( selector.phi_projection_eta_bound_low, selector.phi_projection_eta_bound_high );
        }
        
        projections[i][j] = (TH1F*) correlation2d[i][j]->ProjectionY();
        projections[i][j]->SetName( tmp.c_str() );
        
        if ( restrictDeta ) {
          correlation2d[i][j]->GetXaxis()->SetRange();
        }
        
      }
    }
    
    return projections;
  }
  
  std::vector<std::vector<TH1F*> > ProjectDphiNearMinusFar( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, std::string uniqueID, bool restrictDeta ) {
    
    // build the return vector
    std::vector<std::vector<TH1F*> > projections;
    projections.resize( correlation2d.size() );
    
    // now loop over every 2d histogram and project
    for ( int i = 0; i < correlation2d.size(); ++i ) {
      projections[i].resize( correlation2d[i].size() );
      for ( int j = 0; j < correlation2d[i].size(); ++j ) {
        
        // new name for the projection
        std::string tmp = uniqueID + "_dphi_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
        
        // now find the bins for the subtraction
        selector.SetHistogramBins( correlation2d[i][j] );
        double etaMin = selector.dEtaLow;
        double etaMax = selector.dEtaHigh;
        double etaBins = selector.bindEta;
        double etaBinWidth = ( etaMax - etaMin ) / etaBins;
        
        // now we set the actual range of interest
        // but we can restrict the region of interest if restrictDeta is
        double acceptEtaMin;
        double acceptEtaMax;
        if ( restrictDeta ) {
          acceptEtaMin = selector.dEtaAcceptanceLow+0.4;
          acceptEtaMax = selector.dEtaAcceptanceHigh-0.4;
        }
        else { // if not restricted, the bounds are set by the kinematics of the jet
          acceptEtaMin = selector.dEtaAcceptanceLow;
          acceptEtaMax = selector.dEtaAcceptanceHigh;
        }
        
        // and define our working range between eta min and eta max
        double etaRange = ( acceptEtaMax - acceptEtaMin );
        
        // now we define three regions -
        // region 2 - inner range - the "near side"
        // region 2 has twice the width of range 1 or 3
        // region 1 & 3 - two outer regions - the remainder
        // of the histogram on either side
        // here we find the bin numbers corresponding
        // to these ranges
        
        // region one will be from [etaMin -> etaMin + etaRange/4.0 )
        // region two will be from [etaMin + etaRange -> etaMin + 3.0*etaRange/4.0]
        // region three will be from ( etaMin + 3.0*etaRange/4.0 -> etaMax]
        double bound1 = acceptEtaMin;
        double bound2 = acceptEtaMin + etaRange/4.0;
        double bound3 = acceptEtaMin + 3.0*etaRange/4.0;
        double bound4 = acceptEtaMax;
        
        
        int range1Low = 0;
        int range1High = 0;
        int range2Low = 0;
        int range2High = 0;
        int range3Low = 0;
        int range3High = 0;
        
        // now we will search for each bin
        for ( int i = 1; i <= etaBins; ++i ) {
          double binLowEdge = etaMin + ( i - 1 )*etaBinWidth;
          double binUpEdge = etaMin + ( i )*etaBinWidth;
          
          if ( bound1 >= binLowEdge && bound1 < binUpEdge ) {
            range1Low = i;
          }
          if ( bound2 >= binLowEdge && bound2 < binUpEdge ) {
            range1High = i-1;
            range2Low = i;
          }
          if ( bound3 >= binLowEdge && bound3 < binUpEdge ) {
            range2High = i;
            range3Low = i+1;
          }
          if ( bound4 >= binLowEdge && bound4 < binUpEdge ) {
            range3High = i;
          }
        }
        
        // now do the projections
        correlation2d[i][j]->GetXaxis()->SetRange( range2Low, range2High );
        projections[i][j] = (TH1F*) correlation2d[i][j]->ProjectionY();
        projections[i][j]->SetName( tmp.c_str() );
        
        // build the subtraction histogram
        correlation2d[i][j]->GetXaxis()->SetRange( range1Low, range1High );
        TH1F* sub_tmp = (TH1F*) ((TH1F*)  correlation2d[i][j]->ProjectionY())->Clone();
        correlation2d[i][j]->GetXaxis()->SetRange( range3Low, range3High );
        sub_tmp->Add( (TH1F*) correlation2d[i][j]->ProjectionY() );
        
        // scale the subtraction histogram by the relative number of bins
        sub_tmp->Scale( (range2High-range2Low)/( (range1High-range1Low) + (range3High - range3Low) ) );
        
        // subtract
        projections[i][j]->Add( sub_tmp, -1 );
        
      }
    }
    
    return projections;
    
  }
  
  std::vector<std::vector<TH1F*> > ProjectDeta( std::vector<std::vector<TH2F*> >& correlation2d, binSelector selector, std::string uniqueID, bool restrictDphi ) {
    // build the return vector
    std::vector<std::vector<TH1F*> > projections;
    projections.resize( correlation2d.size() );
    
    // now loop over every 2d histogram and project
    for ( int i = 0; i < correlation2d.size(); ++i ) {
      projections[i].resize( correlation2d[i].size() );
      for ( int j = 0; j < correlation2d[i].size(); ++j ) {
        
        // new name for the projection
        std::string tmp = uniqueID + "_deta_file_" + patch::to_string(i) + "_pt_" + patch::to_string(j);
        
        if ( restrictDphi ) {
          correlation2d[i][j]->GetYaxis()->SetRangeUser( selector.eta_projection_phi_bound_low, selector.eta_projection_phi_bound_high );
        }
        
        projections[i][j] = (TH1F*) correlation2d[i][j]->ProjectionX();
        projections[i][j]->SetName( tmp.c_str() );
        
        if ( restrictDphi ) {
          correlation2d[i][j]->GetYaxis()->SetRange();
        }
        
      }
    }
    
    return projections;

  }
  
  
  
  
  
  // *****************************
  // HISTOGRAM PRINTING AND SAVING
  // *****************************
  
  // Used to print out 2D plots ( correlations, mixed events )
  void Print2DHistograms( std::vector<TH2F*>& histograms, std::string outputDir, std::string analysisName, binSelector selector ) {
    
    // First, make the output directory if it doesnt exist
    boost::filesystem::path dir( outputDir.c_str() );
    boost::filesystem::create_directories( dir );
    
    for ( int i = 0; i < histograms.size(); ++i ) {
      histograms[i]->GetXaxis()->SetTitle("#Delta#eta");
      histograms[i]->GetYaxis()->SetTitle("#Delta#phi");
      histograms[i]->SetTitle( selector.ptBinString[i].c_str() );
      
      std::string tmp = outputDir + "/" + analysisName + "_" + patch::to_string(i) + ".pdf";
      
      TCanvas c1;
      histograms[i]->Draw("surf1");
      c1.SaveAs( tmp.c_str() );
    }
    
  }
  
  // Used to print out 2D plots ( correlations, mixed events )
  // However, this one restricts the eta range shown to what
  // is set in selector
  void Print2DHistogramsEtaRestricted( std::vector<TH2F*>& histograms, std::string outputDir, std::string analysisName, binSelector selector ) {
    
    // First, make the output directory if it doesnt exist
    boost::filesystem::path dir( outputDir.c_str() );
    boost::filesystem::create_directories( dir );
    
    for ( int i = 0; i < histograms.size(); ++i ) {
      histograms[i]->GetXaxis()->SetTitle("#Delta#eta");
      histograms[i]->GetXaxis()->SetRangeUser(selector.phi_projection_eta_bound_low, selector.phi_projection_eta_bound_high );
      histograms[i]->GetYaxis()->SetTitle("#Delta#phi");
      histograms[i]->SetTitle( selector.ptBinString[i].c_str() );
      
      std::string tmp = outputDir + "/" + analysisName + "_" + patch::to_string(i) + ".pdf";
      
      TCanvas c1;
      histograms[i]->Draw("surf1");
      c1.SaveAs( tmp.c_str() );
    }
    
  }
  
  // Print out individual dPhi histograms
  void Print1DHistogramsDphi( std::vector<TH1F*>& histograms, std::string outputDir, std::string analysisName, binSelector selector ) {
  
    // First, make the output directory if it doesnt exist
    boost::filesystem::path dir( outputDir.c_str() );
    boost::filesystem::create_directories( dir );
    
    for ( int i = 0; i < histograms.size(); ++i ) {
      histograms[i]->GetXaxis()->SetTitle("#Delta#phi");
      histograms[i]->SetTitle( selector.ptBinString[i].c_str() );
      
      std::string tmp = outputDir + "/" + analysisName + "_" + patch::to_string(i) + ".pdf";
      
      TCanvas c1;
      histograms[i]->Draw();
      c1.SaveAs( tmp.c_str() );
    }
    
  }
  
  // Print out dPhi histograms, each pt bin
  // Seperately, but each file overlaid
  void Print1DHistogramsOverlayedDphi( std::vector<std::vector<TH1F*> >& histograms, std::string outputDir, std::vector<std::string> analysisName, binSelector ) {
    
    // First, make the output directory if it doesnt exist
    boost::filesystem::path dir( outputDir.c_str() );
    boost::filesystem::create_directories( dir );
    
    TCanvas c1;
    for ( int i = 0; i < histograms[0].size(); ++i ) {
      for ( int j = 0; j < histograms.size(); ++j ) {
      
        histograms[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        histograms[j][i]->SetTitle( selector.ptBinString[i].c_str() );
        histograms[j][i]->SetLineColor( i+1 );
        histograms[j][i]->SetMarkerStyle( i+20 );
        histograms[j][i]->SetMarkerColor( i+1 );
        histograms[j][i]->SetMarkerSize( 2 );
      
        if ( j == 0 ) {
        histograms[j][i]->Draw();
        }
        else {
          histograms[j][i]->Draw("same");
        }
        
      }
      std::string tmp = outputDir + "/" + analysisName + "_" + patch::to_string(i) + ".pdf";
      c1.SaveAs( tmp.c_str() );
    }
    
  }
  
  // Print out individual dEta histograms
  void Print1DHistogramsDeta( std::vector<TH1F*>& histograms, std::string outputDir, std::string analysisName, binSelector selector ) {
    
    // First, make the output directory if it doesnt exist
    boost::filesystem::path dir( outputDir.c_str() );
    boost::filesystem::create_directories( dir );
    
    for ( int i = 0; i < histograms.size(); ++i ) {
      histograms[i]->GetXaxis()->SetTitle("#Delta#eta");
      histograms[i]->SetTitle( selector.ptBinString[i].c_str() );
      
      std::string tmp = outputDir + "/" + analysisName + "_" + patch::to_string(i) + ".pdf";
      
      TCanvas c1;
      histograms[i]->Draw();
      c1.SaveAs( tmp.c_str() );
    }

    
  }
  
  // Print out dPhi histograms, each pt bin
  // Seperately, but each file overlaid
  void Print1DHistogramsOverlayedDeta( std::vector<std::vector<TH1F*> >& histograms, std::string outputDir, std::vector<std::string> analysisName, binSelector ) {
    
    // First, make the output directory if it doesnt exist
    boost::filesystem::path dir( outputDir.c_str() );
    boost::filesystem::create_directories( dir );
    
    TCanvas c1;
    for ( int i = 0; i < histograms[0].size(); ++i ) {
      for ( int j = 0; j < histograms.size(); ++j ) {
        
        histograms[j][i]->GetXaxis()->SetTitle("#Delta#phi");
        histograms[j][i]->SetTitle( selector.ptBinString[i].c_str() );
        histograms[j][i]->SetLineColor( i+1 );
        histograms[j][i]->SetMarkerStyle( i+20 );
        histograms[j][i]->SetMarkerColor( i+1 );
        histograms[j][i]->SetMarkerSize( 2 );
        
        if ( j == 0 ) {
          histograms[j][i]->Draw();
        }
        else {
          histograms[j][i]->Draw("same");
        }
        
      }
      std::string tmp = outputDir + "/" + analysisName + "_" + patch::to_string(i) + ".pdf";
      c1.SaveAs( tmp.c_str() );
    }

    
  }

  
} // end namespace
