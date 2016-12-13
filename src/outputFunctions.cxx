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
  void binSelector::SetHistogramBins( TH3F* h ) {
    
    bindEta = h->GetXaxis()->GetNbins();
    dEtaLow = h->GetXaxis()->GetBinLowEdge(1);
    dEtaHigh = h->GetXaxis()->GetBinUpEdge( h->GetXaxis()->GetNbins() );
    bindPhi = h->GetYaxis()->GetNbins();
    dPhiLow = h->GetYaxis()->GetBinLowEdge(1);
    dPhiHigh = h->GetYaxis()->GetBinUpEdge( h->GetYaxis()->GetNbins() );
    
  }
  
  
  // Function used to read in histograms from
  // the files passed in - it returns the correlations,
  // and the number of events, and selects using the centralities,
  // vz bin range, and aj ranges passed in via binSelector
  void ReadInFiles(std::vector<TFile*>& filesIn, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& leadingCorrelations, std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& subLeadingCorrelations, std::vector<TH3F*>& nEvents, binSelector selector ) {
    
    // loop over all files and all aj, centrality, and vz bins
    // to return a 4D vector of histograms
    for ( int i = 0; i < filesIn.size(); ++i ) {
      // tell the user what is going on
      std::string outMsg = "Reading in file " + patch::to_string(i);
      __OUT(outMsg.c_str() )
      
      // for each file, get the number of events
      nEvents.push_back( (TH3F*) filesIn[i]->Get("nevents") );
      // rename because root can't handle simple crap
      std::string tmpName = "nevents_" + patch::to_string(i);
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
            
            
            // get the correlation histograms
            leadingCorrelations[i][cent_index][vz_index].push_back( (TH3F*) filesIn[i]->Get( leadName.c_str() ) );
            
            subLeadingCorrelations[i][cent_index][vz_index].push_back( (TH3F*) filesIn[i]->Get( subLeadName.c_str() ) );
            
            // check to make sure it was successful
            if ( !leadingCorrelations[i][cent_index][vz_index][aj_index] ) {
              std::string errorMsg = "Couldn't read in leading correlation: " + patch::to_string(i) + " " + patch::to_string(j) + " " + patch::to_string(k) + " " + patch::to_string(l);
              __ERR( errorMsg.c_str() )
              continue;
            }
            if ( !subLeadingCorrelations[i][cent_index][vz_index][aj_index] ) {
              std::string errorMsg = "Couldn't read in subleading correlation: " + patch::to_string(i) + " " + patch::to_string(j) + " " + patch::to_string(k) + " " + patch::to_string(l);
              __ERR( errorMsg.c_str() )
              continue;
            }
            
            //now differentiate the names by file
            leadName = "file_" + patch::to_string(i) + leadName;
            subLeadName = "file_" + patch::to_string(i) + subLeadName;
            
            leadingCorrelations[i][cent_index][vz_index][aj_index]->SetName( leadName.c_str() );
            subLeadingCorrelations[i][cent_index][vz_index][aj_index]->SetName( subLeadName.c_str() );
            
          }
        }
      }
    }
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
  void BuildSingleCorrelation( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& correlations, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& reducedCorrelations, binSelector selector ) {
    
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
                std::string tmp = "corr_file_" + patch::to_string(i) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k);
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
  
  void BuildAjSplitCorrelation( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& correlations, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& reducedCorrelationsHigh, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& reducedCorrelationsLow, binSelector selector, int ajBinSplit ) {
    
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
              if ( l <= ajBinSplit ) {
                if ( !reducedCorrelationsHigh[i][j][k][m] ) {
                  std::string tmp = "corr_aj_low_file_" + patch::to_string(i) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k);
                  reducedCorrelationsHigh[i][j][k][m] = (TH2F*) ((TH2F*) correlations[i][j][k][l]->Project3D("YX"))->Clone();
                  reducedCorrelationsHigh[i][j][k][m]->SetName( tmp.c_str() );
                }
                else {
                  reducedCorrelationsHigh[i][j][k][m]->Add( (TH2F*) correlations[i][j][k][l]->Project3D("YX") );
                }
              }
              else {
                if ( !reducedCorrelationsLow[i][j][k][m] ) {
                  std::string tmp = "corr_aj_high_file_" + patch::to_string(i) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k);
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
  
  
  // For Mixed Events
  // Used to recombine Aj and split in pt
  // to give 2D projections we can turn use
  // to correct the correlations
  std::vector<std::vector<std::vector<std::vector<TH2F*> > > > BuildMixedEvents( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& mixedEvents, binSelector selector ) {
    
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
                std::string tmp = "mix_file_" + patch::to_string(i) + "_cent_" + patch::to_string(j) + "_vz_" + patch::to_string(k) + "_pt_" + patch::to_string(m);
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
  std::vector<std::vector<TH2F*> > RecombineMixedEvents( std::vector<std::vector<std::vector<std::vector<TH3F*> > > >& mixedEvents, binSelector selector ) {
    
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
                  std::string tmp = "mix_file_" + patch::to_string(i) + "_pt_" + patch::to_string(m);
                  combinedMixedEvents[i][m]->SetName( tmp.c_str() );
                }
                else {
                  combinedMixedEvents[i][m]->Add( (TH2F*) mixedEvents[i][j][k][l]->Project3D("YX") );
                }
              }
              else {
                if ( !combinedMixedEvents[i][2] ) {
                  combinedMixedEvents[i][2] = (TH2F*) ((TH2F*) mixedEvents[i][j][k][l]->Project3D("YX"))->Clone();
                  std::string tmp = "mix_file_" + patch::to_string(i) + "_pt_" + patch::to_string(m);
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
          std::cout<<"max: "<< mixedEvents[i][j]->GetMaximum() << std::endl;
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
  std::vector<std::vector<TH2F*> > EventMixingCorrection( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& correlations, std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& mixedEvents, binSelector selector) {
    
    // create holder for the output
    std::vector<std::vector<TH2F*> > correctedCorrelations;
    correctedCorrelations.resize(correlations.size() );
    
    for ( int i = 0; i < correlations.size(); ++i ) {
      correctedCorrelations[i].resize( selector.nPtBins );
      for ( int j = 0; j < correlations[i].size(); ++j ) {
        for ( int k = 0; k < correlations[i][j].size(); ++k ) {
          for ( int l = 0; l < correlations[i][j][k].size(); ++l ) {
            
            if ( !correctedCorrelations[i][l] ) {
              std::string tmp = "corrected_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              if ( correlations[i][j][k][l]->GetName().Contains("low") ) {
                tmp = "corrected_aj_low_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              }
              if ( correlations[i][j][k][l]->GetName().Contains("high") ) {
                tmp = "corrected_aj_high_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              }
              
              if ( correlations[i][j][k][l]->GetEntries() && mixedEvents[i][j][k][l]->GetEntries() ) {
                correctedCorrelations[i][l] = ((TH2F*) ((TH2F*) correlations[i][j][k][l]->Project3D("YX"))->Clone())->Divide( mixedEvents[i][j][k][l] );
                correctedCorrelations[i][l]->SetName( tmp.c_str() );
              }
            }
            else {
              if ( correlations[i][j][k][l]->GetEntries() && mixedEvents[i][j][k][l]->GetEntries() ) {
                correctedCorrelations[i][l]->Add( ((TH2F*) ((TH2F*) correlations[i][j][k][l]->Project3D("YX"))->Clone())->Divide( mixedEvents[i][j][k][l] ) );
              }
            }
          }
        }
      }
    }
    return correctedCorrelations;
  }
  
  std::vector<std::vector<TH2F*> > EventMixingCorrection( std::vector<std::vector<std::vector<std::vector<TH2F*> > > >& correlations, std::vector<std::vector<TH2F*> >& mixedEvents, binSelector selector ) {
    
    // create holder for the output
    std::vector<std::vector<TH2F*> > correctedCorrelations;
    correctedCorrelations.resize(correlations.size() );
    
    for ( int i = 0; i < correlations.size(); ++i ) {
      correctedCorrelations[i].resize( selector.nPtBins );
      for ( int j = 0; j < correlations[i].size(); ++j ) {
        for ( int k = 0; k < correlations[i][j].size(); ++k ) {
          for ( int l = 0; l < correlations[i][j][k].size(); ++l ) {
            
            if ( !correctedCorrelations[i][l] ) {
              std::string tmp = "corrected_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              if ( correlations[i][j][k][l]->GetName().Contains("low") ) {
                tmp = "corrected_aj_low_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              }
              if ( correlations[i][j][k][l]->GetName().Contains("high") ) {
                tmp = "corrected_aj_high_file_" + patch::to_string(i) + "_pt_" + patch::to_string(l);
              }
              
              if ( correlations[i][j][k][l]->GetEntries() && mixedEvents[i][l]->GetEntries() ) {
                correctedCorrelations[i][l] = ((TH2F*) ((TH2F*) correlations[i][j][k][l]->Project3D("YX"))->Clone())->Divide( mixedEvents[i][l] );
                correctedCorrelations[i][l]->SetName( tmp.c_str() );
              }
            }
            else {
              if ( correlations[i][j][k][l]->GetEntries() && mixedEvents[i][l]->GetEntries() ) {
                correctedCorrelations[i][l]->Add( ((TH2F*) ((TH2F*) correlations[i][j][k][l]->Project3D("YX"))->Clone())->Divide( mixedEvents[i][l] ) );
              }
            }
          }
        }
      }
    }
    return correctedCorrelations;
  }
  

  
  
  
  
  // *****************************
  // HISTOGRAM PRINTING AND SAVING
  // *****************************
  
  // Used to print out 2D plots ( correlations, mixed events )
  void Print2DHistograms( std::vector<TH2F*>& histograms, std::string outputDir, std::string analysisName, binSelector selector ) {
    
    // First, make the output directory if it doesnt exist
    boost::filesystem::path dir( outputDir.c_str() );
    if ( !boost::filesystem::create_directories( dir ) ) {
      std::cout << "Success" << "\n";
    }
    
    
  }
  
  
} // end namespace
