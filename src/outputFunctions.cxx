// Implementation of basic functions
// used in the output workflow of
// the jet hadron correlation analysis

#include "outputFunctions.hh"
#include "corrParameters.hh"
#include "histograms.hh"

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
      std::cout<<"got nevents"<<std::endl;
      // push back the vectors
      leadingCorrelations.push_back( std::vector<std::vector<std::vector<TH3F*> > >() );
      subLeadingCorrelations.push_back( std::vector<std::vector<std::vector<TH3F*> > >() );
      std::cout<<"pushing back first"<<std::endl;
      for ( int j = selector.ajLow; j <= selector.ajHigh; ++j ) {
        
        int aj_index = j - selector.ajLow;
        
        // push back the vectors
        leadingCorrelations[i].push_back( std::vector<std::vector<TH3F*> >() );
        subLeadingCorrelations[i].push_back( std::vector<std::vector<TH3F*> >() );
        std::cout<<"pushing back second"<<std::endl;
        for ( int k = selector.centLow; k <= selector.centHigh; ++k ) {
          
          int cent_index = k - selector.centLow;
          
          // push back the vectors
          leadingCorrelations[i][aj_index].push_back( std::vector<TH3F*>() );
          subLeadingCorrelations[i][aj_index].push_back( std::vector<TH3F*>() );
          std::cout<<"pushing back third"<<std::endl;
          for ( int l = selector.vzLow; l <= selector.vzHigh; ++l ) {
            
            int vz_index = l - selector.vzLow;
            
            // build the in-file histogram names
            std::string leadName = "lead_aj_" + patch::to_string(j) + "_cent_" + patch::to_string(k) + "_vz_" + patch::to_string(l);
            std::string subLeadName = "sub_aj_" + patch::to_string(j) + "_cent_" + patch::to_string(k) + "_vz_" + patch::to_string(l);
            std::cout<<"getting histograms"<<std::endl;
            
            // get the correlation histograms
            leadingCorrelations[i][aj_index][cent_index].push_back( (TH3F*) filesIn[i]->Get( leadName.c_str() ) );
            std::cout<<"got lead"<<std::endl;
            subLeadingCorrelations[i][aj_index][cent_index].push_back( (TH3F*) filesIn[i]->Get( subLeadName.c_str() ) );
            std::cout<<"got sub"<<std::endl;
            // check to make sure it was successful
            if ( !leadingCorrelations[i][aj_index][cent_index][vz_index] ) {
              std::string errorMsg = "Couldn't read in leading correlation: " + patch::to_string(i) + " " + patch::to_string(j) + " " + patch::to_string(k) + " " + patch::to_string(l);
              __ERR( errorMsg.c_str() )
              continue;
            }
            if ( !subLeadingCorrelations[i][aj_index][cent_index][vz_index] ) {
              std::string errorMsg = "Couldn't read in subleading correlation: " + patch::to_string(i) + " " + patch::to_string(j) + " " + patch::to_string(k) + " " + patch::to_string(l);
              __ERR( errorMsg.c_str() )
              continue;
            }
            
            //now differentiate the names by file
            leadName = "file_" + patch::to_string(i) + leadName;
            subLeadName = "file_" + patch::to_string(i) + subLeadName;
            
            leadingCorrelations[i][aj_index][cent_index][vz_index]->SetName( leadName.c_str() );
            subLeadingCorrelations[i][aj_index][cent_index][vz_index]->SetName( subLeadName.c_str() );
            
          }
        }
      }
    }
  }
  

  
  
  
} // end namespace
