// testing nick's new correlation implementations

// first test the parameters
#include "corrParameters.hh"
#include "corrFunctions.hh"
#include <iostream>

int main() {
	std::cout<<"vzRange: "<< corrAnalysis::vzRange<<std::endl;
	std::cout<<"vzBins: "<< corrAnalysis::binsVz<<std::endl;
	std::cout<<"dVz: "<< corrAnalysis::dVz<<std::endl;
	std::cout<<"vzMin/Max: "<< corrAnalysis::vzLowEdge<<" "<<corrAnalysis::vzLowEdge<<std::endl;
	
	corrAnalysis::EndSummaryJet( 10.0, 10.0, 5 );
	
	return 0;
}