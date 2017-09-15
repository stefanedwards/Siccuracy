#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
#include <string>
#include <vector>
using namespace Rcpp;



// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
int convert_phase(std::string fin, std::string fout) 
{
  
  
  std::ifstream myfile(fin.c_str()); //Opens the file. c_str is mandatory here so that ifstream accepts the string path
  
  std::string line1;
  std::string line2;
  char id;
  std::vector<int> phase1;
  std::vector<int> phase2;
  
  // For copying to Rcout
  std::stringstream sout;
  
  while (std::getline(myfile, line1)) {
    // Read in pairs of lines
    if (!std::getline(myfile, line2))
      break;
    
    Rcout << std::endl;
    Rcout << "line 1: " << line1 << std::endl;
    Rcout << "line 2: " << line2 << std::endl;
    
    
    std::istringstream iss1(line1);
    std::istringstream iss2(line2);
    iss1 >> id;
    iss2 >> id;
      
    std::vector<int> phase1(( std::istream_iterator<int>( iss1 )), std::istream_iterator<int>() ); // Extra set of parantheses around first argument to count C++'s most vexing parse.
    std::vector<int> phase2(( std::istream_iterator<int>( iss2 )), std::istream_iterator<int>() );
    
    Rcpp::Rcout << "id for phase 2 is "<< id << std::endl;
    Rcpp::Rcout << "something something phase 1: " << phase1.size() << std::endl;
    
    sout.str(std::string()); // reset, or clear??
    Rcpp::Rcout << "phase 1 is "; // << phase1 << std::endl;
    std::copy(phase1.begin(), phase1.end(), std::ostream_iterator<int>(sout, " "));
    Rcout << sout.str() << std::endl;
    sout.clear();
    Rcpp::Rcout << "phase 2 is "; // << phase1 << std::endl;
    std::copy(phase2.begin(), phase2.end(), std::ostream_iterator<int>(sout, " "));
    Rcout << sout.str() << std::endl;
    sout.clear();

  }

  return(0);
}
