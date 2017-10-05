#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
using namespace Rcpp;

//' Convert phase files to genotype files.
//' 
//' Simply sums every pair of rows together.
//' 
//' A phase file has format similar to SNP files, expect the genotype on each allele are listed on two pairs of rows.
//' Each individual has therefore two rows, one for each allele, see following example:
//'
//' \code{Genotype file:}
//' \tabular{lccccc}{
//'  1003 \tab 0 \tab 1 \tab 1 \tab 2 \tab 0 \cr
//'  1004 \tab 1 \tab 1 \tab 0 \tab 1 \tab 0 \cr
//' }
//'
//' \code{Phase file:}
//' \tabular{lccccc}{
//'  1003 \tab 0 \tab 0 \tab 1 \tab 1 \tab 0 \cr
//'  1003 \tab 0 \tab 1 \tab 0 \tab 1 \tab 0 \cr
//'  1004 \tab 0 \tab 1 \tab 0 \tab 0 \tab 0 \cr
//'  1004 \tab 1 \tab 0 \tab 0 \tab 1 \tab 0 \cr
//' }
//'
//' \strong{Missing values:} Values after summing less than 0 or greater than 2 are assumed as missing and replaced with \code{na}.
//' Value of this range can be changed with argument \code{range}.
//'
//' @param fnin Filename of input file, every two rows are for same animal.
//' @param fnout Filename of output file.
//' @param nlines Integer, maximum number of pairs of lines to convert. 
//'               When \code{-1} (default), no maximum.
//' @param na Value to use for missing values.
//' @param range Integer vector of range of allowable values. 
//'              Values in either input or output that are strictly smaller than 
//'              first element or strictly larger than second element are replaced
//'              with \code{na}.
//' @param idwidth Width of ID column.
//' @param precision Number of decimals to print genotypes; 
//'                  when \code{-1} (default) genotypes are printed with mixed 
//'                  precision. When \code{0}, genotypes are rounded to whole integers.
//' @return Number of rows written.
//' @backref src/convert_phases.cpp
//' @export
// [[Rcpp::export]]
long convert_phases(std::string fnin, 
                    std::string fnout, 
                    long nlines = -1, //Rcpp::Nullable<Rcpp::IntegerVector> nlines = R_NilValue,
                    long na = 9,
                    Rcpp::IntegerVector range = Rcpp::IntegerVector::create(0, 2),
                    int idwidth = 4,
                    int precision = -1) 
{
  
  long cnt;
  std::string id1;
  std::string id2;
  float g;
  std::string line1;
  std::string line2;
  std::vector<float> phase1;
  std::vector<float> phase2;
  
  std::ifstream infile(fnin.c_str()); //Opens the file. c_str is mandatory here so that ifstream accepts the string path
  std::ofstream outfile(fnout.c_str()); 
  
  if (!infile.is_open())
    Rcpp::stop("Could not open '%s'.", fnin);
  if (!outfile.is_open())
    Rcpp::stop("Could not open '%s' for writing.", fnout);
    
  cnt = 0;
  
  while(std::getline(infile, line1)) {
    if (!std::getline(infile, line2))
      break;
    cnt++;
    
    // Consider using iss1.reserve  with number of columns. Would require wrapping function to determine number of columns beforehand.
    std::istringstream iss1(line1);
    std::istringstream iss2(line2);
    iss1 >> id1;
    iss2 >> id2;
    
    if (id1 != id2)
      Rcpp::warning("Mismatching ID's (%s and %s) at approx. line %d.", id1, id2, cnt * 2 - 1);
    
    std::vector<float> phase1(( std::istream_iterator<float>( iss1 )), std::istream_iterator<float>() ); // Extra set of parantheses around first argument to count C++'s most vexing parse.
    std::vector<float> phase2(( std::istream_iterator<float>( iss2 )), std::istream_iterator<float>() );
    
    if (phase1.size() != phase2.size()) {
      Rcpp::warning("Mismatching lengths of phases (%d vs %d) for id %s at approx. line %d. Will skip", phase1.size(), phase2.size(), id1, cnt * 2 - 1);
      continue;
    }
    
    
    outfile.width(idwidth);
    outfile.setf(std::ios::left, std::ios::adjustfield);
    outfile << id1 << " ";
    for (size_t  i = 0; i < phase1.size(); i++) {
      if (phase1[i] < range[0] || phase1[i] > range[1]) {
        g = na;
      } else if (phase2[i] < range[0] || phase2[i] > range[1]) {
        g = na;
      } else {
        g = phase1[i] + phase2[i];
        if (g < range[0] || g > range[1])
          g = na;
      }
      if (precision > -1) {
        outfile.precision(precision);
        outfile << std::fixed;
      }
      outfile << g << " ";
    }
    outfile << std::endl;
    
    if (nlines == cnt)
      break;
    
  }
  
  outfile.close();
  infile.close();
  
  return cnt;
}