#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm> // for copy
#include <iterator> // for ostream_iterator
using namespace Rcpp;

//' Convert PLINK recoded A to SNP file.
//' 
//' Facilitates converting a PLINK binary file to simplified SNP file format.
//' Requires using PLINK to recode it to the \code{A} format by using command line \code{plink -bfile <file stem> --recode A}.
//' This function then swiftly strips of first 6 columns (family ID, sample ID, paternal ID, maternal ID, sex, phenotypic record) 
//' and inserts an integer-based ID column. \code{NA}'s outputted from PLINK are replaced with \code{na} argument.
//' 
//' @param fnraw Plink output filename. Most likely \code{plink.raw} if PLINK command line argument \code{--out} is not used.
//' @param fnout Filename of new file.
//' @param newID Integer scalar (default \code{1}) for first integer ID, or vector of new integer IDs.
//' @param na Character. Numeric value to use for missing values, but passed as a character scalar.
//' @inheritParams convert_phases
//' @return Data.frame with columns \code{famID}, \code{sampID}, and \code{newID}.
//' @references 
//' \itemize{
//'  \item PLINK. Purcell and Chang. \url{https://www.cog-genomics.org/plink2}
//'  \item Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ (2015) Second-generation PLINK: rising to the challenge of larger and richer datasets. GigaScience, 4. doi: \href{https://doi.org/10.1186/s13742-015-0047-8}{10.1186/s13742-015-0047-8} \href{http://gigascience.biomedcentral.com/articles/10.1186/s13742-015-0047-8}{link}.
//' }
//' @export
//' @backref src/convert_plinkA.cpp
//' @seealso 
//' \code{\link{convert_plink}} is a direct conversion that does not rely on PLINK.
// [[Rcpp::export]]
long convert_plinkA(std::string fnraw,
                    std::string fnout,
                    Rcpp::IntegerVector newID = Rcpp::IntegerVector::create(1),
                    int nlines = -1,
                    std::string na = "9",
                    int idwidth = 4,
                    int precision = -1)
{
  long cnt;  // counter
  long skipped;
  int lastid;
  std::string plc;  // placeholder
  std::string g; 
  std::string line;
  std::vector<float> genotypes;
  
  std::ifstream infile(fnraw.c_str());
  std::ofstream outfile(fnout.c_str());
  
  if (!infile.is_open())
    Rcpp::stop("Could not open '%s'.", fnraw);
  if (!outfile.is_open())
    Rcpp::stop("Could not open '%s' for writing.", fnout);
  
  cnt = 0;
  skipped = 0;
  lastid = newID[0];
  
  while(std::getline(infile, line)) {
    cnt++;
    
    if (newID.size() > 1 && IntegerVector::is_na(newID[cnt-1])) {
      ++skipped;
      continue;
    }

    std::istringstream iss(line);
    
    // Gobble up first 6 columns corresponding to famID, sampID, father, mother, sex, and phenotype (not that order).
    for (size_t i = 1; i <= 6; i++) {
      iss >> plc;
    }
    
    outfile.width(idwidth);
    outfile.setf(std::ios::left, std::ios::adjustfield);
    
    if (cnt-1 < newID.size()) {
      outfile << newID[cnt-1] << " ";
      lastid = newID[cnt-1];
    } else {
      outfile << ++lastid << " ";
    }

    while (iss >> g) {
      Rcout << g << " ";
      if (g.compare("NA") == 0) {
        g = na;
      }
      if (precision > -1) {
        outfile.precision(precision);
        outfile << std::fixed;
        outfile << strtof((g).c_str(),0) << " ";
      } else {
        outfile << g << " " ;
      }
    }
    outfile << std::endl;
    
    if (nlines == cnt)
      break; 
  }
  
  infile.close();
  outfile.close();
  
  return cnt - skipped;
}
